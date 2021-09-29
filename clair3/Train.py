import logging
import random
import numpy as np
from argparse import ArgumentParser, SUPPRESS
import tensorflow_addons as tfa
import tensorflow as tf
import tables
import os
import sys
from itertools import accumulate

import clair3.model as model_path
from shared.utils import str2bool

logging.basicConfig(format='%(message)s', level=logging.INFO)
tables.set_blosc_max_threads(512)
os.environ['NUMEXPR_MAX_THREADS'] = '64'
os.environ['NUMEXPR_NUM_THREADS'] = '8'


def get_label_task(label, label_shape_cum, task):
    if task == 0:
        return label[:label_shape_cum[task]]
    elif task == len(label_shape_cum) - 1:
        return label[label_shape_cum[task - 1]:]
    else:
        return label[label_shape_cum[task - 1]:label_shape_cum[task]]


def cal_class_weight(samples_per_cls, no_of_classes, beta=0.999):
    effective_num = 1.0 - np.power(beta, samples_per_cls)
    cls_weights = (1.0 - beta) / np.array(effective_num)
    cls_weights = cls_weights / np.sum(cls_weights) * no_of_classes
    return cls_weights


class FocalLoss(tf.keras.losses.Loss):
    """
    updated version of focal loss function, for multi class classification, we remove alpha parameter, which the loss
    more stable, and add gradient clipping to avoid gradient explosion and precision overflow.
    """

    def __init__(self, label_shape_cum, task, effective_label_num=None, gamma=2):
        super(FocalLoss, self).__init__()
        self.gamma = gamma
        self.cls_weights = None
        if effective_label_num is not None:
            task_label_num = get_label_task(effective_label_num, label_shape_cum, task)
            cls_weights = cal_class_weight(task_label_num, len(task_label_num))
            cls_weights = tf.constant(cls_weights, dtype=tf.float32)
            cls_weights = tf.expand_dims(cls_weights, axis=0)
            self.cls_weights = cls_weights

    def call(self, y_true, y_pred):
        y_pred = tf.clip_by_value(y_pred, clip_value_min=1e-9, clip_value_max=1 - 1e-9)
        cross_entropy = -y_true * tf.math.log(y_pred)
        weight = ((1 - y_pred) ** self.gamma) * y_true
        FCLoss = cross_entropy * weight
        if self.cls_weights is not None:
            FCLoss = FCLoss * self.cls_weights
        reduce_fl = tf.reduce_sum(FCLoss, axis=-1)
        return reduce_fl


def get_chunk_list(chunk_offset, train_chunk_num):
    """
    get chunk list for training and validation data. we will randomly split training and validation dataset,
    all training data is directly acquired from various tensor bin files.

    """
    all_shuffle_chunk_list = []
    for bin_idx, chunk_num in enumerate(chunk_offset):
        all_shuffle_chunk_list += [(bin_idx, chunk_idx) for chunk_idx in range(chunk_num)]
    np.random.seed(0)
    np.random.shuffle(all_shuffle_chunk_list)  # keep the same random validate dataset
    return np.array(all_shuffle_chunk_list[:train_chunk_num]), np.array(all_shuffle_chunk_list[train_chunk_num:])


def exist_file_prefix(exclude_training_samples, f):
    for prefix in exclude_training_samples:
        if prefix in f:
            return True
    return False


def train_model(args):
    platform = args.platform
    pileup = args.pileup
    add_indel_length = args.add_indel_length
    exclude_training_samples = args.exclude_training_samples
    exclude_training_samples = set(exclude_training_samples.split(',')) if exclude_training_samples else set()
    add_validation_dataset = args.random_validation or (args.validation_fn is not None)
    validation_fn = args.validation_fn
    ochk_prefix = args.ochk_prefix if args.ochk_prefix is not None else ""
    if pileup:
        import shared.param_p as param
        model = model_path.Clair3_P()
    else:
        import shared.param_f as param
        model = model_path.Clair3_F(add_indel_length=add_indel_length)

    tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape
    label_size, label_shape = param.label_size, param.label_shape
    label_shape_cum = list(accumulate(label_shape))
    batch_size, chunk_size = param.trainBatchSize, param.chunk_size
    chunks_per_batch = batch_size // chunk_size
    random.seed(param.RANDOM_SEED)
    np.random.seed(param.RANDOM_SEED)
    learning_rate = args.learning_rate if args.learning_rate else param.initialLearningRate
    max_epoch = args.maxEpoch if args.maxEpoch else param.maxEpoch
    task_num = 4 if add_indel_length else 2
    TensorShape = (
    tf.TensorShape([None] + tensor_shape), tuple(tf.TensorShape([None, label_shape[task]]) for task in range(task_num)))
    TensorDtype = (tf.int32, tuple(tf.float32 for _ in range(task_num)))

    def populate_dataset_table(file_list, file_path):
        chunk_offset = np.zeros(len(file_list), dtype=int)
        table_dataset_list = []
        for bin_idx, bin_file in enumerate(file_list):
            table_dataset = tables.open_file(os.path.join(file_path, bin_file), 'r')
            table_dataset_list.append(table_dataset)
            chunk_num = (len(table_dataset.root.label) - chunk_size) // chunk_size
            chunk_offset[bin_idx] = chunk_num
        return table_dataset_list, chunk_offset

    bin_list = os.listdir(args.bin_fn)
    # default we exclude sample hg003 and all chr20 for training
    bin_list = [f for f in bin_list if '_20_' not in f and not exist_file_prefix(exclude_training_samples, f)]
    logging.info("[INFO] total {} training bin files: {}".format(len(bin_list), ','.join(bin_list)))

    effective_label_num = None

    table_dataset_list, chunk_offset = populate_dataset_table(bin_list, args.bin_fn)

    if validation_fn:
        val_list = os.listdir(validation_fn)
        logging.info("[INFO] total {} validation bin files: {}".format(len(val_list), ','.join(val_list)))
        validate_table_dataset_list, validate_chunk_offset = populate_dataset_table(val_list, args.validation_fn)

        train_chunk_num = int(sum(chunk_offset))
        train_shuffle_chunk_list, _ = get_chunk_list(chunk_offset, train_chunk_num)

        validate_chunk_num = int(sum(validate_chunk_offset))
        validate_shuffle_chunk_list, _ = get_chunk_list(validate_chunk_offset, validate_chunk_num)
        total_chunks = train_chunk_num + validate_chunk_num
    else:
        total_chunks = int(sum(chunk_offset))
        if add_validation_dataset:
            total_batches = total_chunks // chunks_per_batch
            validate_chunk_num = int(max(1., np.floor(total_batches * (1 - param.trainingDatasetPercentage))) * chunks_per_batch)
            train_chunk_num = int(total_chunks - validate_chunk_num)
        else:
            train_chunk_num = total_chunks
        train_shuffle_chunk_list, validate_shuffle_chunk_list = get_chunk_list(chunk_offset, train_chunk_num)
        train_chunk_num = len(train_shuffle_chunk_list)
        validate_chunk_num = len(validate_shuffle_chunk_list)


    def DataGenerator(x, shuffle_chunk_list, train_flag=True):

        """
        data generator for pileup or full alignment data processing, pytables with blosc:lz4hc are used for extreme fast
        compression and decompression. random chunk shuffling and random start position to increase training model robustness.

        """

        batch_num = len(shuffle_chunk_list) // chunks_per_batch
        position_matrix = np.empty([batch_size] + tensor_shape, np.int32)
        label = np.empty((batch_size, param.label_size), np.float32)

        random_start_position = np.random.randint(0, chunk_size) if train_flag else 0
        if train_flag:
            np.random.shuffle(shuffle_chunk_list)
        for batch_idx in range(batch_num):
            for chunk_idx in range(chunks_per_batch):
                bin_id, chunk_id = shuffle_chunk_list[batch_idx * chunks_per_batch + chunk_idx]
                position_matrix[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.position_matrix[
                        random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                label[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.label[
                        random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]

            if add_indel_length:
                yield position_matrix, (
                        label[:,                   :label_shape_cum[0]],
                        label[:, label_shape_cum[0]:label_shape_cum[1]],
                        label[:, label_shape_cum[1]:label_shape_cum[2]],
                        label[:, label_shape_cum[2]:                  ]
                    )
            else:
                yield position_matrix, (
                        label[:,                   :label_shape_cum[0]],
                        label[:, label_shape_cum[0]:label_shape_cum[1]]
                    )


    train_dataset = tf.data.Dataset.from_generator(
        lambda: DataGenerator(table_dataset_list, train_shuffle_chunk_list, True), TensorDtype,
        TensorShape).prefetch(buffer_size=tf.data.experimental.AUTOTUNE)
    validate_dataset = tf.data.Dataset.from_generator(
        lambda: DataGenerator(validate_table_dataset_list if validation_fn else table_dataset_list, validate_shuffle_chunk_list, False), TensorDtype,
        TensorShape).prefetch(buffer_size=tf.data.experimental.AUTOTUNE)

    total_steps = max_epoch * (train_chunk_num // chunks_per_batch)

    #RectifiedAdam with warmup start
    optimizer = tfa.optimizers.Lookahead(tfa.optimizers.RectifiedAdam(
        lr=learning_rate,
        total_steps=total_steps,
        warmup_proportion=0.1,
        min_lr=learning_rate*0.75,
    ))

    loss_func = [FocalLoss(label_shape_cum, task, effective_label_num) for task in range(task_num)]
    loss_task = {"output_{}".format(task + 1): loss_func[task] for task in range(task_num)}
    metrics = {"output_{}".format(task + 1): tfa.metrics.F1Score(num_classes=label_shape[task], average='micro') for
               task in range(task_num)}

    model.compile(
        loss=loss_task,
        metrics=metrics,
        optimizer=optimizer
    )
    early_stop_callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10, mode="min")
    model_save_callback = tf.keras.callbacks.ModelCheckpoint(ochk_prefix + ".{epoch:02d}", period=1, save_weights_only=False)
    model_best_callback = tf.keras.callbacks.ModelCheckpoint("best_val_loss", monitor='val_loss', save_best_only=True, mode="min")
    train_log_callback = tf.keras.callbacks.CSVLogger("training.log", separator='\t')

    # Use first 20 element to initialize tensorflow model using graph mode
    output = model(np.array(table_dataset_list[0].root.position_matrix[:20]))
    logging.info(model.summary(print_fn=logging.info))

    logging.info("[INFO] The size of dataset: {}".format(total_chunks * chunk_size))
    logging.info("[INFO] The training batch size: {}".format(batch_size))
    logging.info("[INFO] The training learning_rate: {}".format(learning_rate))
    logging.info("[INFO] Total training steps: {}".format(total_steps))
    logging.info("[INFO] Maximum training epoch: {}".format(max_epoch))
    logging.info("[INFO] Start training...")

    validate_dataset = validate_dataset if add_validation_dataset else None
    if args.chkpnt_fn is not None:
        model.load_weights(args.chkpnt_fn)
        logging.info("[INFO] Starting from model {}".format(args.chkpnt_fn))

    train_history = model.fit(x=train_dataset,
                              epochs=max_epoch,
                              validation_data=validate_dataset,
                              callbacks=[early_stop_callback,
                                         model_save_callback,
                                         model_best_callback,
                                         train_log_callback],
                              verbose=1,
                              shuffle=False)

    for table_dataset in table_dataset_list:
        table_dataset.close()

    if validation_fn:
        for table_dataset in validate_table_dataset_list:
            table_dataset.close()

    # show the parameter set with the smallest validation loss
    if 'val_loss' in train_history.history:
        best_validation_epoch = np.argmin(np.array(train_history.history["val_loss"])) + 1
        logging.info("[INFO] Best validation loss at epoch: %d" % best_validation_epoch)
    else:
        best_train_epoch = np.argmin(np.array(train_history.history["loss"])) + 1
        logging.info("[INFO] Best train loss at epoch: %d" % best_train_epoch)


def main():
    parser = ArgumentParser(description="Train a Clair3 model")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--bin_fn', type=str, default="", required=True,
                        help="Binary tensor input generated by Tensor2Bin.py, support multiple bin readers using pytables")

    parser.add_argument('--chkpnt_fn', type=str, default=None,
                        help="Input a model to resume training or for fine-tuning")

    parser.add_argument('--ochk_prefix', type=str, default=None, required=True,
                        help="Prefix for model output after each epoch")

    # options for advanced users
    parser.add_argument('--maxEpoch', type=int, default=None,
                        help="Maximum number of training epochs")

    parser.add_argument('--learning_rate', type=float, default=1e-3,
                        help="Set the initial learning rate, default: %(default)s")


    parser.add_argument('--exclude_training_samples', type=str, default=None,
                        help="Define training samples to be excluded")

    # Internal process control
    ## In pileup training mode or not
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    ## Add indel length for training and calling, default true for full alignment
    parser.add_argument('--add_indel_length', type=str2bool, default=False,
                        help=SUPPRESS)

    # mutually-incompatible validation options
    vgrp = parser.add_mutually_exclusive_group()
    vgrp.add_argument('--random_validation', action='store_true',
                        help="Use random sample of dataset for validation, default: %(default)s")

    vgrp.add_argument('--validation_fn', type=str, default=None,
                        help="Binary tensor input for use in validation: %(default)s")


    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    train_model(args)


if __name__ == "__main__":
    main()
