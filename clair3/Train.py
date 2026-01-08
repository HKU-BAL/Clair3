import logging
import os
import random
import sys
from argparse import ArgumentParser, SUPPRESS
from itertools import accumulate

import h5py
import numpy as np
import torch
try:
    from tqdm import tqdm
except Exception:
    tqdm = None
from torch import nn
from torch.utils.data import DataLoader, Dataset

import clair3.model as model_path
from clair3.utils import ensure_hdf5_plugin_path
from shared.utils import str2bool

logging.basicConfig(format='%(message)s', level=logging.INFO)


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


class FocalLoss(nn.Module):
    """Multi-class focal loss using one-hot labels."""

    def __init__(self, label_shape_cum, task, effective_label_num=None, gamma=2):
        super().__init__()
        self.gamma = gamma
        self.cls_weights = None
        if effective_label_num is not None:
            task_label_num = get_label_task(effective_label_num, label_shape_cum, task)
            cls_weights = cal_class_weight(task_label_num, len(task_label_num))
            cls_weights = torch.tensor(cls_weights, dtype=torch.float32).unsqueeze(0)
            self.cls_weights = cls_weights

    def forward(self, y_true, y_pred):
        y_pred = torch.clamp(y_pred, min=1e-9, max=1 - 1e-9)
        cross_entropy = -y_true * torch.log(y_pred)
        weight = ((1 - y_pred) ** self.gamma) * y_true
        focal_loss = cross_entropy * weight
        if self.cls_weights is not None:
            focal_loss = focal_loss * self.cls_weights.to(y_pred.device)
        return focal_loss.sum(dim=-1)


class BinChunkDataset(Dataset):
    def __init__(self, data, chunk_list, tensor_shape, bin_tensor_shape, chunk_size):
        self.data = data
        self.chunk_list = chunk_list
        self.chunk_size = chunk_size
        self.tensor_shape = list(tensor_shape)
        self.bin_tensor_shape = list(bin_tensor_shape)
        self.random_offset = 0

    def __len__(self):
        return len(self.chunk_list)

    def set_random_offset(self, offset):
        self.random_offset = offset

    def __getitem__(self, index):
        bin_id, chunk_id = self.chunk_list[index]
        start_pos = self.random_offset + chunk_id * self.chunk_size
        end_pos = start_pos + self.chunk_size
        position_matrix = self.data[bin_id]["position_matrix"][start_pos:end_pos]
        label = self.data[bin_id]["label"][start_pos:end_pos]
        if self.bin_tensor_shape[-1] != self.tensor_shape[-1]:
            position_matrix = position_matrix[..., :self.tensor_shape[-1]]
        return position_matrix, label


def collate_chunks(batch):
    position_matrix = np.concatenate([item[0] for item in batch], axis=0)
    label = np.concatenate([item[1] for item in batch], axis=0)
    return position_matrix, label


def get_chunk_list(chunk_offset, train_chunk_num, chunks_per_batch=10, training_dataset_percentage=None):
    """
    get chunk list for training and validation data. we will randomly split training and validation dataset,
    all training data is directly acquired from various tensor bin files.

    """
    need_split_validation_data = training_dataset_percentage is not None
    all_shuffle_chunk_list = []
    training_chunk_list, validation_chunk_list = [], []
    for bin_idx, chunk_num in enumerate(chunk_offset):
        current_chunk_list = [(bin_idx, chunk_idx) for chunk_idx in range(chunk_num)]
        all_shuffle_chunk_list += current_chunk_list
        if need_split_validation_data:
            buffer_chunk_num = chunks_per_batch
            if chunk_num < buffer_chunk_num:
                training_chunk_list += [(bin_idx, chunk_idx) for chunk_idx in range(chunk_num)]
                continue

            training_chunk_num = int((chunk_num - buffer_chunk_num) * training_dataset_percentage)
            validation_chunk_num = int(chunk_num - buffer_chunk_num - training_chunk_num)
            if training_chunk_num > 0:
                training_chunk_list += current_chunk_list[:training_chunk_num]
            if validation_chunk_num > 0:
                validation_chunk_list += current_chunk_list[-validation_chunk_num:]

    if need_split_validation_data:
        return np.array(training_chunk_list), np.array(validation_chunk_list)

    return np.array(all_shuffle_chunk_list[:train_chunk_num]), np.array(all_shuffle_chunk_list[train_chunk_num:])


def exist_file_prefix(exclude_training_samples, f):
    for prefix in exclude_training_samples:
        if prefix in f:
            return True
    return False


def _load_checkpoint(model, checkpoint_path, device):
    checkpoint = torch.load(checkpoint_path, map_location=device)
    if isinstance(checkpoint, dict) and "state_dict" in checkpoint:
        state_dict = checkpoint["state_dict"]
    else:
        state_dict = checkpoint
    model.load_state_dict(state_dict)


def _run_epoch(model, loader, optimizer, loss_funcs, label_shapes, device, train=True, desc=None):
    if train:
        model.train()
    else:
        model.eval()

    epoch_loss = 0.0
    batches = 0
    iterator = loader
    if tqdm is not None:
        iterator = tqdm(loader, desc=desc, unit="batch", leave=False)
    for position_matrix, label in iterator:
        position_matrix = torch.from_numpy(position_matrix).to(device)
        label = torch.from_numpy(label).float().to(device)
        if train:
            optimizer.zero_grad()

        with torch.set_grad_enabled(train):
            outputs = model(position_matrix)
            label_chunks = []
            start_idx = 0
            for chunk_size in label_shapes:
                end_idx = start_idx + chunk_size
                label_chunks.append(label[:, start_idx:end_idx])
                start_idx = end_idx
            task_losses = []
            for task_id, output in enumerate(outputs):
                task_loss = loss_funcs[task_id](label_chunks[task_id], output).mean()
                task_losses.append(task_loss)
            loss = sum(task_losses)

            if train:
                loss.backward()
                optimizer.step()

        epoch_loss += loss.item()
        if tqdm is not None and hasattr(iterator, "set_postfix"):
            iterator.set_postfix(loss=f"{loss.item():.4f}")
        batches += 1

    return epoch_loss / max(1, batches)


def train_model(args):
    platform = args.platform
    pileup = args.pileup
    add_indel_length = args.add_indel_length
    exclude_training_samples = args.exclude_training_samples
    exclude_training_samples = set(exclude_training_samples.split(',')) if exclude_training_samples else set()
    add_validation_dataset = args.random_validation or (args.validation_fn is not None)
    validation_fn = args.validation_fn
    ochk_prefix = args.ochk_prefix if args.ochk_prefix is not None else ""

    ensure_hdf5_plugin_path()
    if pileup:
        import shared.param_p as param
    else:
        import shared.param_f as param

    base_tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape
    tensor_shape = list(base_tensor_shape)
    if args.enable_dwell_time:
        tensor_shape[-1] += 1

    if pileup:
        model = model_path.Clair3_P(add_indel_length=add_indel_length, input_channels=tensor_shape[-1])
    else:
        model = model_path.Clair3_F(add_indel_length=add_indel_length, input_channels=tensor_shape[-1])
    bin_tensor_shape = list(tensor_shape)
    if args.reuse_bin:
        bin_tensor_shape[-1] += 1
        logging.info("[INFO] Reusing %d-channel bin tensors -> %d-channel model inputs (dropping the last channel).",
                     bin_tensor_shape[-1], tensor_shape[-1])

    label_shape = param.label_shape
    label_shapes = label_shape[0:4 if add_indel_length else 2]
    label_shape_cum = param.label_shape_cum
    batch_size, chunk_size = param.trainBatchSize, param.chunk_size
    assert batch_size % chunk_size == 0
    chunks_per_batch = batch_size // chunk_size
    random.seed(param.RANDOM_SEED)
    np.random.seed(param.RANDOM_SEED)
    if param.RANDOM_SEED is not None:
        torch.manual_seed(param.RANDOM_SEED)

    learning_rate = args.learning_rate if args.learning_rate else param.initialLearningRate
    max_epoch = args.maxEpoch if args.maxEpoch else param.maxEpoch
    task_num = 4 if add_indel_length else 2
    mini_epochs = args.mini_epochs

    def populate_dataset_table(file_list, file_path):
        chunk_offset = np.zeros(len(file_list), dtype=int)
        table_dataset_list = []
        for bin_idx, bin_file in enumerate(file_list):
            table_dataset = h5py.File(os.path.join(file_path, bin_file), 'r')
            table_dataset_list.append(table_dataset)
            chunk_num = (len(table_dataset["label"]) - batch_size) // chunk_size
            chunk_offset[bin_idx] = chunk_num
        return table_dataset_list, chunk_offset

    bin_list = os.listdir(args.bin_fn)
    bin_list = [f for f in bin_list if '_20_' not in f and not exist_file_prefix(exclude_training_samples, f)]
    logging.info("[INFO] total %d training bin files: %s", len(bin_list), ','.join(bin_list))

    effective_label_num = None
    table_dataset_list, chunk_offset = populate_dataset_table(bin_list, args.bin_fn)

    if validation_fn:
        val_list = os.listdir(validation_fn)
        logging.info("[INFO] total %d validation bin files: %s", len(val_list), ','.join(val_list))
        validate_table_dataset_list, validate_chunk_offset = populate_dataset_table(val_list, args.validation_fn)

        train_chunk_num = int(sum(chunk_offset))
        train_shuffle_chunk_list, _ = get_chunk_list(chunk_offset, train_chunk_num)

        validate_chunk_num = int(sum(validate_chunk_offset))
        validate_shuffle_chunk_list, _ = get_chunk_list(validate_chunk_offset, validate_chunk_num)
        total_chunks = train_chunk_num + validate_chunk_num
    else:
        total_chunks = int(sum(chunk_offset))
        training_dataset_percentage = param.trainingDatasetPercentage if add_validation_dataset else None
        if add_validation_dataset:
            total_batches = total_chunks // chunks_per_batch
            validate_chunk_num = int(max(1.0, np.floor(total_batches * (1 - training_dataset_percentage))) * chunks_per_batch)
            train_chunk_num = int(total_chunks - validate_chunk_num)
        else:
            train_chunk_num = total_chunks
        train_shuffle_chunk_list, validate_shuffle_chunk_list = get_chunk_list(
            chunk_offset, train_chunk_num, chunks_per_batch, training_dataset_percentage
        )
        train_chunk_num = len(train_shuffle_chunk_list)
        validate_chunk_num = len(validate_shuffle_chunk_list)

    total_steps = max_epoch * (train_chunk_num // chunks_per_batch)

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=learning_rate, weight_decay=param.l2RegularizationLambda
    )

    loss_funcs = [FocalLoss(label_shape_cum, task, effective_label_num) for task in range(task_num)]

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)

    if args.chkpnt_fn is not None:
        _load_checkpoint(model, args.chkpnt_fn, device)
        logging.info("[INFO] Starting from model %s", args.chkpnt_fn)

    logging.info("[INFO] The size of dataset: %d", total_chunks * chunk_size)
    logging.info("[INFO] The training batch size: %d", batch_size)
    logging.info("[INFO] The training learning_rate: %f", learning_rate)
    logging.info("[INFO] Total training steps: %d", total_steps)
    logging.info("[INFO] Maximum training epoch: %d", max_epoch)
    logging.info("[INFO] Mini-epochs per epoch: %d", mini_epochs)
    logging.info("[INFO] Start training...")

    log_path = "training.log"
    with open(log_path, "w") as log_f:
        log_f.write("epoch\tloss\tval_loss\n")

    best_val_loss = None
    patience = 10 * mini_epochs
    patience_count = 0

    total_epochs = max_epoch * mini_epochs
    chunks_per_epoch = max(1, len(train_shuffle_chunk_list) // mini_epochs)

    for epoch in range(total_epochs):
        logging.info("[INFO] Epoch %d/%d", epoch + 1, total_epochs)
        if epoch % mini_epochs == 0:
            np.random.shuffle(train_shuffle_chunk_list)
            random_offset = np.random.randint(0, chunk_size)
        else:
            random_offset = np.random.randint(0, chunk_size)

        start_idx = (epoch % mini_epochs) * chunks_per_epoch
        end_idx = start_idx + chunks_per_epoch
        epoch_chunk_list = train_shuffle_chunk_list[start_idx:end_idx]

        train_dataset = BinChunkDataset(table_dataset_list, epoch_chunk_list, tensor_shape, bin_tensor_shape, chunk_size)
        train_dataset.set_random_offset(random_offset)
        train_loader = DataLoader(
            train_dataset,
            batch_size=chunks_per_batch,
            shuffle=True,
            drop_last=True,
            collate_fn=collate_chunks,
            pin_memory=True,
            num_workers=8,
            prefetch_factor=2,
        )

        train_loss = _run_epoch(
            model,
            train_loader,
            optimizer,
            loss_funcs,
            label_shapes,
            device,
            train=True,
            desc="train",
        )

        val_loss = None
        if add_validation_dataset:
            val_data = validate_table_dataset_list if validation_fn else table_dataset_list
            val_dataset = BinChunkDataset(val_data, validate_shuffle_chunk_list, tensor_shape, bin_tensor_shape, chunk_size)
            val_dataset.set_random_offset(0)
            val_loader = DataLoader(
                val_dataset,
                batch_size=chunks_per_batch,
                shuffle=False,
                drop_last=True,
                collate_fn=collate_chunks,
                pin_memory=True,
                num_workers=8,
                prefetch_factor=2,
            )
            val_loss = _run_epoch(
                model,
                val_loader,
                optimizer,
                loss_funcs,
                label_shapes,
                device,
                train=False,
                desc="val",
            )

        checkpoint_path = f"{ochk_prefix}.{epoch + 1:02d}.pt" if ochk_prefix else f"epoch.{epoch + 1:02d}.pt"
        torch.save(model.state_dict(), checkpoint_path)

        val_loss_str = "" if val_loss is None else f"{val_loss:.6f}"
        with open(log_path, "a") as log_f:
            log_f.write(f"{epoch + 1}\t{train_loss:.6f}\t{val_loss_str}\n")

        if val_loss is not None:
            if best_val_loss is None or val_loss < best_val_loss:
                best_val_loss = val_loss
                patience_count = 0
                torch.save(model.state_dict(), "best_val_loss.pt")
            else:
                patience_count += 1
                if patience_count >= patience:
                    logging.info("[INFO] Early stopping at epoch %d", epoch + 1)
                    break

    for table_dataset in table_dataset_list:
        table_dataset.close()

    if validation_fn:
        for table_dataset in validate_table_dataset_list:
            table_dataset.close()

    if best_val_loss is not None:
        logging.info("[INFO] Best validation loss: %.6f", best_val_loss)


def main():
    parser = ArgumentParser(description="Train a Clair3 model")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--bin_fn', type=str, default="", required=True,
                        help="Binary tensor input generated by Tensor2Bin.py, supports HDF5 bin files")

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

    parser.add_argument('--mini_epochs', type=int, default=1,
                        help="Number of mini-epochs per epoch")

    # Internal process control
    ## In pileup training mode or not
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    ## Add indel length for training and calling, default true for full alignment
    parser.add_argument('--add_indel_length', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--enable_dwell_time', action='store_true',
                        help="Enable dwell time for training, default: %(default)s")
    parser.add_argument('--reuse_bin', action='store_true',
                        help="Reuse 9-channel bin tensors by dropping the last channel before training, default: %(default)s")

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
