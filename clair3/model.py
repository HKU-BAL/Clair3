import warnings
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)
    from tensorflow.python.util import deprecation
    deprecation._PRINT_DEPRECATION_WARNINGS = False
    import tensorflow as tf
import logging
import numpy as np
logging.basicConfig(format='%(message)s', level=logging.INFO)
tf.get_logger().setLevel(logging.ERROR)

from clair3.task.main import GT21, GENOTYPE, VARIANT_LENGTH_1, VARIANT_LENGTH_2
import shared.param_f as param
params = dict(
            float_type=tf.float32,
            task_loss_weights=[
                1,                       # gt21
                1,                       # genotype
                1,                       # variant/indel length 0
                1,                       # variant/indel length 1
                1                        # l2 loss
            ],
            output_shape=GT21.output_label_count + \
                         GENOTYPE.output_label_count + \
                         VARIANT_LENGTH_1.output_label_count + \
                         VARIANT_LENGTH_2.output_label_count,
            output_gt21_shape=GT21.output_label_count,
            output_genotype_shape=GENOTYPE.output_label_count,
            output_indel_length_shape_1=VARIANT_LENGTH_1.output_label_count,
            output_indel_length_shape_2=VARIANT_LENGTH_2.output_label_count,
            output_gt21_entropy_weights=[1] * GT21.output_label_count,
            output_genotype_entropy_weights=[1] * GENOTYPE.output_label_count,
            output_indel_length_entropy_weights_1=[1] * VARIANT_LENGTH_1.output_label_count,
            output_indel_length_entropy_weights_2=[1] * VARIANT_LENGTH_2.output_label_count,
            L3_dropout_rate=0.2,
            L4_num_units=256,
            L4_pileup_num_units=128,
            L4_dropout_rate=0.5,
            L5_1_num_units=128,
            L5_1_dropout_rate=0.2,
            L5_2_num_units=128,
            L5_2_dropout_rate=0.2,
            L5_3_num_units=128,
            L5_3_dropout_rate=0.2,
            L5_4_num_units=128,
            L5_4_dropout_rate=0.2,
            LSTM1_num_units=128,
            LSTM2_num_units=160,
            LSTM1_dropout_rate=0,
            LSTM2_dropout_rate=0.5,
            l2_regularization_lambda=param.l2RegularizationLambda,
        )

add_l2_regulation = True
L2_regularizers = tf.keras.regularizers.l2(params['l2_regularization_lambda']) if add_l2_regulation else None

class Clair3_P(tf.keras.Model):
    # Bi-lstm model for clair3 pileup input
    def __init__(self, add_indel_length=False, predict=False):
        super(Clair3_P, self).__init__()

        # output
        self.output_gt21_shape = params['output_gt21_shape']
        self.output_genotype_shape = params['output_genotype_shape']
        self.output_indel_length_shape_1 = params['output_indel_length_shape_1']
        self.output_indel_length_shape_2 = params['output_indel_length_shape_2']

        self.L3_dropout_rate = params['L3_dropout_rate']
        self.L4_num_units = params['L4_num_units']
        self.L4_pileup_num_units = params['L4_pileup_num_units']
        self.L4_dropout_rate = params['L4_dropout_rate']
        self.L5_1_num_units = params['L5_1_num_units']
        self.L5_1_dropout_rate = params['L5_1_dropout_rate']
        self.L5_2_num_units = params['L5_2_num_units']
        self.L5_2_dropout_rate = params['L5_2_dropout_rate']
        self.L5_3_num_units = params['L5_3_num_units']
        self.L5_3_dropout_rate = params['L5_3_dropout_rate']
        self.L5_4_num_units = params['L5_4_num_units']
        self.L5_4_dropout_rate = params['L5_4_dropout_rate']
        self.LSTM1_num_units = params['LSTM1_num_units']
        self.LSTM2_num_units = params['LSTM2_num_units']
        self.LSTM1_dropout_rate = params['LSTM1_dropout_rate']
        self.LSTM2_dropout_rate = params['LSTM2_dropout_rate']

        self.output_label_split = [
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2
        ]

        self.add_indel_length = add_indel_length
        self.predict = predict

        self.LSTM1 = tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(
            units=self.LSTM1_num_units,
            return_sequences=True,
            kernel_regularizer=L2_regularizers
        ))

        self.LSTM2 = tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(
            units=self.LSTM2_num_units,
            return_sequences=True,
            kernel_regularizer=L2_regularizers
        ))

        self.L3_dropout = tf.keras.layers.Dropout(rate=self.L3_dropout_rate)

        self.L3_dropout_flatten = tf.keras.layers.Flatten()

        self.L4 = tf.keras.layers.Dense(units=self.L4_pileup_num_units, activation='selu',kernel_regularizer=L2_regularizers)

        self.L4_dropout = tf.keras.layers.Dropout(rate=self.LSTM2_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_1 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)

        self.L5_1_dropout = tf.keras.layers.Dropout(rate=self.L5_1_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_2 = tf.keras.layers.Dense(units=self.L5_2_num_units, activation='selu', kernel_regularizer=L2_regularizers)

        self.L5_2_dropout = tf.keras.layers.Dropout(rate=self.L5_2_dropout_rate, seed=param.OPERATION_SEED)

        self.Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)

        self.Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)

        if self.add_indel_length:

            self.L5_3 = tf.keras.layers.Dense(units=self.L5_3_num_units, activation='selu', kernel_regularizer=L2_regularizers)

            self.L5_3_dropout = tf.keras.layers.Dropout(rate=self.L5_3_dropout_rate, seed=param.OPERATION_SEED)

            self.L5_4 = tf.keras.layers.Dense(units=self.L5_4_num_units, activation='selu', kernel_regularizer=L2_regularizers)

            self.L5_4_dropout = tf.keras.layers.Dropout(rate=self.L5_4_dropout_rate, seed=param.OPERATION_SEED)

            self.Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu', kernel_regularizer=L2_regularizers)

            self.Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu', kernel_regularizer=L2_regularizers)

        self.softmax = tf.keras.layers.Softmax()


    def call(self, x,):

        x = tf.cast(x, tf.float32)

        x = self.LSTM1(x)  # (batch_size, inp_seq_len, d_model)

        x = self.LSTM2(x)

        x = self.L3_dropout(x)

        x = self.L3_dropout_flatten(x)

        x = self.L4(x)

        x = self.L4_dropout(x)

        l5_1_dropout = self.L5_1_dropout(self.L5_1(x))

        l5_2_dropout = self.L5_2_dropout(self.L5_2(x))

        y_gt21_logits = self.softmax(self.Y_gt21_logits(l5_1_dropout))

        y_genotype_logits = self.softmax(self.Y_genotype_logits(l5_2_dropout))

        if self.add_indel_length:
            l5_3_dropout = self.L5_3_dropout(self.L5_3(x))

            l5_4_dropout = self.L5_4_dropout(self.L5_4(x))

            y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(l5_3_dropout))

            y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(l5_4_dropout))

            if self.predict:
                return tf.concat([y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2], axis=1)

            return [y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2]

        if self.predict:
            return tf.concat([y_gt21_logits, y_genotype_logits],axis=1)

        return [y_gt21_logits, y_genotype_logits]


class BasicConv2D(tf.keras.layers.Layer):
    def __init__(self, filters, kernel_size, strides, padding, SeparableConv=False):
        super(BasicConv2D, self).__init__()
        conv = tf.keras.layers.SeparableConv2D if SeparableConv else tf.keras.layers.Conv2D
        self.conv = conv(filters=filters,
                                           kernel_size=kernel_size,
                                           strides=strides,
                                           padding=padding,
                                           kernel_regularizer=L2_regularizers)
        self.bn = tf.keras.layers.BatchNormalization()
        self.relu = tf.keras.layers.ReLU()

    def call(self, inputs):
        output = self.conv(inputs)
        output = self.bn(output)
        output = self.relu(output)

        return output

class BasicBlock(tf.keras.layers.Layer):

    def __init__(self, filter_num, stride=1,SeparableConv=False):
        super(BasicBlock, self).__init__()
        conv = tf.keras.layers.SeparableConv2D if SeparableConv else tf.keras.layers.Conv2D

        self.conv1 = conv(filters=filter_num,
                                            kernel_size=(3, 3),
                                            strides=stride,
                                            padding="same",
                                            kernel_regularizer=L2_regularizers)
        self.bn1 = tf.keras.layers.BatchNormalization()
        self.conv2 = conv(filters=filter_num,
                                            kernel_size=(3, 3),
                                            strides=1,
                                            padding="same",
                                            kernel_regularizer=L2_regularizers)
        self.bn2 = tf.keras.layers.BatchNormalization()
        if stride != 1:
            self.downsample = tf.keras.Sequential()
            self.downsample.add(tf.keras.layers.Conv2D(filters=filter_num,
                                                       kernel_size=(1, 1),
                                                       strides=stride,
                                                       kernel_regularizer=L2_regularizers))
            self.downsample.add(tf.keras.layers.BatchNormalization())
        else:
            self.downsample = lambda x: x

    def call(self, inputs):
        residual = self.downsample(inputs)

        x = self.conv1(inputs)
        x = self.bn1(x, )
        x = tf.nn.relu(x)
        x = self.conv2(x)
        x = self.bn2(x, )

        output = tf.nn.relu(tf.keras.layers.add([residual, x]))

        return output

def make_basic_block_layer(filter_num, blocks, stride=1, SeparableConv=False):

    res_block = tf.keras.Sequential()

    res_block.add(BasicBlock(filter_num, stride=stride, SeparableConv=SeparableConv))

    for _ in range(1, blocks):
        res_block.add(BasicBlock(filter_num, stride=1,SeparableConv=SeparableConv))

    return res_block

class PyramidPolling(tf.keras.layers.Layer):
    def __init__(self, spatial_pool_size=(3, 2, 1)):
        super(PyramidPolling, self).__init__()

        self.spatial_pool_size = spatial_pool_size
        self.pool_len = len(self.spatial_pool_size)
        self.window_h = np.empty(self.pool_len, dtype=int)
        self.stride_h = np.empty(self.pool_len, dtype=int)
        self.window_w = np.empty(self.pool_len, dtype=int)
        self.stride_w = np.empty(self.pool_len, dtype=int)

        self.flatten = tf.keras.layers.Flatten()

    def build(self, input_shape):
        height = int(input_shape[1])
        width = int(input_shape[2])

        for i in range(self.pool_len):
            self.window_h[i] = self.stride_h[i] = int(np.ceil(height / self.spatial_pool_size[i]))
            self.window_w[i] = self.stride_w[i] = int(np.ceil(width / self.spatial_pool_size[i]))

    def call(self, x):
        for i in range(self.pool_len):
            max_pool = tf.nn.max_pool(x,
                                      ksize=[1, self.window_h[i], self.window_w[i], 1],
                                      strides=[1, self.stride_h[i], self.stride_w[i], 1],
                                      padding='SAME')
            if i == 0:
                pp = self.flatten(max_pool)

            else:
                pp = tf.concat([pp, self.flatten(max_pool)], axis=-1)

        return pp

class Clair3_F(tf.keras.Model):
    # Residual CNN model for clair3 full alignment input
    def __init__(self, add_indel_length=False, predict=False):
        super(Clair3_F, self).__init__()
        self.output_gt21_shape = params['output_gt21_shape']
        self.output_genotype_shape = params['output_genotype_shape']
        self.output_indel_length_shape_1 = params['output_indel_length_shape_1']
        self.output_indel_length_shape_2 = params['output_indel_length_shape_2']

        self.L3_dropout_rate = params['L3_dropout_rate']
        self.L4_num_units = params['L4_num_units']
        self.L4_dropout_rate = params['L4_dropout_rate']
        self.L5_1_num_units = params['L5_1_num_units']
        self.L5_1_dropout_rate = params['L5_1_dropout_rate']
        self.L5_2_num_units = params['L5_2_num_units']
        self.L5_2_dropout_rate = params['L5_2_dropout_rate']
        self.L5_3_num_units = params['L5_3_num_units']
        self.L5_3_dropout_rate = params['L5_3_dropout_rate']
        self.L5_4_num_units = params['L5_4_num_units']
        self.L5_4_dropout_rate = params['L5_4_dropout_rate']

        self.output_label_split = [
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2
        ]

        self.add_indel_length = add_indel_length
        self.predict = predict

        self.conv1 = BasicConv2D(filters=64,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same",)

        self.res_block1 = make_basic_block_layer(filter_num=64,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv3 = BasicConv2D(filters=128,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block2 = make_basic_block_layer(filter_num=128,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv5 = BasicConv2D(filters=256,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block3 = make_basic_block_layer(filter_num=256,
                                            blocks=1, stride=1)

        self.pyramidpolling = PyramidPolling()

        self.L3_dropout = tf.keras.layers.Dropout(rate=self.L3_dropout_rate)

        self.flatten = tf.keras.layers.Flatten()

        self.L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)

        self.L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_1 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)

        self.L5_1_dropout = tf.keras.layers.Dropout(rate=self.L5_1_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_2 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)

        self.L5_2_dropout = tf.keras.layers.Dropout(rate=self.L5_2_dropout_rate, seed=param.OPERATION_SEED)

        self.Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)

        self.Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)

        if self.add_indel_length:
            self.L5_3 = tf.keras.layers.Dense(units=self.L5_3_num_units, activation='selu', kernel_regularizer=L2_regularizers)

            self.L5_3_dropout = tf.keras.layers.Dropout(rate=self.L5_3_dropout_rate, seed=param.OPERATION_SEED)

            self.L5_4 = tf.keras.layers.Dense(units=self.L5_4_num_units, activation='selu', kernel_regularizer=L2_regularizers)

            self.L5_4_dropout = tf.keras.layers.Dropout(rate=self.L5_4_dropout_rate, seed=param.OPERATION_SEED)

            self.Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)

            self.Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.softmax = tf.keras.layers.Softmax()

    def call(self, inputs):

        x = tf.cast(inputs, tf.float32) / param.NORMALIZE_NUM

        x = self.conv1(x)
        x = self.res_block1(x)
        x = self.conv3(x)
        x = self.res_block2(x)
        x = self.conv5(x)
        x = self.res_block3(x)
        x = self.pyramidpolling(x)
        x = self.flatten(self.L3_dropout(x))

        x = self.L4(x)
        x = self.L4_dropout(x)

        l5_1_dropout = self.L5_1_dropout(self.L5_1(x))

        l5_2_dropout = self.L5_2_dropout(self.L5_2(x))

        y_gt21_logits = self.softmax(self.Y_gt21_logits(l5_1_dropout))

        y_genotype_logits = self.softmax(self.Y_genotype_logits(l5_2_dropout))

        if self.add_indel_length:

            l5_3_dropout = self.L5_3_dropout(self.L5_3(x))

            l5_4_dropout = self.L5_4_dropout(self.L5_4(x))

            y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(l5_3_dropout))

            y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(l5_4_dropout))

            if self.predict:

                return tf.concat([y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2], axis=1)

            return [y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2]

        if self.predict:

            return tf.concat([y_gt21_logits, y_genotype_logits],axis=1)

        return [y_gt21_logits, y_genotype_logits]