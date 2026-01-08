import argparse
import os
import sys

import numpy as np
import torch

import shared.param_p as param_p
import shared.param_f as param_f
from clair3.model import Clair3_P, Clair3_F, SeparableConv2d


def _require_tensorflow():
    try:
        import tensorflow as tf
    except Exception as exc:
        sys.exit("[ERROR] TensorFlow is required for conversion: {}".format(exc))
    return tf


def _resolve_checkpoint(tf, path):
    if os.path.isdir(path):
        ckpt = tf.train.latest_checkpoint(path)
        if ckpt is None:
            sys.exit("[ERROR] No checkpoint found under {}".format(path))
        return ckpt
    return path


def _copy_dense(tf_layer, torch_layer):
    weights = tf_layer.get_weights()
    if not weights:
        return
    kernel, bias = weights
    torch_layer.weight.data.copy_(torch.from_numpy(kernel.T))
    torch_layer.bias.data.copy_(torch.from_numpy(bias))


def _copy_batchnorm(tf_layer, torch_layer):
    weights = tf_layer.get_weights()
    if len(weights) != 4:
        return
    gamma, beta, moving_mean, moving_var = weights
    torch_layer.weight.data.copy_(torch.from_numpy(gamma))
    torch_layer.bias.data.copy_(torch.from_numpy(beta))
    torch_layer.running_mean.data.copy_(torch.from_numpy(moving_mean))
    torch_layer.running_var.data.copy_(torch.from_numpy(moving_var))


def _copy_conv(tf_layer, torch_layer):
    if isinstance(torch_layer, SeparableConv2d):
        depthwise_kernel, pointwise_kernel, bias = tf_layer.get_weights()
        # depthwise: (kh, kw, in_c, depth_mult) -> (in_c*depth_mult, 1, kh, kw)
        kh, kw, in_c, depth_mult = depthwise_kernel.shape
        depthwise_kernel = depthwise_kernel.transpose(2, 3, 0, 1)
        depthwise_kernel = depthwise_kernel.reshape(in_c * depth_mult, 1, kh, kw)
        torch_layer.depthwise.weight.data.copy_(torch.from_numpy(depthwise_kernel))

        # pointwise: (1, 1, in_c*depth_mult, out_c) -> (out_c, in_c*depth_mult, 1, 1)
        pointwise_kernel = pointwise_kernel.transpose(3, 2, 0, 1)
        torch_layer.pointwise.weight.data.copy_(torch.from_numpy(pointwise_kernel))
        if torch_layer.pointwise.bias is not None:
            torch_layer.pointwise.bias.data.copy_(torch.from_numpy(bias))
        return

    weights = tf_layer.get_weights()
    if len(weights) == 2:
        kernel, bias = weights
    else:
        kernel = weights[0]
        bias = None
    kernel = kernel.transpose(3, 2, 0, 1)
    torch_layer.weight.data.copy_(torch.from_numpy(kernel))
    if bias is not None and torch_layer.bias is not None:
        torch_layer.bias.data.copy_(torch.from_numpy(bias))


def _copy_basic_conv(tf_layer, torch_layer):
    _copy_conv(tf_layer.conv, torch_layer.conv)
    _copy_batchnorm(tf_layer.bn, torch_layer.bn)


def _copy_basic_block(tf_block, torch_block):
    _copy_conv(tf_block.conv1, torch_block.conv1)
    _copy_batchnorm(tf_block.bn1, torch_block.bn1)
    _copy_conv(tf_block.conv2, torch_block.conv2)
    _copy_batchnorm(tf_block.bn2, torch_block.bn2)
    if hasattr(tf_block, "downsample") and hasattr(torch_block, "downsample"):
        if hasattr(tf_block.downsample, "layers") and isinstance(torch_block.downsample, torch.nn.Sequential):
            if len(tf_block.downsample.layers) >= 2:
                _copy_conv(tf_block.downsample.layers[0], torch_block.downsample[0])
                _copy_batchnorm(tf_block.downsample.layers[1], torch_block.downsample[1])


def _copy_lstm(tf_layer, torch_lstm, reverse=False):
    kernel, recurrent_kernel, bias = tf_layer.get_weights()
    suffix = "_reverse" if reverse else ""
    weight_ih = getattr(torch_lstm, "weight_ih_l0" + suffix)
    weight_hh = getattr(torch_lstm, "weight_hh_l0" + suffix)
    bias_ih = getattr(torch_lstm, "bias_ih_l0" + suffix)
    bias_hh = getattr(torch_lstm, "bias_hh_l0" + suffix)

    weight_ih.data.copy_(torch.from_numpy(kernel.T))
    weight_hh.data.copy_(torch.from_numpy(recurrent_kernel.T))
    bias_ih.data.copy_(torch.from_numpy(bias))
    bias_hh.data.zero_()


def _copy_bidirectional_lstm(tf_bi, torch_lstm):
    _copy_lstm(tf_bi.forward_layer, torch_lstm, reverse=False)
    _copy_lstm(tf_bi.backward_layer, torch_lstm, reverse=True)


def _convert_pileup(tf_model, torch_model, add_indel_length):
    _copy_bidirectional_lstm(tf_model.LSTM1, torch_model.LSTM1)
    _copy_bidirectional_lstm(tf_model.LSTM2, torch_model.LSTM2)
    _copy_dense(tf_model.L4, torch_model.L4)
    _copy_dense(tf_model.L5_1, torch_model.L5_1)
    _copy_dense(tf_model.L5_2, torch_model.L5_2)
    _copy_dense(tf_model.Y_gt21_logits, torch_model.Y_gt21_logits)
    _copy_dense(tf_model.Y_genotype_logits, torch_model.Y_genotype_logits)
    if add_indel_length:
        _copy_dense(tf_model.L5_3, torch_model.L5_3)
        _copy_dense(tf_model.L5_4, torch_model.L5_4)
        _copy_dense(tf_model.Y_indel_length_logits_1, torch_model.Y_indel_length_logits_1)
        _copy_dense(tf_model.Y_indel_length_logits_2, torch_model.Y_indel_length_logits_2)


def _convert_full_alignment(tf_model, torch_model, add_indel_length):
    _copy_basic_conv(tf_model.conv1, torch_model.conv1)
    _copy_basic_block(tf_model.res_block1.layers[0], torch_model.res_block1[0])
    _copy_basic_conv(tf_model.conv3, torch_model.conv3)
    _copy_basic_block(tf_model.res_block2.layers[0], torch_model.res_block2[0])
    _copy_basic_conv(tf_model.conv5, torch_model.conv5)
    _copy_basic_block(tf_model.res_block3.layers[0], torch_model.res_block3[0])
    _copy_dense(tf_model.L4, torch_model.L4)
    _copy_dense(tf_model.L5_1, torch_model.L5_1)
    _copy_dense(tf_model.L5_2, torch_model.L5_2)
    _copy_dense(tf_model.Y_gt21_logits, torch_model.Y_gt21_logits)
    _copy_dense(tf_model.Y_genotype_logits, torch_model.Y_genotype_logits)
    if add_indel_length:
        _copy_dense(tf_model.L5_3, torch_model.L5_3)
        _copy_dense(tf_model.L5_4, torch_model.L5_4)
        _copy_dense(tf_model.Y_indel_length_logits_1, torch_model.Y_indel_length_logits_1)
        _copy_dense(tf_model.Y_indel_length_logits_2, torch_model.Y_indel_length_logits_2)


def _build_tf_models(tf, model_type, add_indel_length, input_channels):
    class BasicConv2D(tf.keras.layers.Layer):
        def __init__(self, filters, kernel_size, strides, padding, separable=False):
            super().__init__()
            conv = tf.keras.layers.SeparableConv2D if separable else tf.keras.layers.Conv2D
            self.conv = conv(
                filters=filters,
                kernel_size=kernel_size,
                strides=strides,
                padding=padding,
                kernel_regularizer=None,
            )
            self.bn = tf.keras.layers.BatchNormalization()
            self.relu = tf.keras.layers.ReLU()

        def call(self, inputs):
            output = self.conv(inputs)
            output = self.bn(output)
            return self.relu(output)

    class BasicBlock(tf.keras.layers.Layer):
        def __init__(self, filter_num, stride=1, separable=False):
            super().__init__()
            conv = tf.keras.layers.SeparableConv2D if separable else tf.keras.layers.Conv2D
            self.conv1 = conv(filters=filter_num, kernel_size=(3, 3), strides=stride, padding="same")
            self.bn1 = tf.keras.layers.BatchNormalization()
            self.conv2 = conv(filters=filter_num, kernel_size=(3, 3), strides=1, padding="same")
            self.bn2 = tf.keras.layers.BatchNormalization()
            if stride != 1:
                self.downsample = tf.keras.Sequential(
                    [
                        tf.keras.layers.Conv2D(filters=filter_num, kernel_size=(1, 1), strides=stride),
                        tf.keras.layers.BatchNormalization(),
                    ]
                )
            else:
                self.downsample = lambda x: x

        def call(self, inputs):
            residual = self.downsample(inputs)
            x = self.conv1(inputs)
            x = self.bn1(x)
            x = tf.nn.relu(x)
            x = self.conv2(x)
            x = self.bn2(x)
            return tf.nn.relu(tf.keras.layers.add([residual, x]))

    class PyramidPolling(tf.keras.layers.Layer):
        def __init__(self, spatial_pool_size=(3, 2, 1)):
            super().__init__()
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
                max_pool = tf.nn.max_pool(
                    x,
                    ksize=[1, self.window_h[i], self.window_w[i], 1],
                    strides=[1, self.stride_h[i], self.stride_w[i], 1],
                    padding='SAME',
                )
                if i == 0:
                    pp = self.flatten(max_pool)
                else:
                    pp = tf.concat([pp, self.flatten(max_pool)], axis=-1)
            return pp

    def make_basic_block_layer(filter_num, blocks, stride=1, separable=False):
        res_block = tf.keras.Sequential()
        res_block.add(BasicBlock(filter_num, stride=stride, separable=separable))
        for _ in range(1, blocks):
            res_block.add(BasicBlock(filter_num, stride=1, separable=separable))
        return res_block

    class Clair3_P(tf.keras.Model):
        def __init__(self, add_indel_length=False, predict=False):
            super().__init__()
            self.output_gt21_shape = param_p.label_shape[0]
            self.output_genotype_shape = param_p.label_shape[1]
            self.output_indel_length_shape_1 = param_p.label_shape[2]
            self.output_indel_length_shape_2 = param_p.label_shape[3]
            self.add_indel_length = add_indel_length
            self.predict = predict

            self.LSTM1 = tf.keras.layers.Bidirectional(
                tf.keras.layers.LSTM(units=128, return_sequences=True)
            )
            self.LSTM2 = tf.keras.layers.Bidirectional(
                tf.keras.layers.LSTM(units=160, return_sequences=True)
            )
            self.L3_dropout = tf.keras.layers.Dropout(rate=0.2)
            self.L3_dropout_flatten = tf.keras.layers.Flatten()
            self.L4 = tf.keras.layers.Dense(units=128, activation='selu')
            self.L4_dropout = tf.keras.layers.Dropout(rate=0.5)
            self.L5_1 = tf.keras.layers.Dense(units=128, activation='selu')
            self.L5_1_dropout = tf.keras.layers.Dropout(rate=0.2)
            self.L5_2 = tf.keras.layers.Dense(units=128, activation='selu')
            self.L5_2_dropout = tf.keras.layers.Dropout(rate=0.2)
            self.Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu')
            self.Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu')
            if self.add_indel_length:
                self.L5_3 = tf.keras.layers.Dense(units=128, activation='selu')
                self.L5_3_dropout = tf.keras.layers.Dropout(rate=0.2)
                self.L5_4 = tf.keras.layers.Dense(units=128, activation='selu')
                self.L5_4_dropout = tf.keras.layers.Dropout(rate=0.2)
                self.Y_indel_length_logits_1 = tf.keras.layers.Dense(
                    units=self.output_indel_length_shape_1, activation='selu'
                )
                self.Y_indel_length_logits_2 = tf.keras.layers.Dense(
                    units=self.output_indel_length_shape_2, activation='selu'
                )
            self.softmax = tf.keras.layers.Softmax()

        def call(self, x):
            x = tf.cast(x, tf.float32)
            x = self.LSTM1(x)
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
                    return tf.concat(
                        [y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2], axis=1
                    )
                return [y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2]
            if self.predict:
                return tf.concat([y_gt21_logits, y_genotype_logits], axis=1)
            return [y_gt21_logits, y_genotype_logits]

    class Clair3_F(tf.keras.Model):
        def __init__(self, add_indel_length=False, predict=False):
            super().__init__()
            self.output_gt21_shape = param_f.label_shape[0]
            self.output_genotype_shape = param_f.label_shape[1]
            self.output_indel_length_shape_1 = param_f.label_shape[2]
            self.output_indel_length_shape_2 = param_f.label_shape[3]
            self.add_indel_length = add_indel_length
            self.predict = predict

            self.conv1 = BasicConv2D(filters=64, kernel_size=(3, 3), strides=2, padding="same")
            self.res_block1 = make_basic_block_layer(filter_num=64, blocks=1, stride=1)
            self.conv3 = BasicConv2D(filters=128, kernel_size=(3, 3), strides=2, padding="same")
            self.res_block2 = make_basic_block_layer(filter_num=128, blocks=1, stride=1)
            self.conv5 = BasicConv2D(filters=256, kernel_size=(3, 3), strides=2, padding="same")
            self.res_block3 = make_basic_block_layer(filter_num=256, blocks=1, stride=1)
            self.pyramidpolling = PyramidPolling()
            self.L3_dropout = tf.keras.layers.Dropout(rate=0.2)
            self.flatten = tf.keras.layers.Flatten()
            self.L4 = tf.keras.layers.Dense(units=256, activation='selu')
            self.L4_dropout = tf.keras.layers.Dropout(rate=0.5)
            self.L5_1 = tf.keras.layers.Dense(units=128, activation='selu')
            self.L5_1_dropout = tf.keras.layers.Dropout(rate=0.2)
            self.L5_2 = tf.keras.layers.Dense(units=128, activation='selu')
            self.L5_2_dropout = tf.keras.layers.Dropout(rate=0.2)
            self.Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu')
            self.Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu')
            if self.add_indel_length:
                self.L5_3 = tf.keras.layers.Dense(units=128, activation='selu')
                self.L5_3_dropout = tf.keras.layers.Dropout(rate=0.2)
                self.L5_4 = tf.keras.layers.Dense(units=128, activation='selu')
                self.L5_4_dropout = tf.keras.layers.Dropout(rate=0.2)
                self.Y_indel_length_logits_1 = tf.keras.layers.Dense(
                    units=self.output_indel_length_shape_1, activation='selu'
                )
                self.Y_indel_length_logits_2 = tf.keras.layers.Dense(
                    units=self.output_indel_length_shape_2, activation='selu'
                )
            self.softmax = tf.keras.layers.Softmax()

        def call(self, inputs):
            x = tf.cast(inputs, tf.float32) / param_f.NORMALIZE_NUM
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
                    return tf.concat(
                        [y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2], axis=1
                    )
                return [y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2]
            if self.predict:
                return tf.concat([y_gt21_logits, y_genotype_logits], axis=1)
            return [y_gt21_logits, y_genotype_logits]

    if model_type == "pileup":
        model = Clair3_P(add_indel_length=add_indel_length, predict=True)
    else:
        model = Clair3_F(add_indel_length=add_indel_length, predict=True)

    return model


def _convert_single(tf, checkpoint, output_path, model_type, platform, add_indel_length, enable_dwell_time):
    if model_type == "pileup":
        tensor_shape = param_p.input_shape
        input_channels = tensor_shape[-1]
        dummy = tf.zeros((1, tensor_shape[0], input_channels), dtype=tf.float32)
        tf_model = _build_tf_models(tf, "pileup", add_indel_length, input_channels)
        tf_model(dummy)
        try:
            tf_model.load_weights(checkpoint).expect_partial()
        except Exception:
            ckpt = tf.train.Checkpoint(model=tf_model)
            ckpt.restore(checkpoint).expect_partial()

        torch_model = Clair3_P(add_indel_length=add_indel_length, predict=True, input_channels=input_channels)
        _convert_pileup(tf_model, torch_model, add_indel_length)
    else:
        tensor_shape = param_f.ont_input_shape if platform == "ont" else param_f.input_shape
        input_channels = tensor_shape[-1]
        if enable_dwell_time:
            input_channels += 1
        dummy = tf.zeros((1, tensor_shape[0], tensor_shape[1], input_channels), dtype=tf.float32)
        tf_model = _build_tf_models(tf, "full_alignment", add_indel_length, input_channels)
        tf_model(dummy)
        try:
            tf_model.load_weights(checkpoint).expect_partial()
        except Exception:
            ckpt = tf.train.Checkpoint(model=tf_model)
            ckpt.restore(checkpoint).expect_partial()

        torch_model = Clair3_F(add_indel_length=add_indel_length, predict=True, input_channels=input_channels)
        _convert_full_alignment(tf_model, torch_model, add_indel_length)

    torch.save(torch_model.state_dict(), output_path)
    print("[INFO] Saved PyTorch checkpoint to {}".format(output_path))


def main():
    parser = argparse.ArgumentParser(description="Convert Clair3 TensorFlow checkpoints to PyTorch state_dict")
    parser.add_argument("--checkpoint", help="TF checkpoint prefix or directory")
    parser.add_argument("--output", help="Output .pt file path")
    parser.add_argument(
        "--model-type",
        choices=["pileup", "full_alignment"],
        help="Model type to convert",
    )
    parser.add_argument("--checkpoint-dir", help="Directory containing pileup/full_alignment checkpoints")
    parser.add_argument("--output-dir", help="Output directory for both pileup and full-alignment .pt files")
    parser.add_argument("--pileup-prefix", default="pileup", help="Pileup checkpoint prefix name")
    parser.add_argument("--fa-prefix", default="full_alignment", help="Full-alignment checkpoint prefix name")
    parser.add_argument("--platform", default="ont", choices=["ont", "hifi", "ilmn"], help="Platform")
    parser.add_argument("--add-indel-length", action="store_true", help="Include indel length heads")
    parser.add_argument("--pileup-add-indel-length", action="store_true", help="Include indel length heads for pileup")
    parser.add_argument("--fa-add-indel-length", action="store_true", help="Include indel length heads for full-alignment")
    parser.add_argument("--enable-dwell-time", action="store_true", help="Enable dwell time channel (full alignment)")

    args = parser.parse_args()

    if args.output_dir:
        if args.checkpoint_dir:
            ckpt_dir = args.checkpoint_dir
        elif args.checkpoint and os.path.isdir(args.checkpoint):
            ckpt_dir = args.checkpoint
        else:
            sys.exit("[ERROR] --output-dir requires --checkpoint-dir or a directory passed via --checkpoint")

        os.makedirs(args.output_dir, exist_ok=True)
        tf = _require_tensorflow()

        pileup_ckpt = _resolve_checkpoint(tf, os.path.join(ckpt_dir, args.pileup_prefix))
        fa_ckpt = _resolve_checkpoint(tf, os.path.join(ckpt_dir, args.fa_prefix))

        pileup_add_indel = args.pileup_add_indel_length
        fa_add_indel = args.fa_add_indel_length or args.add_indel_length

        _convert_single(
            tf,
            pileup_ckpt,
            os.path.join(args.output_dir, f"{args.pileup_prefix}.pt"),
            "pileup",
            args.platform,
            pileup_add_indel,
            False,
        )
        _convert_single(
            tf,
            fa_ckpt,
            os.path.join(args.output_dir, f"{args.fa_prefix}.pt"),
            "full_alignment",
            args.platform,
            fa_add_indel,
            args.enable_dwell_time,
        )
        return

    if not args.checkpoint or not args.output or not args.model_type:
        sys.exit("[ERROR] --checkpoint, --output, and --model-type are required unless --output-dir is set")

    tf = _require_tensorflow()
    checkpoint = _resolve_checkpoint(tf, args.checkpoint)
    add_indel_length = args.add_indel_length
    if args.model_type == "pileup" and args.pileup_add_indel_length:
        add_indel_length = True
    if args.model_type == "full_alignment" and args.fa_add_indel_length:
        add_indel_length = True

    _convert_single(
        tf,
        checkpoint,
        args.output,
        args.model_type,
        args.platform,
        add_indel_length,
        args.enable_dwell_time,
    )


if __name__ == "__main__":
    main()
