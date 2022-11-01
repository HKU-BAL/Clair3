
# install and import the necessary packages
import sys
import logging
import os
import sys
import tables
import numpy as np
from argparse import ArgumentParser, SUPPRESS
import tensorflow as tf
import tensorflow.keras as keras
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm

current_path = os.getcwd()
sys.path.append(os.path.join(current_path, 'Clair3'))
import shared.param_f as param
from clair3.model_GP import Clair3_F

#disable cuda for computation
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
logging.basicConfig(format='%(message)s', level=logging.INFO)

label_list = ['Reference base', 'Alternative base', 'Strand Info', 'Mapping quality', \
              'Base quality', 'Candidate proportion', 'Insertion base', 'Phasing Info']

class FocalLoss(tf.keras.losses.Loss):
    """
    updated version of focal loss function, for multi class classification, we remove alpha parameter, which the loss
    more stable, and add gradient clipping to avoid gradient explosion and precision overflow.
    """

    def __init__(self, label_shape_cum, task, effective_label_num=None, gamma=2):
        super(FocalLoss, self).__init__()
        self.gamma = gamma
        self.cls_weights = None

    def call(self, y_true, y_pred):
        y_pred = tf.clip_by_value(y_pred, clip_value_min=1e-9, clip_value_max=1 - 1e-9)
        cross_entropy = -y_true * tf.math.log(y_pred)
        weight = ((1 - y_pred) ** self.gamma) * y_true
        FCLoss = cross_entropy * weight
        if self.cls_weights is not None:
            FCLoss = FCLoss * self.cls_weights
        reduce_fl = tf.reduce_sum(FCLoss, axis=-1)
        return reduce_fl

@tf.custom_gradient
def guidedRelu(x):
    def grad(dy):
        return tf.cast(dy > 0, "float32") * tf.cast(x > 0, "float32") * dy

    return tf.nn.relu(x), grad



def vis_feature(args, model, input_tensor):
    if len(input_tensor.shape) == 3:
        input_tensor = input_tensor[np.newaxis, :, :, :]
    with tf.GradientTape() as tape:
        inputs = tf.cast(input_tensor, tf.float32)
        tape.watch(inputs)
        outputs = model(inputs)
    all_grads = tape.gradient(outputs, inputs)
    grads = all_grads[0]

    tensor_fa = grads
    nb_cols = 4
    fig = plt.figure(figsize=(36, 40))
    gs = gridspec.GridSpec(2, nb_cols)
    axes1 = [fig.add_subplot(gs[0, col], aspect="equal") for col in range(nb_cols)]
    axes2 = [fig.add_subplot(gs[1, col], aspect="equal") for col in range(nb_cols)]
    axes = [axes1, axes2]
    plt.xticks((0, 16, 32))
    for i in range(2):
        for col, ax in enumerate(axes[i]):
            idx = i * nb_cols + col
            im = ax.pcolormesh(grads[:, :, idx], cmap=cm.RdBu, )
            y = list(range(0, 90, 10)) + [89]
            ax.set_yticks(y)
            ax.set_ylim(ax.get_ylim()[::-1])
            if col > 0:
                ax.yaxis.set_visible(False)
            ax.set_xlabel(label_list[idx], fontsize=36)
            ax.set_xticks((0, 16, 32))
            ax.xaxis.set_tick_params(labelsize=25)
            ax.yaxis.set_tick_params(labelsize=25)

    plt.subplots_adjust(hspace=0.08, wspace=0.08)
    cb = fig.colorbar(im, ax=axes1, pad=0.04, shrink=0.98, aspect=40)
    cb.ax.tick_params(labelsize=25)
    cb1 = fig.colorbar(im, ax=axes2, pad=0.04, shrink=0.98, aspect=40)
    cb1.ax.tick_params(labelsize=25)
    # plt.show()
    plt.savefig("{}/guided_backpropagation_feature_maps.png".format(args.output_dir), dpi=100)

def vis_input(args, input_tensor):
    nb_cols = 4
    fig = plt.figure(figsize=(36, 40))
    gs = gridspec.GridSpec(2, nb_cols)
    axes1 = [fig.add_subplot(gs[0, col], aspect="equal") for col in range(nb_cols)]
    axes2 = [fig.add_subplot(gs[1, col], aspect="equal") for col in range(nb_cols)]
    axes = [axes1, axes2]
    plt.xticks((0, 16, 32))
    for i in range(2):
        for col, ax in enumerate(axes[i]):
            idx = i * nb_cols + col
            im = ax.pcolormesh(input_tensor[:, :, idx], cmap=cm.RdBu, vmin=-100, vmax=100)
            y = list(range(0, 90, 10)) + [89]
            ax.set_yticks(y)
            ax.set_ylim(ax.get_ylim()[::-1])
            if col > 0:
                ax.yaxis.set_visible(False)
            ax.set_xlabel(label_list[idx], fontsize=36)
            ax.set_xticks((0, 16, 32))
            ax.xaxis.set_tick_params(labelsize=25)
            ax.yaxis.set_tick_params(labelsize=25)

    plt.subplots_adjust(hspace=0.08, wspace=0.08)
    cb = fig.colorbar(im, ax=axes1, pad=0.04, shrink=0.98, aspect=40)
    cb.ax.tick_params(labelsize=25)
    cb1 = fig.colorbar(im, ax=axes2, pad=0.04, shrink=0.98, aspect=40)
    cb1.ax.tick_params(labelsize=25)
    plt.savefig("{}/tensor.png".format(args.output_dir), dpi=100)


def Run(args):

    # parameters
    platform = args.platform
    tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape
    input_tensor_shape = tuple([None] + tensor_shape)
    label_shape = param.label_shape
    full_alignment_model_path = os.path.join(args.model_path, 'full_alignment')

    # build Clair3 full-alignment model in graph mode
    m = Clair3_F(add_indel_length=True, predict=True)

    m.build(input_shape=input_tensor_shape)

    m.load_weights(full_alignment_model_path)
    task_num = 4
    m.build(input_shape=input_tensor_shape)
    loss_func = [FocalLoss(label_shape, task, None) for task in range(task_num)]
    loss_task = {"output_{}".format(task + 1): loss_func[task] for task in range(task_num)}
    m.compile(
        loss=loss_task,
        optimizer='adam'
    )
    network = keras.Sequential([
        m
    ])

    network.build(input_shape=input_tensor_shape)

    tmp_matrix = np.ones(shape=[10] + tensor_shape)
    prediction = m.predict_on_batch(tmp_matrix)
    output = m(tmp_matrix)

    layer_dict = [layer for layer in m.layers[:] if hasattr(layer, 'activation')]
    for layer in layer_dict:
        if layer.activation == tf.keras.activations.relu:
            layer.activation = guidedRelu

    bin_file = tables.open_file(args.input_bin)
    tensor = np.array(bin_file.root.position_matrix[0])
    bin_file.close()
    print("[INFO] Visualize guided backpropagation feature maps")
    vis_feature(args, m, tensor)
    print("[INFO] Visualize input tensor")
    vis_input(args, tensor)



def main():
    parser = ArgumentParser(description="Clair3 feature visualization")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--model_path', type=str, default=None,
                        help="Normal tensor input, required")

    parser.add_argument('--input_bin', type=str, default=None,
                        help="Tumor tensor input, required")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Truth variants list input, required")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()
