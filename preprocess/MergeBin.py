import sys
import logging
import numpy as np
from argparse import ArgumentParser, SUPPRESS
import h5py

import clair3.utils as utils

logging.basicConfig(format='%(message)s', level=logging.INFO)

def Run(args):
    in_fn_list = args.in_fn
    out_fn = args.out_fn
    platform = args.platform
    pileup = args.pileup
    enable_dwell_time = args.enable_dwell_time
    global param
    float_type = 'int32'
    if pileup:
        import shared.param_p as param
    else:
        import shared.param_f as param
        float_type = 'int8'

    tensor_shape = list(param.ont_input_shape if platform == 'ont' else param.input_shape)
    if enable_dwell_time:
        tensor_shape[2] += 1

    # select all match prefix if file path not exists
    utils.ensure_hdf5_plugin_path()
    compression_kwargs = utils._hdf5_compression_kwargs()
    table_file = h5py.File(out_fn, mode='w')
    chunk_rows = 500
    table_file.create_dataset(
        "position_matrix",
        shape=(0,) + tuple(tensor_shape),
        maxshape=(None,) + tuple(tensor_shape),
        chunks=(chunk_rows,) + tuple(tensor_shape),
        dtype=np.dtype(float_type),
        **compression_kwargs,
    )
    table_file.create_dataset(
        "position",
        shape=(0, 1),
        maxshape=(None, 1),
        chunks=(chunk_rows, 1),
        dtype=f"S{param.no_of_positions + 50}",
        **compression_kwargs,
    )
    table_file.create_dataset(
        "label",
        shape=(0, param.label_size),
        maxshape=(None, param.label_size),
        chunks=(chunk_rows, param.label_size),
        dtype=np.dtype(float_type),
        **compression_kwargs,
    )
    table_file.create_dataset(
        "alt_info",
        shape=(0, 1),
        maxshape=(None, 1),
        chunks=(chunk_rows, 1),
        dtype="S5000",
        **compression_kwargs,
    )

    table_dict = utils.update_table_dict()
    total_compressed = 0

    for f in in_fn_list:
        print("[INFO] Merging file {}".format(f))
        fi = h5py.File(f, mode='r')
        assert (len(fi["label"]) == len(fi["position"]) == len(fi["position_matrix"]) == len(fi["alt_info"]))
        for index in range(len(fi["label"])):
            table_dict['label'].append(fi["label"][index])
            table_dict['position'].append(fi["position"][index])
            table_dict['position_matrix'].append(fi["position_matrix"][index])
            table_dict['alt_info'].append(fi["alt_info"][index])

            total_compressed += 1

            if total_compressed % 500 == 0 and total_compressed > 0:
                table_dict = utils.write_table_file(table_file, table_dict, tensor_shape, param.label_size, float_type)

            if total_compressed % 50000 == 0:
                print("[INFO] Compressed %d tensor" % (total_compressed), file=sys.stderr)
        fi.close()

    if total_compressed % 500 != 0 and total_compressed > 0:
        table_dict = utils.write_table_file(table_file, table_dict, tensor_shape, param.label_size, float_type)
        print("[INFO] Compressed %d tensor" % (total_compressed), file=sys.stderr)

    table_file.close()


def main():
    parser = ArgumentParser(description="Combine tensor binaries into a single file")

    parser.add_argument('in_fn', type=str, nargs='+',
                        help="Tensor input files, required")

    parser.add_argument('--out_fn', type=str, default=None, required=True,
                        help="Output a binary tensor file, required")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    ## In pileup mode or not (full alignment mode), default: False
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)
    parser.add_argument('--enable_dwell_time', action='store_true',
                        help="Enable dwell time, default: False")
    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()
