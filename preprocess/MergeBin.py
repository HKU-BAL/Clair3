import sys
import logging
import numpy as np
from argparse import ArgumentParser, SUPPRESS
import tables

import clair3.utils as utils

logging.basicConfig(format='%(message)s', level=logging.INFO)

def Run(args):
    in_fn_list = args.in_fn
    out_fn = args.out_fn
    platform = args.platform
    pileup = args.pileup

    global param
    float_type = 'int32'
    if pileup:
        import shared.param_p as param
    else:
        import shared.param_f as param
        float_type = 'int8'

    tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape

    # select all match prefix if file path not exists
    tables.set_blosc_max_threads(64)
    int_atom = tables.Atom.from_dtype(np.dtype(float_type))
    string_atom = tables.StringAtom(itemsize=param.no_of_positions + 50)
    long_string_atom = tables.StringAtom(itemsize=5000)  # max alt_info length
    table_file = tables.open_file(out_fn, mode='w', filters=utils.FILTERS)
    table_file.create_earray(where='/', name='position_matrix', atom=int_atom, shape=[0] + tensor_shape,
                             filters=utils.FILTERS)
    table_file.create_earray(where='/', name='position', atom=string_atom, shape=(0, 1), filters=utils.FILTERS)
    table_file.create_earray(where='/', name='label', atom=int_atom, shape=(0, param.label_size), filters=utils.FILTERS)
    table_file.create_earray(where='/', name='alt_info', atom=long_string_atom, shape=(0, 1), filters=utils.FILTERS)

    table_dict = utils.update_table_dict()
    total_compressed = 0

    for f in in_fn_list:
        print("[INFO] Merging file {}".format(f))
        fi = tables.open_file(f, model='r')
        assert (len(fi.root.label) == len(fi.root.position) == len(fi.root.position_matrix) == len(fi.root.alt_info))
        for index in range(len(fi.root.label)):
            table_dict['label'].append(fi.root.label[index])
            table_dict['position'].append(fi.root.position[index])
            table_dict['position_matrix'].append(fi.root.position_matrix[index])
            table_dict['alt_info'].append(fi.root.alt_info[index])

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

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()


