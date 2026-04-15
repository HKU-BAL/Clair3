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

    # Detect phasing datasets from the first input file
    has_phasing = False
    if len(in_fn_list) > 0:
        with h5py.File(in_fn_list[0], mode='r') as probe:
            has_phasing = "phasing_matrix" in probe

    if has_phasing:
        # Probe shapes from first file
        with h5py.File(in_fn_list[0], mode='r') as probe:
            pm_shape = probe["phasing_matrix"].shape[1:]   # (matrix_depth, MAX_NEARBY_HETE_SNPS)
            hp_shape = probe["hp_labels"].shape[1:]         # (matrix_depth,)
        table_file.create_dataset(
            "phasing_matrix",
            shape=(0,) + pm_shape, maxshape=(None,) + pm_shape,
            chunks=(chunk_rows,) + pm_shape, dtype=np.int8, **compression_kwargs)
        table_file.create_dataset(
            "hp_labels",
            shape=(0,) + hp_shape, maxshape=(None,) + hp_shape,
            chunks=(chunk_rows,) + hp_shape, dtype=np.int8, **compression_kwargs)
        table_file.create_dataset(
            "phasing_num_variants",
            shape=(0,), maxshape=(None,),
            chunks=(chunk_rows,), dtype=np.int16, **compression_kwargs)
        print("[INFO] Phasing datasets detected, will merge phasing_matrix/hp_labels/phasing_num_variants")

    table_dict = utils.update_table_dict()
    if has_phasing:
        table_dict['phasing_matrix'] = []
        table_dict['hp_labels'] = []
        table_dict['phasing_num_variants'] = []
    total_compressed = 0

    def flush_to_hdf5():
        nonlocal table_dict
        table_dict = utils.write_table_file(table_file, table_dict, tensor_shape, param.label_size, float_type)
        if has_phasing:
            for key in ('phasing_matrix', 'hp_labels', 'phasing_num_variants'):
                if len(table_dict.get(key, [])) == 0 and key in table_file:
                    pass  # already flushed by write_table_file or handled below
            # Flush phasing buffers that write_table_file doesn't know about
            for key in ('phasing_matrix', 'hp_labels', 'phasing_num_variants'):
                buf = table_dict.get(key, [])
                if len(buf) > 0:
                    arr = np.array(buf)
                    ds = table_file[key]
                    old_len = ds.shape[0]
                    ds.resize(old_len + len(arr), axis=0)
                    ds[old_len:] = arr
                table_dict[key] = []

    for f in in_fn_list:
        print("[INFO] Merging file {}".format(f))
        fi = h5py.File(f, mode='r')
        assert (len(fi["label"]) == len(fi["position"]) == len(fi["position_matrix"]) == len(fi["alt_info"]))
        file_has_phasing = has_phasing and "phasing_matrix" in fi
        if has_phasing and not file_has_phasing:
            print("[WARN] File {} missing phasing datasets, skipping phasing for this file".format(f))
        for index in range(len(fi["label"])):
            table_dict['label'].append(fi["label"][index])
            table_dict['position'].append(fi["position"][index])
            table_dict['position_matrix'].append(fi["position_matrix"][index])
            table_dict['alt_info'].append(fi["alt_info"][index])
            if file_has_phasing:
                table_dict['phasing_matrix'].append(fi["phasing_matrix"][index])
                table_dict['hp_labels'].append(fi["hp_labels"][index])
                table_dict['phasing_num_variants'].append(fi["phasing_num_variants"][index])

            total_compressed += 1

            if total_compressed % 500 == 0 and total_compressed > 0:
                flush_to_hdf5()

            if total_compressed % 50000 == 0:
                print("[INFO] Compressed %d tensor" % (total_compressed), file=sys.stderr)
        fi.close()

    if total_compressed % 500 != 0 and total_compressed > 0:
        flush_to_hdf5()
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
