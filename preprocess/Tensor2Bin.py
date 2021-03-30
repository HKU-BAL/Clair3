import sys
import logging
from argparse import ArgumentParser

import clair3.utils as utils

logging.basicConfig(format='%(message)s', level=logging.INFO)

def Run(args):

    
    utils.setup_environment()
    logging.info("Loading the dataset ...")
    
    utils.get_training_array(
        tensor_fn=args.tensor_fn,
        var_fn=args.var_fn,
        bed_fn=args.bed_fn,
        bin_fn=args.bin_fn,
        shuffle=args.shuffle,
        is_allow_duplicate_chr_pos=args.allow_duplicate_chr_pos,
        chunk_id=args.chunk_id-1, # 1-base to 0-base
        chunk_num=args.chunk_num,
        pileup=args.pileup,
        platform=args.platform,
        maximum_non_variant_ratio=args.maximum_non_variant_ratio,
        alt_fn_prefix=args.alt_fn_prefix)
    logging.info("Finish!")


def main():
    parser = ArgumentParser(description="Generate a binary format input tensor")

    parser.add_argument('--tensor_fn', type=str, default="vartensors",
                        help="Tensor input")

    parser.add_argument('--var_fn', type=str, default="truthvars",
                        help="Truth variants list input")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="High confident genome regions input in the BED format")

    parser.add_argument('--bin_fn', type=str, default=None,
                        help="Output a binary tensor file")

    parser.add_argument('--shuffle', action='store_true',
                        help="Shuffle on building bin")

    parser.add_argument('--allow_duplicate_chr_pos', action='store_true',
                        help="Allow duplicate chromosome:position in tensor input")

    parser.add_argument('--calling', action='store_true',
                        help="will create tf record if in calling process")

    parser.add_argument('--chunk_id', type=int, default=None,
                        help="Specific chunk id works with total chunk_num for parallel execution.")

    parser.add_argument('--chunk_num', type=int, default=None,
                        help="Total chunk number for parallel execution. Each chunk refer to a smaller reference regions.")

    parser.add_argument('--platform', type=str, default='pb',
                        help="Select specific platform for variant calling. Optional: 'ont,pb,illumina', default: %(default)s")

    parser.add_argument('--pileup', action='store_true',
                        help="Whether in pileup mode. Define two bin generation mode, pileup or full alignment, default: False")

    # for training data subsample
    parser.add_argument('--maximum_non_variant_ratio', type=float, default=None,
                        help="subsample ratio for non-variant training data, optional")

    parser.add_argument('--alt_fn_prefix', type=str, default=None,
                        help="Path of alternative variants, works with maximum_non_variant_ratio in training, optional")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()
