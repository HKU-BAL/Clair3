"""
BuildPangenomeTensor - Generate pangenome haplotype tensors for FA training bins.

Runs AFTER MergeBin. Reads merged HDF5 bin files, extracts candidate positions,
queries a pangenome wave VCF for haplotype alleles, and produces a companion
HDF5 file with pangenome tensors that correspond 1:1 to the original bin records.

Usage:
    python3 clair3.py BuildPangenomeTensor \
        --bin_fn merged_bin.h5 \
        --out_fn pangenome_tensor.h5 \
        --pangenome_vcf hprc_wave.vcf.gz \
        --ref_fn reference.fa \
        --num_haplotypes 32
"""

import sys
import os
import logging
import numpy as np
from argparse import ArgumentParser
from collections import defaultdict

logging.basicConfig(format='%(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


def Run(args):
    import h5py
    from clair3.utils import ensure_hdf5_plugin_path, _hdf5_compression_kwargs, _append_hdf5_dataset
    from clair3.pangenome import PangenomeDB
    import shared.param_f as param

    ensure_hdf5_plugin_path()

    bin_fn = args.bin_fn
    out_fn = args.out_fn
    pangenome_vcf = args.pangenome_vcf
    ref_fn = args.ref_fn
    num_haplotypes = args.num_haplotypes
    platform = args.platform

    no_of_positions = param.no_of_positions
    flanking_base_num = param.flankingBaseNum
    channel_size = param.channel_size
    if args.enable_dwell_time:
        channel_size += 1

    # --- Step 1: Open input bin file and read positions ---
    logger.info("[INFO] Reading positions from %s", bin_fn)
    bin_file = h5py.File(bin_fn, 'r')
    positions = bin_file['position'][:]  # [N, 1] bytes
    N = len(positions)
    logger.info("[INFO] Total records: %d", N)

    if N == 0:
        logger.info("[INFO] Empty bin file, skipping")
        bin_file.close()
        return

    # --- Step 2: Parse positions -> (chrom, pos_1based) ---
    # Position format: b"chr22:15238453:ACGTACGT...ACGT"
    chrom_pos_list = []
    for i in range(N):
        pos_str = positions[i][0]
        if isinstance(pos_str, bytes):
            pos_str = pos_str.decode('utf-8', 'ignore')
        parts = pos_str.split(':')
        chrom = parts[0]
        pos_1based = int(parts[1])
        chrom_pos_list.append((chrom, pos_1based))

    # --- Step 3: Group by chromosome for efficient batch processing ---
    chrom_groups = defaultdict(list)  # chrom -> [(local_idx, global_idx, pos_1based)]
    for global_idx, (chrom, pos) in enumerate(chrom_pos_list):
        chrom_groups[chrom].append((len(chrom_groups[chrom]), global_idx, pos))

    # --- Step 4: Initialize PangenomeDB ---
    logger.info("[INFO] Initializing PangenomeDB from %s", pangenome_vcf)
    pg_db = PangenomeDB(pangenome_vcf, ref_fasta_path=ref_fn, num_haplotypes=num_haplotypes)

    # --- Step 5: Create output HDF5 ---
    compression_kwargs = _hdf5_compression_kwargs()
    chunk_rows = 500
    pg_shape = (num_haplotypes, no_of_positions, channel_size)

    out_file = h5py.File(out_fn, 'w')
    out_file.create_dataset(
        "pangenome_matrix",
        shape=(N,) + pg_shape,
        dtype=np.int8,
        chunks=(min(chunk_rows, N),) + pg_shape,
        **compression_kwargs,
    )
    # Copy position field for verification of 1:1 correspondence
    out_file.create_dataset(
        "position",
        data=positions,
        dtype=bin_file['position'].dtype,
    )

    # --- Step 6: Build pangenome tensors per chromosome ---
    total_processed = 0
    for chrom in sorted(chrom_groups.keys()):
        entries = chrom_groups[chrom]
        candidate_positions = [pos for _, _, pos in entries]
        global_indices = [gi for _, gi, _ in entries]

        logger.info("[INFO] Processing %s: %d candidates", chrom, len(candidate_positions))

        # Batch build for this chromosome
        pg_tensor = pg_db.build_batch_haplotype_tensor(
            chrom=chrom,
            candidate_positions_1based=candidate_positions,
            no_of_positions=no_of_positions,
            channel_size=channel_size,
            flanking_base_num=flanking_base_num,
        )
        # pg_tensor shape: [len(entries), num_haplotypes, 33, channel_size]

        # Write to output HDF5 at correct global indices
        for local_idx, global_idx in enumerate(global_indices):
            out_file['pangenome_matrix'][global_idx] = pg_tensor[local_idx]

        total_processed += len(entries)
        if total_processed % 50000 == 0 or chrom == sorted(chrom_groups.keys())[-1]:
            logger.info("[INFO] Processed %d/%d records", total_processed, N)

    # --- Cleanup ---
    pg_db.close()
    bin_file.close()
    out_file.close()

    logger.info("[INFO] Pangenome tensor written to %s (%d records, %d haplotypes)",
                out_fn, N, num_haplotypes)


def main():
    parser = ArgumentParser(
        description="Build pangenome haplotype tensors from merged FA bin files")

    parser.add_argument('--bin_fn', type=str, required=True,
                        help="Input merged HDF5 bin file from MergeBin")

    parser.add_argument('--out_fn', type=str, required=True,
                        help="Output pangenome tensor HDF5 file")

    parser.add_argument('--pangenome_vcf', type=str, required=True,
                        help="Path to tabix-indexed pangenome wave VCF (e.g., HPRC)")

    parser.add_argument('--ref_fn', type=str, required=True,
                        help="Path to reference FASTA file")

    parser.add_argument('--num_haplotypes', type=int, default=32,
                        help="Number of pangenome haplotype rows per candidate (default: %(default)d)")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform: ont, hifi, ilmn (default: %(default)s)")

    parser.add_argument('--enable_dwell_time', action='store_true',
                        help="Match dwell-time-enabled tensor channel count")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()
