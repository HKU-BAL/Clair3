"""
CreateTrainingTensorDirect - Single-process training tensor creation using CFFI/C.

Replaces the two-process PyPy+CPython text-pipe architecture with a single CPython
process that calls C code via CFFI to generate tensors, then directly matches labels
and writes HDF5 output.

Supports both full-alignment and pileup modes.

Key improvements over CreateTrainingTensor:
  1. No PyPy process, no samtools subprocess, no text serialization
  2. C code (htslib) reads BAM directly → numpy array via np.frombuffer (zero-copy)
  3. Label matching and HDF5 writing in the same process
  4. Output HDF5 format is identical to the old pipeline → MergeBin unchanged
"""

import sys
import os
import copy
import shlex
import logging
import numpy as np
from argparse import ArgumentParser, SUPPRESS

import libclair3
from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import (
    subprocess_popen, file_path_from,
    IUPAC_base_to_num_dict as BASE2NUM,
    IUPAC_base_to_ACGT_base_dict as BASE2BASE,
    str2bool,
)
from clair3.task.main import output_labels_from_vcf_columns, output_labels_from_reference

logging.basicConfig(format='%(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Shared utilities (label loading, HDF5 writing, non-variant filtering)
# ---------------------------------------------------------------------------

def variant_map_from(var_fn, tree, is_tree_empty):
    """Load truth variants from var_fn. Replicates clair3.utils.variant_map_from."""
    Y = {}
    truth_alt_dict = {}
    miss_variant_set = set()
    if var_fn is None:
        return Y, miss_variant_set, truth_alt_dict

    f = subprocess_popen(shlex.split("gzip -fdc %s" % (var_fn)))
    for row in f.stdout:
        if row[0] == "#":
            continue
        columns = row.strip().split()
        ctg_name, position_str, ref_base, alt_base, genotype1, genotype2 = columns
        key = ctg_name + ":" + position_str
        if genotype1 == '-1' or genotype2 == '-1':
            miss_variant_set.add(key)
            continue
        if not (is_tree_empty or is_region_in(tree, ctg_name, int(position_str))):
            continue

        Y[key] = output_labels_from_vcf_columns(columns)

        # Parse alt bases for read support checking
        ref_base_list, alt_base_list = _decode_alt(ref_base, alt_base)
        truth_alt_dict[int(position_str)] = (ref_base_list, alt_base_list)

    f.stdout.close()
    f.wait()
    return Y, miss_variant_set, truth_alt_dict


def _decode_alt(ref_base, alt_base):
    """Decode possibly multi-allelic alt bases."""
    if ',' not in alt_base:
        return [ref_base], [alt_base]
    alt_base = alt_base.split(',')
    ref_base_list, alt_base_list = [], []
    for ab in alt_base:
        rb, ab = _remove_common_suffix(ref_base, [ab])
        ref_base_list.append(rb)
        alt_base_list.append(ab[0])
    return ref_base_list, alt_base_list


def _remove_common_suffix(ref_base, alt_base):
    """Remove common suffix between ref and alt bases."""
    min_length = min(len(ref_base) - 1, min([len(item) - 1 for item in alt_base]))
    prefix = ref_base[::-1]
    for string in alt_base:
        string = string[::-1]
        while string[:len(prefix)] != prefix and prefix:
            prefix = prefix[:len(prefix) - 1]
        if not prefix:
            break
    res_length = len(prefix)
    if res_length > min_length:
        return ref_base, alt_base
    return ref_base[:len(ref_base) - res_length], [item[:len(item) - res_length] for item in alt_base]


def find_read_support(pos, truth_alt_dict, alt_info):
    """Check if truth variant has read support in the candidate. Replicates clair3.utils.find_read_support."""
    alt_info_parts = alt_info.rstrip().split('-')
    seqs = alt_info_parts[1].split(' ') if len(alt_info_parts) > 1 else ''
    seq_alt_bases_dict = dict(zip(seqs[::2], [int(item) for item in seqs[1::2]])) if len(seqs) else {}

    pos = int(pos)
    if pos not in truth_alt_dict:
        return None
    ref_base_list, alt_base_list = truth_alt_dict[pos]
    found = 0
    for alt_type in seq_alt_bases_dict:
        if '*' in alt_type or '#' in alt_type or 'R' in alt_type:
            continue
        if alt_type[0] == 'X':
            if alt_type[1] in alt_base_list:
                found += 1
        elif alt_type[0] == 'I':
            if alt_type[1:] in alt_base_list:
                found += 1
        elif alt_type[0] == 'D':
            del_cigar = alt_type[1:]
            for rb, ab in zip(ref_base_list, alt_base_list):
                if rb[1:] == del_cigar and len(ab) == 1:
                    found += 1
    if found >= len(alt_base_list):
        return True
    return False


def _filter_non_variants(X_keys, is_ref_flags, maximum_non_variant_ratio):
    """
    Filter non-variant samples to maintain variant:non-variant ratio.
    Returns a boolean mask of which samples to keep.
    """
    keep = np.ones(len(X_keys), dtype=bool)
    ref_indices = [i for i, is_ref in enumerate(is_ref_flags) if is_ref]
    non_variant_num = len(ref_indices)
    variant_num = len(X_keys) - non_variant_num
    if variant_num > 0 and non_variant_num > variant_num * maximum_non_variant_ratio:
        non_variant_keep_fraction = maximum_non_variant_ratio * variant_num / (1.0 * non_variant_num)
        probabilities = np.random.random_sample((non_variant_num,))
        for idx, p in zip(ref_indices, probabilities):
            if p > non_variant_keep_fraction:
                keep[idx] = False
    return keep


PREFIX_CHAR_STR = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _init_hdf5(bin_fn, tensor_shape, float_type, param):
    """Create HDF5 file with standard datasets. Format identical to Tensor2Bin output."""
    import h5py
    from clair3.utils import _hdf5_compression_kwargs, ensure_hdf5_plugin_path
    ensure_hdf5_plugin_path()
    compression_kwargs = _hdf5_compression_kwargs()

    chunk_rows = 500
    table_file = h5py.File(bin_fn, mode='w')
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
        dtype="S{}".format(param.no_of_positions + 50),
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
    return table_file


def _flush_phasing_to_hdf5(table_file, phasing_matrices, hp_labels, num_variants):
    """Write a batch of phasing data to HDF5."""
    if len(phasing_matrices) == 0:
        return
    from clair3.utils import _append_hdf5_dataset
    pm = np.array(phasing_matrices, dtype=np.int8)
    hp = np.array(hp_labels, dtype=np.int8)
    nv = np.array(num_variants, dtype=np.int16)
    _append_hdf5_dataset(table_file["phasing_matrix"], pm)
    _append_hdf5_dataset(table_file["hp_labels"], hp)
    _append_hdf5_dataset(table_file["phasing_num_variants"], nv)


def _flush_to_hdf5(table_file, tensors, labels, positions, alt_infos, tensor_shape, label_size, float_type):
    """Write a batch of tensors to HDF5. Matches clair3.utils.write_table_file behavior."""
    if len(tensors) == 0:
        return
    from clair3.utils import _append_hdf5_dataset

    position_matrix = np.array(tensors, dtype=np.dtype(float_type)).reshape([-1] + list(tensor_shape))
    position_dtype = table_file["position"].dtype
    alt_info_dtype = table_file["alt_info"].dtype
    position_arr = np.array(positions, dtype=position_dtype).reshape(-1, 1)
    alt_info_arr = np.array(alt_infos, dtype=alt_info_dtype).reshape(-1, 1)
    label_arr = np.array(labels, dtype=np.dtype(float_type)).reshape(-1, label_size)

    _append_hdf5_dataset(table_file["position_matrix"], position_matrix)
    _append_hdf5_dataset(table_file["position"], position_arr)
    _append_hdf5_dataset(table_file["alt_info"], alt_info_arr)
    _append_hdf5_dataset(table_file["label"], label_arr)


# ---------------------------------------------------------------------------
# Full-Alignment Direct Tensor Creation
# ---------------------------------------------------------------------------

def create_training_tensor_fa(args):
    """
    Create full-alignment training tensors using CFFI/C code directly.
    Replaces: CreateTrainingTensor (PyPy CreateTensorFullAlignment | CPython Tensor2Bin)
    """
    import shared.param_f as param

    ctg_name = args.ctgName
    bam_file_path = file_path_from(args.bam_fn)
    if bam_file_path is None or bam_file_path == "":
        logger.warning("[WARNING] Skip full-alignment variant calling for empty BAM path")
        return
    ref_fn = file_path_from(args.ref_fn, exit_on_not_found=True)
    var_fn = file_path_from(args.var_fn, exit_on_not_found=True)
    bin_fn = args.bin_fn
    bed_fn = file_path_from(args.bed_fn)
    extend_bed = file_path_from(args.extend_bed)
    full_aln_regions = file_path_from(args.full_aln_regions)
    phased_vcf_fn = args.phased_vcf_fn

    platform = args.platform
    enable_dwell_time = args.enable_dwell_time
    enable_long_indel = args.enable_long_indel
    add_no_phasing_data_training = args.add_no_phasing_data_training
    output_variant_matrix = getattr(args, 'output_variant_matrix', False)
    is_allow_duplicate_chr_pos = args.allow_duplicate_chr_pos
    shuffle = args.shuffle
    maximum_non_variant_ratio = args.maximum_non_variant_ratio

    min_mapping_quality = args.minMQ
    min_base_quality = args.minBQ

    no_phasing_for_fa = getattr(args, 'no_phasing_for_fa', False)
    need_haplotagging = not no_phasing_for_fa

    # Tensor dimensions
    matrix_depth = param.matrix_depth_dict[platform]
    BASE_CHANNEL_SIZE = param.channel_size
    channel_size = BASE_CHANNEL_SIZE + 1 if enable_dwell_time else BASE_CHANNEL_SIZE
    no_of_positions = param.no_of_positions
    extend_bp = param.extend_bp
    flanking_base_num = param.flankingBaseNum
    max_indel_length = param.maximum_variant_length_that_need_infer if not enable_long_indel else param.maximum_variant_length_that_need_infer_include_long_indel

    tensor_shape = [matrix_depth, no_of_positions, channel_size]
    float_type = 'int8'

    # --- Step 1: Load truth labels ---
    tree = bed_tree_from(bed_file_path=bed_fn)
    is_tree_empty = len(tree.keys()) == 0
    Y_true_var, miss_variant_set, truth_alt_dict = variant_map_from(var_fn, tree, is_tree_empty)
    Y = copy.deepcopy(Y_true_var)

    # --- Step 2: Parse full_aln_regions BED → candidate positions ---
    candidates_set = set()
    ctg_start = float('inf')
    ctg_end = 0

    if not full_aln_regions:
        logger.info("[INFO] No full_aln_regions provided, skipping")
        return

    with open(full_aln_regions, 'r') as f:
        for row in f:
            row = row.rstrip().split('\t')
            if row[0] != ctg_name:
                continue
            position = int(row[1]) + 1
            end = int(row[2]) + 1
            ctg_start = min(position, ctg_start)
            ctg_end = max(end, ctg_end)

            if len(row) > 3:  # hete snp positions (for haplotagging)
                pass  # handled by the C code via Variant structs
            else:
                if position == 1:
                    center = end - flanking_base_num - 2
                else:
                    center = position + (end - position) // 2 - 1
                candidates_set.add(center)

    if ctg_start == float('inf') or ctg_end == 0:
        logger.info("[INFO] No candidates found for %s, skipping", ctg_name)
        return

    # --- Step 3: Parse phased VCF → Variant structs for C haplotagging ---
    variant_list = []
    if need_haplotagging and phased_vcf_fn and os.path.exists(phased_vcf_fn):
        unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (phased_vcf_fn)))
        for row in unzip_process.stdout:
            row = row.rstrip()
            if row[0] == '#':
                continue
            columns = row.strip().split('\t')
            contig_name = columns[0]
            if ctg_name and contig_name != ctg_name:
                continue
            pos = int(columns[1])
            ref_base = columns[3]
            alt_base = columns[4]
            # Only single-base SNPs can be used for haplotagging (Variant struct uses char)
            if len(ref_base) != 1 or len(alt_base) != 1:
                continue
            genotype_info = columns[9].split(':')
            genotype, phase_set = genotype_info[0], genotype_info[-1]
            if '|' not in genotype:
                continue
            genotype = ('1' if genotype == '0|1' else '2')
            variant_list.append(libclair3.ffi.new(
                "struct Variant *",
                [pos - 1, ref_base.encode(), alt_base.encode(), int(genotype), int(phase_set)]
            ))
        variant_num = len(variant_list)
        Variants = libclair3.ffi.new("struct Variant *[]", variant_list)
    else:
        Variants = libclair3.ffi.new("struct Variant *[]", 1)
        variant_num = 0

    # 1-index to 0-index for C code
    candidates_list = sorted([item - 1 for item in candidates_set
                              if item >= ctg_start and item <= ctg_end])
    candidate_num = len(candidates_list)
    if candidate_num == 0:
        logger.info("[INFO] No candidates after filtering for %s, skipping", ctg_name)
        return

    # --- Step 4: Call C code to generate tensors ---
    region_str = '{}:{}-{}'.format(ctg_name, ctg_start, ctg_end).encode()
    candidates = libclair3.ffi.new("size_t [{}]".format(candidate_num), candidates_list)

    fa_data = libclair3.lib.calculate_clair3_full_alignment(
        region_str,
        bam_file_path.encode(),
        ref_fn.encode(),
        Variants,
        variant_num,
        candidates,
        candidate_num,
        need_haplotagging,
        min_mapping_quality,
        min_base_quality,
        matrix_depth,
        max_indel_length,
        enable_dwell_time,
        output_variant_matrix and need_haplotagging
    )

    # --- Step 5: Zero-copy numpy array from C buffer ---
    ffi = libclair3.ffi
    _dtype = np.int8
    size_sizet = np.dtype(_dtype).itemsize
    np_fa_data = np.frombuffer(
        ffi.buffer(
            fa_data.matrix,
            size_sizet * matrix_depth * no_of_positions * channel_size * candidate_num
        ),
        dtype=_dtype
    ).reshape(candidate_num, matrix_depth, no_of_positions, channel_size).copy()

    # --- Step 6: Extract position_info and alt_info from C data ---
    # Load reference for 33bp ref_seq construction (to match old pipeline's position format)
    import pysam
    ref_fasta = pysam.FastaFile(ref_fn)

    all_position_info = []  # "ctg:coord:33bp_ref_seq"
    all_alt_info = []
    all_coord_str = []  # just the coord number as string
    all_center_ref = []  # single center ref base
    for idx in range(candidate_num):
        alt_info_string = ffi.string(fa_data.all_alt_info[idx]).decode('utf8', 'ignore')
        alt_info = alt_info_string.rstrip().split('-')
        pos, depth, center_ref_base = alt_info[0], alt_info[1], alt_info[2]
        alt = alt_info[3] if len(alt_info) > 3 else ''

        # Fetch 33bp reference sequence around this position (1-based pos)
        pos_int = int(pos)
        ref_seq = ref_fasta.fetch(ctg_name, pos_int - flanking_base_num - 1,
                                  pos_int + flanking_base_num).upper()

        all_position_info.append(ctg_name + ':' + pos + ':' + ref_seq)
        all_alt_info.append(depth + '-' + alt)
        all_coord_str.append(pos)
        all_center_ref.append(center_ref_base)

    ref_fasta.close()

    # --- Step 6b: Extract variant matrix from C output (if requested) ---
    MAX_NEARBY_HETE_SNPS = 64  # must match C constant
    np_allele_matrix = None
    np_hp_labels = None
    np_num_nearby_variants = None
    if fa_data.has_variant_matrix:
        np_allele_matrix = np.frombuffer(
            ffi.buffer(fa_data.allele_matrix,
                       candidate_num * matrix_depth * MAX_NEARBY_HETE_SNPS),
            dtype=np.int8
        ).reshape(candidate_num, matrix_depth, MAX_NEARBY_HETE_SNPS).copy()

        np_hp_labels = np.frombuffer(
            ffi.buffer(fa_data.hp_labels, candidate_num * matrix_depth),
            dtype=np.int8
        ).reshape(candidate_num, matrix_depth).copy()

        np_num_nearby_variants = np.frombuffer(
            ffi.buffer(fa_data.num_nearby_variants, candidate_num * 2),
            dtype=np.int16
        ).reshape(candidate_num).copy()

    libclair3.lib.destroy_fa_data(fa_data)

    # --- Step 7: Label matching + filtering ---
    table_file = _init_hdf5(bin_fn, tensor_shape, float_type, param)

    # Create phasing datasets if variant matrix is available
    if np_allele_matrix is not None:
        from clair3.utils import _hdf5_compression_kwargs, ensure_hdf5_plugin_path
        ensure_hdf5_plugin_path()
        compression_kwargs = _hdf5_compression_kwargs()
        chunk_rows = 500
        table_file.create_dataset(
            "phasing_matrix", shape=(0, matrix_depth, MAX_NEARBY_HETE_SNPS),
            maxshape=(None, matrix_depth, MAX_NEARBY_HETE_SNPS),
            chunks=(chunk_rows, matrix_depth, MAX_NEARBY_HETE_SNPS),
            dtype=np.int8, **compression_kwargs)
        table_file.create_dataset(
            "hp_labels", shape=(0, matrix_depth),
            maxshape=(None, matrix_depth), chunks=(chunk_rows, matrix_depth),
            dtype=np.int8, **compression_kwargs)
        table_file.create_dataset(
            "phasing_num_variants", shape=(0,),
            maxshape=(None,), chunks=(chunk_rows,),
            dtype=np.int16, **compression_kwargs)

    # Buffers for batched HDF5 writing
    buf_tensors = []
    buf_labels = []
    buf_positions = []
    buf_alt_infos = []
    total_compressed = 0

    # Track keys seen for duplicate handling
    seen_keys = set()

    # Collect samples with label matching, then filter by non-variant ratio
    batch_entries = []  # list of (tensor, label, pos_str, alt_info, is_ref)
    batch_keys = []
    batch_is_ref = []

    for idx in range(candidate_num):
        pos_info = all_position_info[idx]
        alt_info = all_alt_info[idx]
        coord = all_coord_str[idx]
        center_ref = all_center_ref[idx]

        if center_ref not in 'ACGT':
            continue

        key = ctg_name + ":" + coord

        # Check bed region
        if not (is_tree_empty or is_region_in(tree, ctg_name, int(coord))):
            continue

        # Skip positions in miss_variant_set
        if key in miss_variant_set:
            continue

        # Check read support
        have_read_support = find_read_support(pos=coord, truth_alt_dict=truth_alt_dict, alt_info=alt_info)
        if have_read_support is not None and not have_read_support:
            miss_variant_set.add(key)
            continue

        is_reference = key not in Y_true_var

        # Handle duplicate keys
        actual_key = key
        if key in seen_keys:
            if is_allow_duplicate_chr_pos:
                for character in PREFIX_CHAR_STR:
                    tmp_key = character + key
                    if tmp_key not in seen_keys:
                        actual_key = tmp_key
                        break
                else:
                    continue  # all prefix chars exhausted
            else:
                continue
        seen_keys.add(actual_key)

        # Get label
        if is_reference and key not in Y:
            Y[key] = output_labels_from_reference(BASE2BASE[center_ref])
        label = Y.get(key)
        if label is None:
            continue

        # The C code already does center padding internally (sort_read_name_by_haplotype),
        # so np_fa_data[idx] is already [matrix_depth, no_of_positions, channel_size]
        tensor = np_fa_data[idx]

        pos_str = pos_info  # "ctg:coord:ref_seq"

        # Phasing data for this candidate (if available)
        phasing_entry = None
        if np_allele_matrix is not None:
            phasing_entry = (np_allele_matrix[idx], np_hp_labels[idx], np_num_nearby_variants[idx])

        batch_entries.append((tensor, label, pos_str, alt_info, is_reference, phasing_entry))
        batch_is_ref.append(is_reference)

        # Also generate unphased copy if requested
        if add_no_phasing_data_training:
            unphased_tensor = tensor.copy()
            # Channel 7 = phasing_info; set to 60 (unknown/unphased) for all reads
            # But only for rows that have actual read data (non-zero rows)
            unphased_tensor[:, :, 7] = np.where(
                np.any(unphased_tensor != 0, axis=(1, 2), keepdims=True).squeeze(axis=(1, 2)),
                60,
                0
            ).reshape(-1, 1)
            # Unphased copy: zero out phasing data
            unphased_phasing = None
            if phasing_entry is not None:
                zero_pm = np.zeros_like(phasing_entry[0])
                zero_hp = np.zeros_like(phasing_entry[1])
                unphased_phasing = (zero_pm, zero_hp, 0)

            unphased_key = None
            for character in PREFIX_CHAR_STR:
                tmp_key = character + key
                if tmp_key not in seen_keys:
                    unphased_key = tmp_key
                    break
            if unphased_key is not None:
                seen_keys.add(unphased_key)
                batch_entries.append((unphased_tensor, label, pos_str, alt_info, is_reference, unphased_phasing))
                batch_is_ref.append(is_reference)

    # Apply non-variant ratio filtering
    if maximum_non_variant_ratio is not None and len(batch_entries) > 0:
        keep = _filter_non_variants(
            [e[2] for e in batch_entries],
            batch_is_ref,
            maximum_non_variant_ratio
        )
    else:
        keep = np.ones(len(batch_entries), dtype=bool)

    # Shuffle if requested
    indices = np.where(keep)[0]
    if shuffle:
        np.random.shuffle(indices)

    # Write to HDF5
    buf_pm, buf_hp, buf_nv = [], [], []
    has_phasing_output = np_allele_matrix is not None

    for i in indices:
        tensor, label, pos_str, alt_info, _, phasing_entry = batch_entries[i]
        buf_tensors.append(tensor)
        buf_labels.append(label)
        buf_positions.append(pos_str)
        buf_alt_infos.append(alt_info)
        if has_phasing_output and phasing_entry is not None:
            buf_pm.append(phasing_entry[0])
            buf_hp.append(phasing_entry[1])
            buf_nv.append(phasing_entry[2])
        total_compressed += 1

        if total_compressed % 500 == 0 and total_compressed > 0:
            _flush_to_hdf5(table_file, buf_tensors, buf_labels, buf_positions,
                           buf_alt_infos, tensor_shape, param.label_size, float_type)
            if has_phasing_output and buf_pm:
                _flush_phasing_to_hdf5(table_file, buf_pm, buf_hp, buf_nv)
                buf_pm, buf_hp, buf_nv = [], [], []
            buf_tensors, buf_labels, buf_positions, buf_alt_infos = [], [], [], []

        if total_compressed % 50000 == 0:
            logger.info("[INFO] Compressed %d tensors", total_compressed)

    # Flush remaining
    if buf_tensors:
        _flush_to_hdf5(table_file, buf_tensors, buf_labels, buf_positions,
                       buf_alt_infos, tensor_shape, param.label_size, float_type)
    if has_phasing_output and buf_pm:
        _flush_phasing_to_hdf5(table_file, buf_pm, buf_hp, buf_nv)

    table_file.close()
    logger.info("[INFO] %s: wrote %d tensors to %s", ctg_name, total_compressed, bin_fn)


# ---------------------------------------------------------------------------
# Pileup Direct Tensor Creation
# ---------------------------------------------------------------------------

def create_training_tensor_pileup(args):
    """
    Create pileup training tensors using CFFI/C code directly.
    Replaces: CreateTrainingTensor --pileup (PyPy CreateTensorPileup | CPython Tensor2Bin)
    """
    import shared.param_p as param
    from preprocess.CreateTensorPileupFromCffi import (
        pileup_counts_clair3, BAMHandler
    )
    from preprocess.medaka_utils import Region

    ctg_name = args.ctgName
    bam_file_path = file_path_from(args.bam_fn, exit_on_not_found=True)
    ref_fn = file_path_from(args.ref_fn, exit_on_not_found=True)
    var_fn = file_path_from(args.var_fn, exit_on_not_found=True)
    bin_fn = args.bin_fn
    bed_fn = file_path_from(args.bed_fn)
    extend_bed = file_path_from(args.extend_bed)
    platform = args.platform

    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num

    minimum_snp_af_for_candidate = args.snp_min_af
    minimum_indel_af_for_candidate = args.indel_min_af
    min_coverage = args.minCoverage
    min_mapping_quality = args.minMQ

    is_allow_duplicate_chr_pos = args.allow_duplicate_chr_pos
    shuffle = args.shuffle
    maximum_non_variant_ratio = args.maximum_non_variant_ratio
    enable_long_indel = args.enable_long_indel

    flanking_base_num = param.flankingBaseNum
    no_of_positions = param.no_of_positions
    channel_size = param.channel_size

    fast_mode = getattr(args, 'fast_mode', False)
    call_snp_only = getattr(args, 'call_snp_only', False)

    tensor_shape = list(param.ont_input_shape if platform == 'ont' else param.input_shape)
    float_type = 'int32'

    # --- Step 1: Load truth labels ---
    tree = bed_tree_from(bed_file_path=bed_fn)
    is_tree_empty = len(tree.keys()) == 0
    Y_true_var, miss_variant_set, truth_alt_dict = variant_map_from(var_fn, tree, is_tree_empty)
    Y = copy.deepcopy(Y_true_var)

    # --- Step 2: Determine region ---
    confident_bed_fn = file_path_from(extend_bed)
    is_confident_bed_file_given = confident_bed_fn is not None
    confident_bed_tree, bed_start, bed_end = bed_tree_from(
        bed_file_path=extend_bed, contig_name=ctg_name, return_bed_region=True
    )

    fai_fn = file_path_from(ref_fn, suffix=".fai", exit_on_not_found=True, sep='.')

    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd

    fast_mode = platform == 'ont' and fast_mode
    minimum_snp_af_for_candidate = max(minimum_snp_af_for_candidate, param.min_af_dict[platform]) if fast_mode else minimum_snp_af_for_candidate
    min_coverage = max(min_coverage, 4) if fast_mode else min_coverage
    max_indel_length = param.maximum_variant_length_that_need_infer if not enable_long_indel else param.maximum_variant_length_that_need_infer_include_long_indel

    if not is_confident_bed_file_given and chunk_id is not None:
        contig_length = 0
        with open(fai_fn, 'r') as fai_fp:
            for row in fai_fp:
                columns = row.strip().split("\t")
                if columns[0] != ctg_name:
                    continue
                contig_length = int(columns[1])
        chunk_size = contig_length // chunk_num + 1 if contig_length % chunk_num else contig_length // chunk_num
        ctg_start = chunk_size * chunk_id
        ctg_end = ctg_start + chunk_size

    if is_confident_bed_file_given and chunk_id is not None:
        chunk_size = (bed_end - bed_start) // chunk_num + 1 if (bed_end - bed_start) % chunk_num else (bed_end - bed_start) // chunk_num
        ctg_start = bed_start + 1 + chunk_size * chunk_id
        ctg_end = ctg_start + chunk_size

    if ctg_start is None or ctg_end is None:
        logger.info("[INFO] No region determined for %s, skipping", ctg_name)
        return

    ctg_start = max(1, ctg_start)
    extend_start = max(1, ctg_start - no_of_positions)
    extend_end = ctg_end + no_of_positions

    region_str = "{}:{}-{}".format(ctg_name, extend_start, extend_end)
    region = Region.from_string(region_str)

    confident_bed_tree_for_filter = bed_tree_from(
        bed_file_path=confident_bed_fn, contig_name=ctg_name,
        bed_ctg_start=extend_start, bed_ctg_end=extend_end
    )

    # --- Step 3: Call C pileup code ---
    chunk_result, all_alt_info_list, _ = pileup_counts_clair3(
        region,
        bam=bam_file_path,
        fasta=ref_fn,
        min_depth=min_coverage,
        min_snp_af=minimum_snp_af_for_candidate,
        min_indel_af=minimum_indel_af_for_candidate,
        min_mq=min_mapping_quality,
        max_indel_length=max_indel_length,
        call_snp_only=call_snp_only,
        max_depth=param.max_depth,
        gvcf=False,
    )

    # --- Step 4: Extract per-candidate tensors from pileup result ---
    # Load reference for 33bp ref_seq construction
    import pysam
    ref_fasta = pysam.FastaFile(ref_fn)

    np_pileup_data = []
    all_position_info = []  # "ctg:coord:33bp_ref_seq"
    all_alt_info = []
    all_coord = []
    all_center_ref = []

    for idx, (pos, pos_info, alt_info) in enumerate(all_alt_info_list):
        pos = int(pos)
        pass_confident_bed = not is_confident_bed_file_given or is_region_in(
            tree=confident_bed_tree_for_filter, contig_name=ctg_name,
            region_start=pos - 1, region_end=pos + 1
        )
        if not pass_confident_bed:
            continue

        start, end = pos - flanking_base_num, pos + flanking_base_num + 1
        for result in chunk_result:
            if start - 1 >= result[1][0][0] and end <= result[1][-1][0]:
                offset = start - result[1][0][0] - 1
                tensor = result[0][offset: offset + no_of_positions]
                if tensor.shape != (no_of_positions, channel_size):
                    continue
                # Check for empty columns
                if np.sum(np.sum(tensor == 0, axis=1) == channel_size) > 0:
                    continue
                np_pileup_data.append(tensor)
                # Fetch 33bp reference sequence (pos is 1-based)
                ref_seq = ref_fasta.fetch(ctg_name, pos - flanking_base_num - 1,
                                          pos + flanking_base_num).upper()
                center_ref = ref_seq[flanking_base_num] if len(ref_seq) > flanking_base_num else ''
                all_position_info.append(ctg_name + ':' + str(pos) + ':' + ref_seq)
                all_alt_info.append(alt_info)
                all_coord.append(str(pos))
                all_center_ref.append(center_ref)

    ref_fasta.close()

    if len(np_pileup_data) == 0:
        logger.info("[INFO] No pileup tensors for %s (chunk %s/%s)", ctg_name, chunk_id, chunk_num)
        return

    np_pileup_data = np.array(np_pileup_data, dtype=np.int32)

    # --- Step 5: Label matching + filtering ---
    table_file = _init_hdf5(bin_fn, tensor_shape, float_type, param)

    seen_keys = set()
    batch_entries = []
    batch_is_ref = []

    for idx in range(len(np_pileup_data)):
        pos_info = all_position_info[idx]
        alt_info = all_alt_info[idx]
        coord = all_coord[idx]
        center_ref = all_center_ref[idx]

        if center_ref not in 'ACGT':
            continue

        key = ctg_name + ":" + coord

        if not (is_tree_empty or is_region_in(tree, ctg_name, int(coord))):
            continue

        if key in miss_variant_set:
            continue

        have_read_support = find_read_support(pos=coord, truth_alt_dict=truth_alt_dict, alt_info=alt_info)
        if have_read_support is not None and not have_read_support:
            miss_variant_set.add(key)
            continue

        is_reference = key not in Y_true_var

        actual_key = key
        if key in seen_keys:
            if is_allow_duplicate_chr_pos:
                for character in PREFIX_CHAR_STR:
                    tmp_key = character + key
                    if tmp_key not in seen_keys:
                        actual_key = tmp_key
                        break
                else:
                    continue
            else:
                continue
        seen_keys.add(actual_key)

        if is_reference and key not in Y:
            Y[key] = output_labels_from_reference(BASE2BASE[center_ref])
        label = Y.get(key)
        if label is None:
            continue

        tensor = np_pileup_data[idx]
        batch_entries.append((tensor, label, pos_info, alt_info, is_reference))
        batch_is_ref.append(is_reference)

    # Apply non-variant ratio filtering
    if maximum_non_variant_ratio is not None and len(batch_entries) > 0:
        keep = _filter_non_variants(
            [e[2] for e in batch_entries],
            batch_is_ref,
            maximum_non_variant_ratio
        )
    else:
        keep = np.ones(len(batch_entries), dtype=bool)

    indices = np.where(keep)[0]
    if shuffle:
        np.random.shuffle(indices)

    # Write to HDF5
    buf_tensors, buf_labels, buf_positions, buf_alt_infos = [], [], [], []
    total_compressed = 0

    for i in indices:
        tensor, label, pos_str, alt_info, _ = batch_entries[i]
        buf_tensors.append(tensor)
        buf_labels.append(label)
        buf_positions.append(pos_str)
        buf_alt_infos.append(alt_info)
        total_compressed += 1

        if total_compressed % 500 == 0 and total_compressed > 0:
            _flush_to_hdf5(table_file, buf_tensors, buf_labels, buf_positions,
                           buf_alt_infos, tensor_shape, param.label_size, float_type)
            buf_tensors, buf_labels, buf_positions, buf_alt_infos = [], [], [], []

        if total_compressed % 50000 == 0:
            logger.info("[INFO] Compressed %d tensors", total_compressed)

    if buf_tensors:
        _flush_to_hdf5(table_file, buf_tensors, buf_labels, buf_positions,
                       buf_alt_infos, tensor_shape, param.label_size, float_type)

    table_file.close()
    logger.info("[INFO] %s (chunk %s/%s): wrote %d pileup tensors to %s",
                ctg_name, chunk_id, chunk_num, total_compressed, bin_fn)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def Run(args):
    if args.pileup:
        create_training_tensor_pileup(args)
    else:
        create_training_tensor_fa(args)


def main():
    parser = ArgumentParser(
        description="Create training tensor binaries directly using CFFI/C (no PyPy, no text pipe)")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="bam.bam", required=True,
                        help="BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa", required=True,
                        help="Reference fasta file input, required")

    parser.add_argument('--var_fn', type=str, default=None, required=True,
                        help="Truth variant file from GetTruth, required")

    parser.add_argument('--bin_fn', type=str, default=None, required=True,
                        help="Output HDF5 binary file, required")

    parser.add_argument('--bed_fn', type=str, nargs='?', action="store", default=None,
                        help="High-confidence BED regions")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="Contig name to process, required")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="1-based start position")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="1-based end position")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to samtools, default: %(default)s")

    # AF thresholds
    parser.add_argument('--min_af', type=float, default=None,
                        help="Minimum allele frequency, default: %(default)s")

    parser.add_argument('--snp_min_af', type=float, default=0.08,
                        help="Minimum SNP AF, default: %(default)f")

    parser.add_argument('--indel_min_af', type=float, default=0.08,
                        help="Minimum Indel AF, default: %(default)f")

    parser.add_argument('--minCoverage', type=int, default=2,
                        help="Minimum coverage, default: %(default)d")

    parser.add_argument('--minMQ', type=int, default=5,
                        help="Minimum mapping quality, default: %(default)d")

    parser.add_argument('--minBQ', type=int, default=0,
                        help="Minimum base quality, default: %(default)d")

    # Training options
    parser.add_argument('--maximum_non_variant_ratio', default=None, type=float,
                        help='Max ratio of non-variants to variants')

    parser.add_argument('--shuffle', action='store_true',
                        help="Shuffle output tensors")

    parser.add_argument('--allow_duplicate_chr_pos', action='store_true',
                        help="Allow duplicate chr:pos in output")

    parser.add_argument('--pileup', action='store_true',
                        help="Pileup mode (default: full-alignment mode)")

    parser.add_argument('--enable_long_indel', action='store_true',
                        help="Enable long indel support")

    # FA-specific options
    parser.add_argument('--extend_bed', nargs='?', action="store", type=str, default=None,
                        help="Extended BED regions for tensor creation")

    parser.add_argument('--full_aln_regions', type=str, default=None,
                        help="Full-alignment candidate regions BED")

    parser.add_argument('--phased_vcf_fn', type=str, default=None,
                        help="Phased VCF for haplotagging")

    parser.add_argument('--no_phasing_for_fa', type=str2bool, default=False,
                        help="Disable haplotagging for FA tensors")

    parser.add_argument('--phasing_info_in_bam', action='store_true',
                        help="Use HP tags from BAM for phasing")

    parser.add_argument('--add_no_phasing_data_training', action='store_true',
                        help="Generate additional unphased tensor copies")

    parser.add_argument('--enable_dwell_time', action='store_true',
                        help="Enable nanopore dwell time channel")

    parser.add_argument('--output_variant_matrix', action='store_true',
                        help="Output per-read allele state matrix at nearby het SNPs (for neural phasing training)")

    # Chunking (for parallel execution)
    parser.add_argument('--chunk_num', type=int, default=None, help=SUPPRESS)
    parser.add_argument('--chunk_id', type=int, default=None, help=SUPPRESS)

    # Pileup-specific
    parser.add_argument('--fast_mode', type=str2bool, default=False, help=SUPPRESS)
    parser.add_argument('--call_snp_only', type=str2bool, default=False, help=SUPPRESS)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()
