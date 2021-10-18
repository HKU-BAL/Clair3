import sys
import gc
import copy
import shlex
import os
import tables
import numpy as np
from functools import partial

from clair3.task.main import *
from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import subprocess_popen, IUPAC_base_to_ACGT_base_dict as BASE2BASE, IUPAC_base_to_num_dict as BASE2NUM

FILTERS = tables.Filters(complib='blosc:lz4hc', complevel=5)
shuffle_bin_size = 50000
PREFIX_CHAR_STR = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"


def setup_environment():
    gc.enable()


def batches_from(iterable, item_from, batch_size=1):
    iterable = iter(iterable)
    while True:
        chunk = []
        for _ in range(batch_size):
            try:
                chunk.append(item_from(next(iterable)))
            except StopIteration:
                yield chunk
                return
        yield chunk


def tensor_generator_from(tensor_file_path, batch_size, pileup, platform):
    global param
    float_type = 'int32'
    if pileup:
        import shared.param_p as param
    else:
        import shared.param_f as param
        float_type = 'int8'

    if tensor_file_path != "PIPE":
        f = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, tensor_file_path)))
        fo = f.stdout
    else:
        fo = sys.stdin

    processed_tensors = 0
    tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape
    prod_tensor_shape = np.prod(tensor_shape)

    def item_from(row):
        chrom, coord, seq, tensor, alt_info = row.split("\t")
        if pileup:
            tensor = np.array(tensor.split(), dtype=np.dtype(float_type))
            depth = int(alt_info.split('-', maxsplit=1)[0])
            max_depth = param.max_depth_dict[platform]
            # for extreme high coverage data, make sure we could have a truncated coverage
            if depth > 0 and depth > max_depth * 1.5:
                scale_factor = depth / max_depth
                tensor = tensor / scale_factor
        else:
            # need add padding if depth is lower than maximum depth.
            tensor = [int(item) for item in tensor.split()]
            tensor_depth = len(tensor) // tensor_shape[1] // tensor_shape[2]
            padding_depth = tensor_shape[0] - tensor_depth
            prefix_padding_depth = int(padding_depth / 2)
            suffix_padding_depth = padding_depth - int(padding_depth / 2)
            prefix_zero_padding = [0] * prefix_padding_depth * tensor_shape[1] * tensor_shape[2]
            suffix_zero_padding = [0] * suffix_padding_depth * tensor_shape[1] * tensor_shape[2]
            tensor = prefix_zero_padding + tensor + suffix_zero_padding
            tensor = np.array(tensor, dtype=np.dtype(float_type))

        pos = chrom + ":" + coord + ":" + seq
        return tensor, pos, seq, alt_info

    for batch in batches_from(fo, item_from=item_from, batch_size=batch_size):
        tensors = np.empty(([batch_size, prod_tensor_shape]), dtype=np.dtype(float_type))
        positions = []
        alt_info_list = []
        for tensor, pos, seq, alt_info in batch:
            if seq[param.flankingBaseNum] not in BASE2NUM:
                continue
            tensors[len(positions)] = tensor
            positions.append(pos)
            alt_info_list.append(alt_info)

        current_batch_size = len(positions)
        X = np.reshape(tensors, ([batch_size] + tensor_shape))

        if processed_tensors > 0 and processed_tensors % 20000 == 0:
            print("Processed %d tensors" % processed_tensors, file=sys.stderr)

        processed_tensors += current_batch_size

        if current_batch_size <= 0:
            continue
        yield X[:current_batch_size], positions[:current_batch_size], alt_info_list[:current_batch_size]

    if tensor_file_path != "PIPE":
        fo.close()
        f.wait()


def remove_common_suffix(ref_base, alt_base):
    min_length = min(len(ref_base) - 1, min([len(item) - 1 for item in alt_base]))  # keep at least one base
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

    return ref_base[-min_length], [item[-min_length] for item in alt_base]


def decode_alt(ref_base, alt_base):
    if ',' not in alt_base:
        return [ref_base], [alt_base]
    alt_base = alt_base.split(',')
    ref_base_list, alt_base_list = [], []
    for ab in alt_base:
        rb,ab = remove_common_suffix(ref_base, [ab])
        ref_base_list.append(rb)
        alt_base_list.append(ab[0])
    return ref_base_list, alt_base_list


def variant_map_from(var_fn, tree, is_tree_empty):
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
        ref_base_list, alt_base_list = decode_alt(ref_base, alt_base)
        truth_alt_dict[int(position_str)] = (ref_base_list, alt_base_list)
    f.stdout.close()
    f.wait()
    return Y, miss_variant_set, truth_alt_dict

def find_read_support(pos, truth_alt_dict, alt_info):
    alt_info = alt_info.rstrip().split('-')
    seqs = alt_info[1].split(' ') if len(alt_info) > 1 else ''
    seq_alt_bases_dict = dict(zip(seqs[::2], [int(item) for item in seqs[1::2]])) if len(seqs) else {}

    pos = int(pos)
    if pos not in truth_alt_dict:
        # candidate position not in the truth vcf or unified truth vcf
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
    # return False if we find any alternative bases missed in subsampled bam, then remove the position from training
    return False

def write_table_dict(table_dict, string, label, pos, total, alt_info, tensor_shape, pileup):
    """
    Write pileup or full alignment tensor into a dictionary.compressed bin file.
    table_dict: dictionary include all training information (tensor position, label, altnative bases).
    string: input tensor string, need add padding to meet the depth requirement.
    label: include gt21 genotype, indel length 1, indel length 2.
    alt_info: altnative information for querying variant.
    """

    if len(string) == 1:
        string = string[0]
    position_matrix = string
    position_matrix = position_matrix.split()

    if pileup:
        table_dict['position_matrix'].append(position_matrix)
    else:
        tensor_depth = len(position_matrix) // tensor_shape[1] // tensor_shape[2]
        padding_depth = tensor_shape[0] - tensor_depth
        prefix_padding_depth = int(padding_depth / 2)
        suffix_padding_depth = padding_depth - int(padding_depth / 2)
        prefix_zero_padding = ['0'] * prefix_padding_depth * tensor_shape[1] * tensor_shape[2]
        suffix_zero_padding = ['0'] * suffix_padding_depth * tensor_shape[1] * tensor_shape[2]
        table_dict['position_matrix'].append(prefix_zero_padding + position_matrix + suffix_zero_padding)

    table_dict['position'].append(pos)
    table_dict['label'].append(label)
    table_dict['alt_info'].append(alt_info)

    return total + 1


def update_table_dict():
    table_dict = {}
    table_dict['position_matrix'] = []
    table_dict['alt_info'] = []
    table_dict['position'] = []
    table_dict['label'] = []
    return table_dict


def write_table_file(table_file, table_dict, tensor_shape, label_size, float_type):
    """
    Write pileup or full alignment tensor into compressed bin file.
    table_dict: dictionary include all training information (tensor position, label, altnative bases).
    string: input tensor string, need add padding to meet the depth requirement.
    tree: dictionary(contig name : intervaltree) for quick region querying.
    miss_variant_set:  sometimes there will have true variant missing after downsampling reads.
    is_allow_duplicate_chr_pos: whether allow duplicate positions when training, if there exists downsampled data, lower depth will add a random prefix character.
    non_variant_subsample_ratio: define a maximum non variant ratio for training, we always expect use more non variant data, while it would greatly increase training
    time, especially in ont data, here we usually use 1:1 or 1:2 for variant candidate: non variant candidate.
    """

    position_matrix = np.array(table_dict['position_matrix'], np.dtype(float_type)).reshape([-1] + tensor_shape)
    table_file.root.position_matrix.append(position_matrix)

    table_file.root.alt_info.append(np.array(table_dict['alt_info']).reshape(-1, 1))
    table_file.root.position.append(np.array(table_dict['position']).reshape(-1, 1))
    table_file.root.label.append(np.array(table_dict['label'], np.dtype(float_type)).reshape(-1, label_size))
    table_dict = update_table_dict()

    return table_dict


def print_bin_size(path, prefix=None):
    import tables
    import os
    total = 0
    for file_name in os.listdir(path):
        if prefix and not file_name.startswith(prefix):
            continue
        table = tables.open_file(os.path.join(path, file_name), 'r')
        print("[INFO] {} size is: {}".format(file_name, len(table.root.label)))
        total += len(table.root.label)
    print('[INFO] total: {}'.format(total))


def bin_reader_generator_from(tensor_fn, Y_true_var, Y, is_tree_empty, tree, miss_variant_set, truth_alt_dict, is_allow_duplicate_chr_pos=False, maximum_non_variant_ratio=None):

    """
    Bin reader generator for bin file generation.
    tensor_fn: tensor file.
    Y_true_var: dictionary (contig name: label information) containing all true variant information (should not be changed).
    Y: dictionary (contig name: label information) to store all variant and non variant information.
    tree: dictionary(contig name : intervaltree) for quick region querying.
    miss_variant_set:  sometimes there will have true variant missing after downsampling reads.
    truth_alt_dict: unified truth reference base and alternative bases to find read support.
    is_allow_duplicate_chr_pos: whether allow duplicate positions when training, if there exists downsampled data, lower depth will add a random prefix character.
    maximum_non_variant_ratio: define a maximum non variant ratio for training, we always expect use more non variant data, while it would greatly increase training
    time, especially in ont data, here we usually use 1:1 or 1:2 for variant candidate: non variant candidate.
    """

    X = {}
    ref_list = []
    total = 0
    variant_set_with_read_support = set()
    variants_without_read_support = 0
    for row_idx, row in enumerate(tensor_fn):
        chrom, coord, seq, string, alt_info = row.split("\t")
        alt_info = alt_info.rstrip()
        if not (is_tree_empty or is_region_in(tree, chrom, int(coord))):
            continue
        seq = seq.upper()
        if seq[param.flankingBaseNum] not in 'ACGT':
            continue
        key = chrom + ":" + coord
        is_reference = key not in Y_true_var

        if key in miss_variant_set:
            continue

        have_read_support = find_read_support(pos=coord, truth_alt_dict=truth_alt_dict, alt_info=alt_info)
        if have_read_support is not None and not have_read_support:
            miss_variant_set.add(key)
            variants_without_read_support += 1
            continue

        variant_set_with_read_support.add(key)
        if key not in X:
            X[key] = (string, alt_info, seq)
            if is_reference:
                ref_list.append(key)
        elif is_allow_duplicate_chr_pos:
            new_key = ""
            for character in PREFIX_CHAR_STR:
                tmp_key = character + key
                if tmp_key not in X:
                    new_key = tmp_key
                    break
            if len(new_key) > 0:
                X[new_key] = (string, alt_info, seq)
            if is_reference:
                ref_list.append(new_key)

        if is_reference and key not in Y:
            Y[key] = output_labels_from_reference(BASE2BASE[seq[param.flankingBaseNum]])

        if len(X) == shuffle_bin_size:
            if maximum_non_variant_ratio is not None:
                _filter_non_variants(X, ref_list, maximum_non_variant_ratio)
            yield X, total, False
            X = {}
            ref_list = []
        total += 1
        if total % 100000 == 0:
            print("[INFO] Processed %d tensors" % total, file=sys.stderr)

    print("[INFO] Variants with read support/variants without read support: {}/{}".format(len(variant_set_with_read_support), variants_without_read_support))
    if maximum_non_variant_ratio is not None:
        _filter_non_variants(X, ref_list, maximum_non_variant_ratio)
    yield X, total, True


def _filter_non_variants(X, ref_list, maximum_non_variant_ratio):
    non_variant_num = len(ref_list)
    variant_num = len(X) - non_variant_num
    if non_variant_num > variant_num * maximum_non_variant_ratio:
        non_variant_keep_fraction = maximum_non_variant_ratio * variant_num / (1. * non_variant_num)
        probabilities = np.random.random_sample((non_variant_num,))
        for key, p in zip(ref_list, probabilities):
            if p > non_variant_keep_fraction:
                X.pop(key)


def get_training_array(tensor_fn, var_fn, bed_fn, bin_fn, shuffle=True, is_allow_duplicate_chr_pos=True, chunk_id=None,
                       chunk_num=None, platform='ont', pileup=False, maximum_non_variant_ratio=None, candidate_details_fn_prefix=None):

    """
    Generate training array for training. here pytables with blosc:lz4hc are used for extreme fast compression and decompression,
    which can meet the requirement of gpu utilization. lz4hc decompression allows speed up training array decompression 4~5x compared
    with tensorflow tfrecord file format, current gpu utilization could reach over 85% with only 10G memory.
    tensor_fn: string format tensor acquired from CreateTensorPileup or CreateTensorFullAlign, include contig name position, tensor matrix, alternative information.
    var_fn: simplified variant(vcf) format from GetTruths, which include contig name, position, reference base, alternative base, genotype.
    bin_fn: pytables format output bin file name.
    shuffle: whether apply index shuffling when generating training data, default True, which would promote robustness.
    is_allow_duplicate_chr_pos: whether allow duplicate positions when training, if there exists downsampled data, lower depth will add a random prefix character.
    chunk_id: specific chunk id works with total chunk_num for parallel execution. Here will merge all tensor file with sampe prefix.
    chunk_num: total chunk number for parallel execution. Each chunk refer to a smaller reference regions.
    platform: platform for tensor shape, ont give a larger maximum depth compared with pb and illumina.
    pileup: whether in pileup mode. Define two calling mode, pileup or full alignment.
    maximum_non_variant_ratio: define a maximum non variant ratio for training, we always expect use more non variant data, while it would greatly increase training
    time, especially in ont data, here we usually use 1:1 or 1:2 for variant candidate: non variant candidate.
    candidate_details_fn_prefix: a counter to calculate total variant and non variant from the information in alternative file.
    """

    tree = bed_tree_from(bed_file_path=bed_fn)
    is_tree_empty = len(tree.keys()) == 0
    Y_true_var, miss_variant_set, truth_alt_dict = variant_map_from(var_fn, tree, is_tree_empty)
    Y = copy.deepcopy(Y_true_var)

    global param
    float_type = 'int32'
    if pileup:
        import shared.param_p as param
    else:
        import shared.param_f as param
        float_type = 'int8'

    tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape

    subprocess_list = []
    if tensor_fn == 'PIPE':
        subprocess_list.append(sys.stdin)
    elif os.path.exists(tensor_fn):
        subprocess_list.append(subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, tensor_fn))).stdout)
    # select all match prefix if file path not exists
    else:
        tensor_fn = tensor_fn.split('/')
        directry, file_prefix = '/'.join(tensor_fn[:-1]), tensor_fn[-1]
        all_file_name = []
        for file_name in os.listdir(directry):
            if file_name.startswith(file_prefix + '_') or file_name.startswith(
                    file_prefix + '.'):  # add '_.' to avoid add other prefix chr
                all_file_name.append(file_name)
        all_file_name = sorted(all_file_name)
        if chunk_id is not None:
            chunk_size = len(all_file_name) // chunk_num if len(all_file_name) % chunk_num == 0 else len(
                all_file_name) // chunk_num + 1
            chunk_start = chunk_size * chunk_id
            chunk_end = chunk_start + chunk_size
            all_file_name = all_file_name[chunk_start:chunk_end]
        if not len(all_file_name):
            print("[INFO] chunk_id exceed total file number, skip chunk", file=sys.stderr)
            return 0
        for file_name in all_file_name:
            subprocess_list.append(
                subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, os.path.join(directry, file_name)))).stdout)

    tables.set_blosc_max_threads(64)
    int_atom = tables.Atom.from_dtype(np.dtype(float_type))
    string_atom = tables.StringAtom(itemsize=param.no_of_positions + 50)
    long_string_atom = tables.StringAtom(itemsize=5000)  # max alt_info length
    table_file = tables.open_file(bin_fn, mode='w', filters=FILTERS)
    table_file.create_earray(where='/', name='position_matrix', atom=int_atom, shape=[0] + tensor_shape,
                             filters=FILTERS)
    table_file.create_earray(where='/', name='position', atom=string_atom, shape=(0, 1), filters=FILTERS)
    table_file.create_earray(where='/', name='label', atom=int_atom, shape=(0, param.label_size), filters=FILTERS)
    table_file.create_earray(where='/', name='alt_info', atom=long_string_atom, shape=(0, 1), filters=FILTERS)

    table_dict = update_table_dict()

    # generator to avoid high memory occupy
    bin_reader_generator = partial(bin_reader_generator_from, 
                                   Y_true_var=Y_true_var,
                                   Y=Y,
                                   is_tree_empty=is_tree_empty,
                                   tree=tree,
                                   miss_variant_set=miss_variant_set,
                                   truth_alt_dict=truth_alt_dict,
                                   is_allow_duplicate_chr_pos=is_allow_duplicate_chr_pos,
                                   maximum_non_variant_ratio=maximum_non_variant_ratio)

    total_compressed = 0
    for fin in subprocess_list:
        bin_g = bin_reader_generator(tensor_fn=fin)
        completed = False
        while not completed:
            try:
                X, total, completed = next(bin_g)
            except StopIteration:
                completed = True
        
            if X is None or not len(X):
                break
            all_chr_pos = sorted(X.keys())
            if shuffle == True:
                np.random.shuffle(all_chr_pos)
            for key in all_chr_pos:

                string, alt_info, seq = X[key]
                del X[key]
                label = None
                if key in Y:
                    label = Y[key]
                    pos = key + ':' + seq
                    if not is_allow_duplicate_chr_pos:
                        del Y[key]
                elif is_allow_duplicate_chr_pos:
                    tmp_key = key[1:]
                    label = Y[tmp_key]
                    pos = tmp_key + ':' + seq
                if label is None:
                    print(key)
                    continue
                total_compressed = write_table_dict(table_dict, string, label, pos, total_compressed, alt_info,
                                                    tensor_shape, pileup)

                if total_compressed % 500 == 0 and total_compressed > 0:
                    table_dict = write_table_file(table_file, table_dict, tensor_shape, param.label_size, float_type)

                if total_compressed % 50000 == 0:
                    print("[INFO] Compressed %d tensor" % (total_compressed), file=sys.stderr)
        fin.close()

    if total_compressed % 500 != 0 and total_compressed > 0:
        table_dict = write_table_file(table_file, table_dict, tensor_shape, param.label_size, float_type)

    table_file.close()
    print("[INFO] Compressed %d/%d tensor" % (total_compressed, total), file=sys.stderr)

