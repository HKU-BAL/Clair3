from collections import namedtuple

from clair3.task.genotype import Genotype, genotype_enum_from, genotype_enum_for_task
from clair3.task.gt21 import *
from clair3.task.variant_length import VariantLength

OutputLabelNamedTuple = namedtuple(
    'BasePredictNamedTuple', ['output_label_count', 'y_start_index', 'y_end_index']
)
GT21 = OutputLabelNamedTuple(
    output_label_count=21,
    y_start_index=0,
    y_end_index=21,
)
GENOTYPE = OutputLabelNamedTuple(
    output_label_count=3,
    y_start_index=GT21.y_end_index,
    y_end_index=GT21.y_end_index + 3,
)
VARIANT_LENGTH_1 = OutputLabelNamedTuple(
    output_label_count=VariantLength.output_label_count,
    y_start_index=GENOTYPE.y_end_index,
    y_end_index=GENOTYPE.y_end_index + VariantLength.output_label_count,
)
VARIANT_LENGTH_2 = OutputLabelNamedTuple(
    output_label_count=VariantLength.output_label_count,
    y_start_index=VARIANT_LENGTH_1.y_end_index,
    y_end_index=VARIANT_LENGTH_1.y_end_index + VariantLength.output_label_count,
)


def min_max(value, minimum, maximum):
    return max(min(value, maximum), minimum)


def output_labels_from_reference(reference_base):
    gt21_vec = [0] * GT21.output_label_count
    gt21_vec[gt21_enum_from_label(reference_base + reference_base)] = 1

    genotype_vec = [0] * GENOTYPE.output_label_count
    genotype_vec[Genotype.homo_reference] = 1

    variant_length_vec_1 = [0] * VARIANT_LENGTH_1.output_label_count
    variant_length_vec_2 = [0] * VARIANT_LENGTH_2.output_label_count
    variant_length_vec_1[0 + VariantLength.index_offset] = 1
    variant_length_vec_2[0 + VariantLength.index_offset] = 1

    return gt21_vec + genotype_vec + variant_length_vec_1 + variant_length_vec_2


def output_labels_from_vcf_columns(columns, homo_calling=False, haplotype=None):
    reference, alternate = columns[2], columns[3]
    genotype_1, genotype_2 = int(columns[4]), int(columns[5])

    alternate_arr = alternate.split(',')
    if len(alternate_arr) == 1:
        alternate_arr = (
            [reference if genotype_1 == 0 or genotype_2 == 0 else alternate_arr[0]] +
            alternate_arr
        )

    gt21 = gt21_enum_from(reference, alternate, genotype_1, genotype_2, alternate_arr)
    gt21_vec = [0] * GT21.output_label_count
    gt21_vec[gt21] = 1

    genotype = genotype_enum_from(genotype_1, genotype_2)
    genotype_for_task = genotype_enum_for_task(genotype)
    genotype_vec = [0] * GENOTYPE.output_label_count
    genotype_vec[genotype_for_task] = 1

    variant_lengths = [
        min_max(len(alt) - len(reference), VariantLength.min, VariantLength.max)
        for alt in alternate_arr
    ]
    variant_lengths.sort()
    variant_length_vec_1 = [0] * VARIANT_LENGTH_1.output_label_count
    variant_length_vec_2 = [0] * VARIANT_LENGTH_2.output_label_count
    variant_length_vec_1[variant_lengths[0] + VariantLength.index_offset] = 1
    variant_length_vec_2[variant_lengths[1] + VariantLength.index_offset] = 1

    return gt21_vec + genotype_vec + variant_length_vec_1 + variant_length_vec_2

def output_labels_from_reference_new(reference_base, base_idx):
    gt21_vec = [0] * GT21.output_label_count
    gt21_vec[gt21_enum_from_label(reference_base + reference_base)] = 1

    genotype_vec = [0] * GENOTYPE.output_label_count + [0]
    if base_idx == '2':
        genotype_vec[Genotype.homo_reference] = 1
    elif base_idx == '1':
        genotype_vec[3] = 1
    variant_length_vec_1 = [0] * VARIANT_LENGTH_1.output_label_count
    variant_length_vec_2 = [0] * VARIANT_LENGTH_2.output_label_count
    variant_length_vec_1[0 + VariantLength.index_offset] = 1
    variant_length_vec_2[0 + VariantLength.index_offset] = 1

    return gt21_vec + genotype_vec + variant_length_vec_1 + variant_length_vec_2


def output_labels_from_vcf_columns_new(columns):
    reference, alternate = columns[2], columns[3]
    genotype_1, genotype_2 = int(columns[4]), int(columns[5])

    alternate_arr = alternate.split(',')
    if len(alternate_arr) == 1:
        alternate_arr = (
            [reference if genotype_1 == 0 or genotype_2 == 0 else alternate_arr[0]] +
            alternate_arr
        )

    gt21 = gt21_enum_from(reference, alternate, genotype_1, genotype_2, alternate_arr)
    gt21_vec = [0] * GT21.output_label_count
    gt21_vec[gt21] = 1

    genotype = genotype_enum_from(genotype_1, genotype_2)
    genotype_for_task = genotype_enum_for_task(genotype)
    genotype_vec = [0] * GENOTYPE.output_label_count
    genotype_vec[genotype_for_task] = 1

    genotype_vec += [0]

    variant_lengths = [
        min_max(len(alt) - len(reference), VariantLength.min, VariantLength.max)
        for alt in alternate_arr
    ]
    variant_lengths.sort()
    variant_length_vec_1 = [0] * VARIANT_LENGTH_1.output_label_count
    variant_length_vec_2 = [0] * VARIANT_LENGTH_2.output_label_count
    variant_length_vec_1[variant_lengths[0] + VariantLength.index_offset] = 1
    variant_length_vec_2[variant_lengths[1] + VariantLength.index_offset] = 1

    return gt21_vec + genotype_vec + variant_length_vec_1 + variant_length_vec_2
