from collections import namedtuple

VariantLengthNamedTuple = namedtuple(
    'VariantLengthNamedTuple', ['index_offset', 'min', 'max', 'output_label_count']
)
variant_length_index_offset = 16
VariantLength = VariantLengthNamedTuple(
    index_offset=variant_length_index_offset,
    min=-variant_length_index_offset,
    max=variant_length_index_offset,
    output_label_count=variant_length_index_offset * 2 + 1,
)
