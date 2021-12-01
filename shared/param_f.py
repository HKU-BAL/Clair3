# Clair3 full alignment parameters
REPO_NAME = "Clair3"
from itertools import accumulate

zstd='zstd'
default_optimizer = "Radam"
default_loss_function = "FocalLoss"
min_af = 0.08
min_af_dict = {'ont':0.15, 'hifi':min_af, 'ilmn':min_af }
matrix_depth_dict = {'ont': 89, 'hifi': 55, 'ilmn': 55}
max_depth = 144
maximum_variant_length_that_need_infer = 50
maximum_variant_length_that_need_infer_include_long_indel = 100000
cal_precise_long_indel_af = False
long_indel_distance_proportion = 0.1
min_mq = 5
min_bq = 0
min_coverage = 2

# Full alignment input feature list
channel = (
'reference_base', 'alternative_base', 'mapping_quality', 'base_quality', 'strand_info', 'variant_type', 'insert_base',
'phasing_info')  # phasing info if add_phasing
channel_size = len(channel)
flankingBaseNum = 16
no_of_positions = 2 * flankingBaseNum + 1
input_shape = [matrix_depth_dict['hifi'], no_of_positions, channel_size]
ont_input_shape = [matrix_depth_dict['ont'], no_of_positions, channel_size]
label_shape = [21, 3, no_of_positions, no_of_positions]
label_size = sum(label_shape)
label_shape_cum = list(accumulate(label_shape))
expandReferenceRegion = 1000
SAMTOOLS_VIEW_FILTER_FLAG = 2316
NORMALIZE_NUM = 100

# Realignment parameters
partition_size = 500000
realign_chunk_size = 5000
phasing_window_size = 100000
illumina_phasing_window_size = 10000
max_phasing_depth = 15
min_phasing_read_coverage = 2
split_region_size = 1000
extend_bp = 10

# Training hyperparameters
chunk_size = 200
trainBatchSize = 2000
predictBatchSize = 200
initialLearningRate = 1e-3
l2RegularizationLambda = 1e-7
trainingDatasetPercentage = 0.9
maxEpoch = 30
OPERATION_SEED = None
RANDOM_SEED = None
