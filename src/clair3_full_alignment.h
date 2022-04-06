#ifndef _CLAIR3_FULL_ALIGNMENT_H
#define _CLAIR3_FULL_ALIGNMENT_H

#define HAP_UNPHASED 0
#define HAP_1 1
#define HAP_2 2

#define normalize_mq(x) ((int)(x < 60 ? 100 * x / 60.0 : 100))
#define normalize_bq(x) ((int)(x < 40 ? 100 * x / 40.0 : 100))
#define normalize_af(x) ((int)(x < 1.0 ? 100 * x : 100))
#define normalize_strand(x) (x == true ? 50 : 100)

static const int8_t HAP_TYPE[3] = {60, 30, 90};
#define normalize_hap(x) (HAP_TYPE[x])

static const size_t overhang = 10;
static const char *RN = "\0";
static const size_t min_haplotag_mq = 20;
static const size_t expand_reference_region = 2000000;
static const size_t flanking_base_num = 16;
static const size_t no_of_positions = 33;
static const size_t channel_size = 8;
static const size_t min_bq = 0;
static const size_t SAMTOOLS_VIEW_FILTER_FLAG = 2316;
static const size_t MAX_READ_COUNT = 1000;
static const size_t MAX_INDEL_LENGTH = 50;
static const char ACGT[] = "ACGT";

// convert 16bit IUPAC (+16 for strand) to plp_bases index
// {
//  ,  A,  C,   ,  G,   ,   ,   ,
// T,   ,   ,   ,   ,   ,   ,   ,
//  ,  a,  c,   ,  g,   ,   ,   ,
// t,  ,    ,   ,   ,   ,   ,   ,
// }
static const int8_t num2countbase_fa[32] = {
    100, 0, 25, -100, 0, 0, 75, 0, // abcdefgh
    -50, 0, 0, 0, 0, 100, 0, 0,    // ijklmnop
    0, 0, 0, 50, 0, 0, 0, 0,       // qrstuvwx
    0, 0, 0, 0, 0, 0, 0, 0,        // vz
};

// convert A-Z character to 0-index offset
// ACGT: 0123
// non-ACGT: 0
static const int8_t acgt2num[32] = {
    0, 0, 1, 0, 0, 0, 2, 0, // abcdefgh
    0, 0, 0, 0, 0, 0, 0, 0, // ijklmnop
    0, 0, 0, 3, 0, 0, 0, 0, // qrstuvwx
    0, 0, 0, 0, 0, 0, 0, 0, // vz
};

/*! @typedef
 @abstract Structure for full-alignment data
 @field matrix  int array of (total candidate number * matrix depth * no of flanking position * feature channel)
 @field alt_info  alternative information string with all candidates, including all SNPs, insertions and deletions
 @field alt_info_length  length of the alternative information string
 */
typedef struct _fa_data
{

    int8_t *matrix;
    char **all_alt_info;
    size_t candidates_num;
} _fa_data;

typedef _fa_data *fa_data;

/*! @typedef
 @abstract Structure for matrix level alternative information
 @field ins_bases  the char string storing all insertion bases in current position of an alignment
 @field alt_base  alternative base other than reference base in query sequence
 @field del_length  deletion length in current position of an alignment
 @field has_alt_info  true if any of alternative information exists, false for reference base and deletion bases(#*)

 @ by default we only allocate a maximum `matrix depth` struct array and reset all field in each candidate iteration, we
   need to calculate each read candidate proportion in given candidate overlapped region
 */
typedef struct Alt_info
{
    char *ins_bases;
    char alt_base;
    size_t del_length;
    bool has_alt_info;
} Alt_info;

/*! @typedef
 @abstract Structure for matrix level alternative information
 @field read_index  the read start offset of each read, the index is sorted by read start
 @field haplotype information of read, 0: unphased or not phasable 1|2: haplotype1|2
 */
typedef struct HAP
{
    size_t read_index;
    size_t haplotype;
} HAP;

/*! @typedef
 @abstract Structure of a phased heterozygous pileup SNP variant
 @field position  variant start position 0-index
 @field ref_base  reference base tag in VCF
 @field alt_base  alternative base tag in VCF
 @field genotype  phased heterozygous genotype, 0|1 : 1,  1|0: 2
 @field phase_set phase set tag in VCF, which is acquired from whatshap or longphase

 @ in this release, we only store heterozygous SNP info
 */
typedef struct Variant
{
    int position;
    char ref_base;
    char alt_base;
    int genotype;
    int phase_set;
} Variant;

typedef struct Variants_info
{
    Variant **variants;
    size_t variant_num;
    size_t variant_current_pos;
} Variants_info;

/*! @typedef
 @abstract Structure for matrix level alternative information
 @field ins_bases  the char string storing all insertion bases in current position of an alignment
 @field ins_length  length the stored insertion bases
 @field alt_base  alternative base in htslib int format
 @field del_length  deletion length in current position of an alignment
 @field bq phred quality score of given bases

 @ we use the htslib format int alt_base than char as we need to mark the '#*' into -1, for bq field, we only store
   reference base and alternative base quality and skip the insertion quality as there are only one base quality channel
 */
typedef struct Pos_info
{
    char *ins_bases;
    size_t ins_length;
    int alt_base;
    size_t del_length;
    int8_t bq;
} Pos_info;

/*! @typedef
 @abstract Structure for the alignment information
 @field read_start  read start position of alignment, 0-index
 @field q_name  read name
 @field read_end  alignment read end compared with the reference sequence, CIGAR length sum of X=MDN
 @field cigartuples  alignment CIGAR int pointer from htslib bam_get_cigar function
 @field qual  base quality int pointer from htslib core alignment
 @field mq  normalized mapping quality value (0-100)
 @field n_cigar  number of CIGAR operations
 @field l_qseq  length of the read query sequence
 @field haplotype  haplotype information of read, 0: unphased or not phasable 1|2: haplotype1|2
 @field strand  normalized strand value forward: 50 reverse: 100
 @field pos_info  structure array of overlapped flanking candidates information
 @field overlap_candidates_num  number of overlapped flanking candidates between read start and read end, including flanking bases
 @field flanking_start  the first overlapped candidate index (0 index is the the first candidate - 16 by default)

 @note that the seqi and qual pointer information will be released after each htslib sam_itr_next iterator
 */
typedef struct Read
{
    size_t read_start;
    char *q_name;
    size_t read_end;
    uint32_t *cigartuples;
    uint8_t *seqi;
    uint8_t *qual;
    int8_t mq;
    size_t n_cigar;
    uint32_t l_qseq;
    size_t haplotype;
    int8_t strand;
    Pos_info *pos_info;
    size_t overlap_candidates_num;
    size_t flanking_start;
} Read;

/** Destroys a full-alignment data structure.
 *
 *  @param data the full-alignment data object to cleanup.
 *  @returns void.
 *
 */
void destroy_fa_data(fa_data data);

/** Sort overlapped reads of a candidate based on hapltoype information and read start
 *
 *  @param read_hap_array  struct array of all overlap reads
 *  @param matrix_read_index_array  the return reference of the read index array, -1 for padding
 *  @param n  number of overlapped reads
 *  @returns void.
 *
 */
void sort_read_name_by_haplotype(HAP *read_hap_array, int *matrix_read_index_array, size_t matrix_depth, size_t n);

/** get all overlapped flanking candidates number and start position based on read start and read end
 *
 *  @param read_start  read start, 0-index
 *  @param read_end  read end, 0-index
 *  @param candidate_current_index  the first flanking candidate index >= read start
 *  @param flanking_candidates  int array of all flanking candidates, sorted by start position
 *  @returns number of the total overlapped flanking candidates within read start and read end.
 *
 */
size_t get_overlap_candidate_num(size_t read_start, size_t read_end, size_t candidate_current_index, size_t flanking_candidates_num, size_t *flanking_candidates);

/** get the substring of a reference sequence based on start and end
 *
 *  @param ref_seq  a string store all reference sequence from ref_start(0-index)
 *  @param start  sequence query start, 0-index
 *  @param end  sequence query end, 0-index
 *  @returns string of the queried region of reference sequence
 *
 */
char *get_ref_seq(char *ref_seq, size_t start, size_t end);

/** get the substring of a query sequence based on start and end
 *
 *  @param seqi  a htslib format pointer stores all query sequence(0-index)
 *  @param start  query start, 0-index
 *  @param end  query end, 0-index
 *  @returns string of the queried sequence
 *
 */
char *get_query_seq(uint8_t *seqi, size_t start, size_t end);

/** C implement of whatshap hapltagging
 *
 */
void cigar_prefix_length(uint32_t *cigartuples, size_t reference_bases, size_t *ref_bases, size_t *query_bases, size_t left_cigar_index, size_t right_cigar_index, size_t consumed, bool reverse);

int realign_read(Variant *variant, Read *read, size_t i, size_t consumed, size_t query_pos, char *reference, size_t ref_start);

int haplotag_read(Variants_info *variants_info, Read *read, char *ref_seq, size_t ref_start);

/** C implement of clair3-style full-alignment feature data and alternative information in a given region of a bam.
 *
 *  @param region  1-based region string
 *  @param bam_path  input alignment file
 *  @param fasta_path  input reference file
 *  @param variants  C structure pointer of all phased heterozygous pileup SNP variants
 *  @param variant_num  total variants number
 *  @param candidates  int array of all low-quality pileup candidates need to process (0-index)
 *  @param candidate_num total candidates number
 *  @returns a full-alignment data pointer, including the data matrix and all candidates alternative information
 *
 *  The return value can be freed with destroy_fa_data
 *
 */
fa_data calculate_clair3_full_alignment(const char *region, const char *bam_path, const char *fasta_path, Variant **variants, size_t variant_num, size_t *candidates, size_t candidate_num, bool need_haplotagging, size_t min_mq, size_t min_bq, size_t matrix_depth, size_t max_indel_length);

#endif
