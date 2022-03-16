#ifndef _CLAIR3_PILEUP_H
#define _CLAIR3_PILEUP_H

// medaka-style feature data
typedef struct _plp_data {
    size_t buffer_cols;
    size_t num_dtypes;
    size_t num_homop;
    size_t n_cols;
    size_t *matrix;
    size_t *major;
    size_t *minor;
    char **all_alt_info;
    size_t candidates_num;
    size_t* pos_ref_count;
    size_t* pos_total_count;
} _plp_data;
typedef _plp_data *plp_data;


// convert 16bit IUPAC (+16 for strand) to plp_bases index
// {
//  ,  A,  C,   ,  G,   ,   ,   , 
// T,   ,   ,   ,   ,   ,   ,   ,
//  ,  a,  c,   ,  g,   ,   ,   ,
// t,  ,    ,   ,   ,   ,   ,   ,
// }
static const int num2countbase[32] = {
  -1,  4,  5, -1,  6, -1, -1, -1,
   7, -1, -1, -1, -1, -1, -1, -1,
  -1,  0,  1, -1,  2, -1, -1, -1,
   3, -1, -1, -1, -1, -1, -1, -1,
};


static const int base2index[32] = {
    0, 0, 1, 0, 0, 0, 2, 0, // abcdefgh
    0, 0, 0, 0, 0, 0, 0, 0, // ijklmnop
    0, 0, 0, 3, 0, 0, 0, 0, // qrstuvwx
    0, 0, 0, 0, 0, 0, 0, 0, // vz
};


// convert 16bit IUPAC (+16 for strand) to plp_bases clair3 index
//  first i: all insertions
// second i: most common insertion
//  first d: all first base deletion  (actually a reference base)
// second d: most common deletion     (actually a reference base)
//  third d: non-first base deletion  (the deleted bases)
static const char plp_bases_clair3[] = "ACGTIIDDDacgtiiddd";
static const size_t featlenclair3 = 18;   // len of the above
static const size_t c3_fwd_ins_all = 4;     
static const size_t c3_fwd_ins_best = 5;
static const size_t c3_fwd_del_all = 6;   // (preceding ref position)
static const size_t c3_fwd_del_best = 7;  // (preceding ref position)
static const size_t c3_fwd_del = 8;       // (actually deleted base)
static const size_t c3_rev_ins_all = 13;     
static const size_t c3_rev_ins_best = 14;
static const size_t c3_rev_del_all = 15;  // (preceding ref position)
static const size_t c3_rev_del_best = 16; // (preceding ref position)
static const size_t c3_rev_del = 17;      // (actually deleted base)
static const size_t reverse_pos_start = 9;  // position of reverse position start
static const size_t mpileup_expand_reference_region = 1000;
static const size_t pileup_flanking_base_num = 16;

static const int num2countbaseclair3[32] = {
 -1,  0,  1, -1,  2, -1, -1, -1,
  3, -1, -1, -1, -1, -1, -1, -1,
 -1,  9, 10, -1, 11, -1, -1, -1,
 12, -1, -1, -1, -1, -1, -1, -1,
};


/** Constructs a pileup data structure.
 *
 *  @param n_cols number of pileup columns.
 *  @param buffer_cols number of pileup columns.
 *  @param num_dtypes number of datatypes in pileup.
 *  @param num_homop maximum homopolymer length to consider.
 *  @param fixed_size if not zero data matrix is allocated as fixed_size * n_cols, ignoring other arguments
 *  @see destroy_plp_data
 *  @returns a plp_data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data create_plp_data(size_t n_cols, size_t buffer_cols, size_t feature_length, size_t num_dtypes, size_t num_homop, size_t fixed_size);


/** Destroys a pileup data structure.
 *
 *  @param data the object to cleanup.
 *  @returns void.
 *
 */
void destroy_plp_data(plp_data data, bool gvcf);

/** C implement of clair3-style pileup feature data and alternative information in a given region of a bam.
 *
 *  @param region  1-based region string
 *  @param bam_set  bam handler of input bam
 *  @param fasta_path  input reference file
 *  @param min_depth  minimum coverage required to call a variant
 *  @param min_snp_af  minimum snp allele frequency for a site to be considered as a candidate site
 *  @param min_indel_af  minimum indel allele frequency for a site to be considered as a candidate site
 *  @param min_mq  minimum mapping quality for read to use for calling
 *  @param max_indel_length  maximum indel length to format into alternative string stream
 *  @returns a pileup data pointer, including the data matrix and all candidates alternative information
 *
 *  The return value can be freed with destroy_plp_data
 *
 */
plp_data calculate_clair3_pileup(const char *region, const bam_fset* bam_set, const char * fasta_path, size_t min_depth, float min_snp_af, float min_indel_af, size_t min_mq, size_t max_indel_length, bool call_snp_only, size_t max_depth, bool gvcf);

#endif
