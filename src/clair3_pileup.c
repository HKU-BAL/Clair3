#define _GNU_SOURCE
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "kvec.h"
#include "medaka_bamiter.h"
#include "medaka_common.h"
#include "clair3_pileup.h"
#include "medaka_khcounter.h"

#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_seqi(s, i) (bam_seqi((s), (i)))
#define bam_nt16_rev_table seq_nt16_str
#define bam_nt16_table seq_nt16_table


/** Constructs a pileup data structure.
 *
 *  @param n_cols number of pileup columns.
 *  @param buffer_cols number of pileup columns.
 *  @param feature_length length of feature vector.
 *  @param num_dtypes number of datatypes in pileup.
 *  @param num_homop maximum homopolymer length to consider.
 *  @param fixed_size if not zero data matrix is allocated as fixed_size * n_cols, ignoring other arguments
 *  @see destroy_plp_data
 *  @returns a plp_data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data create_plp_data(size_t n_cols, size_t buffer_cols, size_t feature_length, size_t num_dtypes, size_t num_homop, size_t fixed_size) {
    assert(buffer_cols >= n_cols);
    plp_data data = xalloc(1, sizeof(_plp_data), "plp_data");
    data->buffer_cols = buffer_cols;
    data->num_dtypes = num_dtypes;
    data->num_homop = num_homop;
    data->n_cols = n_cols;
    if (fixed_size != 0) {
        assert(buffer_cols == n_cols);
        data->matrix = xalloc(fixed_size * n_cols, sizeof(int), "matrix");
    } else {
        data->matrix = xalloc(feature_length * num_dtypes * buffer_cols * num_homop, sizeof(size_t), "matrix");
    }
    data->major = xalloc(buffer_cols, sizeof(size_t), "major");
    data->minor = xalloc(buffer_cols, sizeof(size_t), "minor");
    data->all_alt_info = NULL;
    data->pos_ref_count = NULL;
    data->pos_total_count = NULL;
    return data;
}


/** Enlarge the internal buffers of a pileup data structure.
 *
 *  @param pileup a plp_data pointer.
 *  @param buffer_cols number of pileup columns for which to allocate memory
 *
 */
void enlarge_plp_data(plp_data pileup, size_t buffer_cols, size_t feature_length) {
    assert(buffer_cols > pileup->buffer_cols);
    size_t old_size = feature_length * pileup->num_dtypes * pileup->num_homop * pileup->buffer_cols;
    size_t new_size = feature_length * pileup->num_dtypes * pileup->num_homop * buffer_cols;

    pileup->matrix = xrealloc(pileup->matrix, new_size * sizeof(size_t), "matrix");
    pileup->major = xrealloc(pileup->major, buffer_cols * sizeof(size_t), "major");
    pileup->minor = xrealloc(pileup->minor, buffer_cols * sizeof(size_t), "minor");
    // zero out new part of matrix
    for (size_t i = old_size; i < new_size; ++i) {
        pileup->matrix[i] = 0;
    }
    pileup->buffer_cols = buffer_cols;
}


/** Destroys a pileup data structure.
 *
 *  @param data the object to cleanup.
 *  @returns void.
 *
 */
void destroy_plp_data(plp_data data, bool gvcf) {
    free(data->matrix);
    free(data->major);
    free(data->minor);
    for (size_t i = 0; i < data->candidates_num; i++) {
       free(data->all_alt_info[i]);
    }
    if (gvcf == true) {
        free(data->pos_ref_count);
        free(data->pos_total_count);
    }

    free(data->all_alt_info);
    free(data);
}

/** Generates clair3-style pileup feature data in a region of a bam.
 *
 *  @param region 1-based region string.
 *  @param bam_file input aligment file.
 *  @param tag_value by which to filter data.
 *  @param keep_missing alignments which do not have tag.
 *  @param weibull_summation use predefined bam tags to perform homopolymer partial counts.
 *  @returns a pileup data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 *  If num_dtypes is 1, dtypes should be NULL; all reads in the bam will be
 *  treated equally. If num_dtypes is not 1, dtypes should be an array of
 *  strings, these strings being prefixes of query names of reads within the
 *  bam file. Any read not matching the prefixes will cause exit(1).
 *
 *  If tag_name is not NULL alignments are filtered by the (integer) tag value.
 *  When tag_name is given the behaviour for alignments without the tag is
 *  determined by keep_missing.
 *
 */

/**
 * The pileup input is 594 integers – 33 genome positions wide with 18 features at each position –
 *
 * A+, C+, G+, T+, I_S+, I^1 S+, D_S+, D^1_S+, D_R+, A-, C-, G-, T-, I_S-, I^1_S-, D_S-, D^1_S-, and D_R-
 *
 * A, C, G, T, I, D, +, - means the count of read support of the four nucleotides: insertion,
 * deletion, positive strand, and negative strand. Superscript “1” means only the indel with the
 * highest read support is counted (i.e., all indels are counted if without “1“). Subscript “S”/“R” means
 * the starting/non-starting position of an indel. For example, a 3bp deletion with the most reads support
 * will have the first deleted base counted in either D1_S+ or D1_S-, and the second and third deleted bases
 * counted in either D_R+ or D_R-. The design was determined experimentally, but the rationale is that for
 * 1bp indels that are easy to call, look into the differences between the “S” counts, but reduce the
 * quality if the “R” counts and discrepancy between positions increase.
 *
 */
plp_data calculate_clair3_pileup(const char *region, const bam_fset* bam_set, const char * fasta_path, size_t min_depth, float min_snp_af, float min_indel_af, size_t min_mq, size_t max_indel_length, bool call_snp_only, size_t max_depth, bool gvcf) {
    // extract `chr`:`start`-`end` from `region`
    //   (start is one-based and end-inclusive),
    //   hts_parse_reg below sets return value to point
    //   at ":", copy the input then set ":" to null terminator
    //   to get `chr`.
    int start, end;
    char *chr = xalloc(strlen(region) + 1, sizeof(char), "chr");
    strcpy(chr, region);
    char *reg_chr = (char *) hts_parse_reg(chr, &start, &end);
    // start and end now zero-based end exclusive
    if (reg_chr) {
        *reg_chr = '\0';
    } else {
        fprintf(stderr, "Failed to parse region: '%s'.\n", region);
    }

    // open bam etc.
    // this is all now deferred to the caller
    htsFile *fp = bam_set->fp;
    hts_idx_t *idx = bam_set->idx;
    sam_hdr_t *hdr = bam_set->hdr;
    // setup bam interator

    mplp_data *data = xalloc(1, sizeof(mplp_data), "pileup init data");
    data->fp = fp; data->hdr = hdr; data->iter = bam_itr_querys(idx, hdr, region);
    data->min_mapQ = min_mq;

    bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void **)& data);

    const bam_pileup1_t **plp = xalloc(1, sizeof(bam_pileup1_t *), "pileup");
    int ret, pos, tid, n_plp;

    int n_cols = 0;
    size_t buffer_cols = end - start;
    plp_data pileup = create_plp_data(n_cols, buffer_cols, featlenclair3, 1, 1, 0);

    // get counts
    size_t major_col = 0;  // index into `pileup` corresponding to pos
    n_cols = 0;            // number of processed columns (including insertions, which clair3 doesn't have ;))

    faidx_t* fai = fai_load(fasta_path);
    int len = 0;
    char *ref_seq = NULL;
//    printf("pos: %s %i %i\n", chr, start, end);
    size_t ref_start = max(0, start - mpileup_expand_reference_region);
    size_t ref_end = max(0, end + mpileup_expand_reference_region);
    ref_seq = faidx_fetch_seq(fai, chr, ref_start, ref_end, &len);

    size_t candidates_num = 0;
    size_t alt_info_p_size = 512;
    char ** alt_info_p = xalloc(alt_info_p_size, sizeof(char*), "alt_info_p");
    for (size_t i = 0; i < alt_info_p_size; i++)
        alt_info_p[i] = NULL;

    size_t pre_pos = 0;
    size_t contiguous_flanking_num = 0;

    if (gvcf == true) {
        pileup->pos_ref_count = xalloc(buffer_cols, sizeof(size_t), "pos_ref_count");
        pileup->pos_total_count = xalloc(buffer_cols, sizeof(size_t), "pos_total_count");
        memset(pileup->pos_ref_count, 0, buffer_cols * sizeof(size_t));
        memset(pileup->pos_total_count, 0, buffer_cols * sizeof(size_t));
    }

    while ((ret=bam_mplp_auto(mplp, &tid, &pos, &n_plp, plp) > 0)) {

        size_t depth = 0;
        size_t alt_count = 0;
        size_t ref_count = 0;
        size_t del_count = 0;
        size_t ins_count = 0;

        bool pass_af = false;
        bool pass_snp_af = false;
        bool pass_indel_af = false;

        const char *c_name = data->hdr->target_name[tid];
        if (strcmp(c_name, chr) != 0) continue;
        if (pos < start) continue;
        if (pos >= end) break;
        n_cols++;


        if (pre_pos + 1 != pos || pre_pos == 0)
            contiguous_flanking_num = 0;
        else
            contiguous_flanking_num++;
        pre_pos = pos;

        //update the deletion buffer in each interation
        size_t del_buf_size = 32;
        size_t* dels_f = xalloc(del_buf_size, sizeof(size_t), "dels_f");
        size_t* dels_r = xalloc(del_buf_size, sizeof(size_t), "dels_r");

        memset(dels_f, 0, del_buf_size * sizeof(size_t));
        memset(dels_r, 0, del_buf_size * sizeof(size_t));

        // we still need this as positions might not be contiguous
        pileup->major[major_col / featlenclair3] = pos;
        pileup->minor[major_col / featlenclair3] = 0;

        // counters for insertion strings
        khash_t(KH_COUNTER) *ins_counts_f = kh_init(KH_COUNTER);
        khash_t(KH_COUNTER) *ins_counts_r = kh_init(KH_COUNTER);
        khash_t(KH_COUNTER) *ins_counts_all = kh_init(KH_COUNTER);
        // loop through all reads at this position
        for (int i = 0; i < n_plp; ++i) {
            const bam_pileup1_t *p = plp[0] + i;
            if (p->is_refskip) continue;

            if (p->indel < 0) {
                // there's a deletion starting on next genomic position,
                // record the length here and finalise after the read loop
                //  - actually deleted bases get recorded in next block
                size_t d = (size_t) -1 * p->indel;

                if (d >= del_buf_size) {
                    size_t new_size = max(d, 2 * del_buf_size);
                    dels_f = xrealloc(dels_f, new_size*sizeof(size_t), "dels_f");
                    memset(dels_f+del_buf_size, 0, (new_size-del_buf_size) * sizeof(size_t));
                    dels_r = xrealloc(dels_r, new_size*sizeof(size_t), "dels_r");
                    memset(dels_r+del_buf_size, 0, (new_size-del_buf_size) * sizeof(size_t));
                    del_buf_size = new_size;
                }
                if (bam_is_rev(p->b)) {
                    dels_r[d - 1] += 1;
                } else {
                    dels_f[d - 1] += 1;
                }
            }

            // handle ref_base/sub/del
            int base_i;
            if (p->is_del) {
                // we've been deleted, +1 to DR
                base_i = bam_is_rev(p->b) ? c3_rev_del : c3_fwd_del;
                depth++;
            } else {
                // just a base
                int base_j = bam1_seqi(bam1_seq(p->b), p->qpos);
                if bam_is_rev(p->b) { base_j += 16; }
                base_i = num2countbaseclair3[base_j];
                depth++;
            }
            pileup->matrix[major_col + base_i] += 1;

            // handle insertion
            //  - build insert string then hash
            if (p->indel > 0) {
                size_t first = p->is_del ? 0 : 1;
                char* indel = (char*) xalloc(p->indel + 1, sizeof(char), "indel");
                for (size_t i = 0, j = first; j < p->indel + first; ++i, ++j) {
                    indel[i] = seq_nt16_str[bam1_seqi(bam1_seq(p->b), p->qpos + j)];
                }
                indel[p->indel] = '\0';
                if (bam_is_rev(p->b)) {
                    kh_counter_increment(ins_counts_r, indel);
                } else {
                    kh_counter_increment(ins_counts_f, indel);
                }
                kh_counter_increment(ins_counts_all, indel);
                free(indel);
            }
        }

        // finalise deletions: DS (all) and D1S (best)
        //
        // forward
        size_t best_count = 0;
        size_t all_count = 0;
        for (size_t i = 0; i < del_buf_size; ++i) {
            size_t d = dels_f[i];
            all_count += d;
            best_count = max(best_count, d);
        }
        pileup->matrix[major_col + c3_fwd_del_all] = all_count;
        pileup->matrix[major_col + c3_fwd_del_best] = best_count;
        del_count += all_count;
        // reverse
        best_count = 0;
        all_count = 0;
        for (size_t i = 0; i < del_buf_size; ++i) {
            size_t d = dels_r[i];
            all_count += d;
            best_count = max(best_count, d);
        }
        pileup->matrix[major_col + c3_rev_del_all] = all_count;
        pileup->matrix[major_col + c3_rev_del_best] = best_count;
        del_count += all_count;

        // finalise IS and I1S
        // forward
        kh_counter_stats_t stats = kh_counter_stats(ins_counts_f);
        pileup->matrix[major_col + c3_fwd_ins_all] = stats.sum;
        pileup->matrix[major_col + c3_fwd_ins_best] = stats.max;
        ins_count += stats.sum;

        kh_counter_destroy(ins_counts_f);
        // reverse
        stats = kh_counter_stats(ins_counts_r);
        pileup->matrix[major_col + c3_rev_ins_all] = stats.sum;
        pileup->matrix[major_col + c3_rev_ins_best] = stats.max;
        ins_count += stats.sum;

        kh_counter_destroy(ins_counts_r);
        int offset = pos - ref_start;
        char ref_base = upper_base(ref_seq[offset]);
        int ref_offset_forward = base2index[ref_base - 'A'];
        int ref_offset_reverse = ref_offset_forward + reverse_pos_start;
        char major_alt_base = '\0';
        size_t forward_sum = 0;
        size_t reverse_sum = 0;
        size_t all_alt_count = 0;
        for (size_t i = 0; i < 4; i++) {
            forward_sum += pileup->matrix[major_col + i];
            reverse_sum += pileup->matrix[major_col + i + reverse_pos_start];
            if (i == ref_offset_forward) {
                ref_count = pileup->matrix[major_col + i] + pileup->matrix[major_col + i + reverse_pos_start];
            } else {
                size_t current_count = pileup->matrix[major_col + i] + pileup->matrix[major_col + i + reverse_pos_start];
                if (current_count > alt_count) {
                    alt_count = current_count;
                    major_alt_base = plp_bases_clair3[i];
                    all_alt_count += alt_count;
                }
            }
        }

        pileup->matrix[major_col + ref_offset_forward] = -1 * forward_sum;
        pileup->matrix[major_col + ref_offset_reverse] = -1 * reverse_sum;

        // calculate candidate allele frequency and apply filtering
        depth = max(1, depth);
        bool pass_min_depth = depth >= min_depth;
        bool pass_ref_base_in_acgt = ref_base == 'A' || ref_base == 'C' || ref_base == 'G' || ref_base == 'T';
        bool non_ref_base_majority = ref_count < alt_count || ref_count < ins_count || ref_count < del_count;
        bool ref_alt_equal_majority = (ref_count > 0 && ref_count == alt_count && ref_base - major_alt_base < 0);
        if (call_snp_only == true) {
            pass_af = alt_count / (float)depth >= min_snp_af;
        } else {
            pass_af = non_ref_base_majority || ref_alt_equal_majority || (alt_count / (float)depth >= min_snp_af);
            pass_af = pass_af || (del_count / (float)depth >= min_indel_af) || (ins_count / (float)depth >= min_indel_af);
        }

        pass_af = pass_af && pass_min_depth && pass_ref_base_in_acgt;
        pass_af = pass_af && (contiguous_flanking_num >= pileup_flanking_base_num);
        // move to next position
        if (pass_af) {

            if (candidates_num + 1 >= alt_info_p_size) {
                alt_info_p_size = alt_info_p_size << 1;
                alt_info_p = xrealloc(alt_info_p, alt_info_p_size * sizeof(char*), "alt_info_p");
            }

            size_t max_alt_length = 64;
            char *alt_info_str = xalloc(max_alt_length, sizeof(char), "alt_info_str");

            sprintf(alt_info_str, "%i-%i-%c-", pos+1, depth, ref_base);
            //snp
            for (size_t i = 0; i < 4; i++) {
                forward_sum += pileup->matrix[major_col + i];
                reverse_sum += pileup->matrix[major_col + i + reverse_pos_start];
                size_t alt_sum = pileup->matrix[major_col + i] + pileup->matrix[major_col + i + reverse_pos_start];

                if (alt_sum > 0 && i != ref_offset_forward)
                    sprintf(alt_info_str + strlen(alt_info_str), "X%c %i ", plp_bases_clair3[i], alt_sum);
            }
            //del
            for (size_t i = 0; i < del_buf_size; i++) {
                size_t d = dels_f[i] + dels_r[i];
                if (d > 0 && i+1 <= max_indel_length) {
                    // 32 bytes is a safe number for integer to string
                    if (strlen(alt_info_str) + i + 32 >= max_alt_length) {
                        while (strlen(alt_info_str) + i + 32 >= max_alt_length)
                            max_alt_length = max_alt_length << 1;
                         alt_info_str = xrealloc(alt_info_str, max_alt_length*sizeof(char), "alt_info_str");
                    }
                    sprintf(alt_info_str + strlen(alt_info_str), "D%.*s %i ", i+1,ref_seq+offset+1, d);
                }

            }
//            //ins
            for (khiter_t k = kh_begin(ins_counts_all); k != kh_end(ins_counts_all); ++k) {
                if (kh_exist(ins_counts_all, k)) {
                    const char *key = kh_key(ins_counts_all, k);
                    size_t val = kh_val(ins_counts_all, k);
                    if (strlen(key) <= max_indel_length) {
                         if (strlen(alt_info_str) + strlen(key) + 32 >= max_alt_length) {
                             while (strlen(alt_info_str) + strlen(key) + 32 >= max_alt_length)
                                 max_alt_length = max_alt_length << 1;
                             alt_info_str = xrealloc(alt_info_str, max_alt_length *sizeof(char), "alt_info_str");
                        }
                        sprintf(alt_info_str + strlen(alt_info_str), "I%c%s %i ", ref_base, key, val);
                    }
                }
            }
            // update the alternative information for current candidates here
            alt_info_p[candidates_num++] = alt_info_str;
        }

        if (gvcf == true) {
            pileup->pos_ref_count[pos-start] = ref_count;
            pileup->pos_total_count[pos-start] = ref_count + all_alt_count + del_count + ins_count;
        }

        free(dels_f);
        free(dels_r);
        kh_counter_destroy(ins_counts_all);
        major_col += featlenclair3;
    }


    pileup->all_alt_info = alt_info_p;
    pileup->candidates_num = candidates_num;
    pileup->n_cols = n_cols;

    bam_itr_destroy(data->iter);
    bam_mplp_destroy(mplp);
    fai_destroy(fai);
    free(data);
    free(plp);
    free(chr);

    return pileup;
}

int main()
{
    return 0;
}
