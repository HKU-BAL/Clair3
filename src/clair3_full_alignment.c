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
#include "khash.h"
#include "kvec.h"
#include "medaka_bamiter.h"
#include "medaka_common.h"
#include "medaka_khcounter.h"
#include "clair3_full_alignment.h"
#include "levenshtein.h"

typedef struct Pos_alt_info
{

    khash_t(KH_COUNTER) * ins_counter;
    khash_t(KH_INT_COUNTER) * del_counter;
    size_t acgt_count[4];
    size_t depth;

} Pos_alt_info;

int com_func(const void *a, const void *b)
{

    return (*(size_t *)a - *(size_t *)b);
}

int hap_cmp(const void *x, const void *y)
{

    HAP a = *(HAP *)x;
    HAP b = *(HAP *)y;
    if (a.haplotype < b.haplotype)
        return -1;
    else if (a.haplotype > b.haplotype)
        return 1;
    else
        return (a.read_index - b.read_index);
}

void destroy_fa_data(fa_data data)
{

    free(data->matrix);
    for (size_t i = 0; i < data->candidates_num; i++)
    {
       free(data->all_alt_info[i]);
    }
    free(data->all_alt_info);
    free(data);
}

void sort_read_name_by_haplotype(HAP *read_hap_array, int *matrix_read_index_array, size_t matrix_depth, size_t n)
{

    size_t read_num = min(n, matrix_depth);
    if (n > matrix_depth)
    {
        // shuffle the read index array with the same random seed
        for (size_t i = 0; i < n - 1; i++)
        {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
            size_t tmp_read_index = read_hap_array[j].read_index;
            size_t tmp_haplotype = read_hap_array[j].haplotype;
            read_hap_array[j].read_index = read_hap_array[i].read_index;
            read_hap_array[j].haplotype = read_hap_array[i].haplotype;
            read_hap_array[i].read_index = tmp_read_index;
            read_hap_array[i].haplotype = tmp_haplotype;
        }
    }

    qsort(read_hap_array, read_num, sizeof(HAP), hap_cmp);

    // if need padding (overlap read num < matrix depth), add -1 to fill the depth
    if (n < matrix_depth)
    {
        size_t padding_depth = matrix_depth - read_num;
        size_t prefix_padding_depth = padding_depth >> 1;
        size_t suffix_padding_depth = padding_depth - prefix_padding_depth;
        for (size_t i = 0; i < prefix_padding_depth; i++)
            matrix_read_index_array[i] = -1;
        for (size_t i = 0; i < read_num; i++)
            matrix_read_index_array[i + prefix_padding_depth] = read_hap_array[i].read_index;
        for (size_t i = 0; i < suffix_padding_depth; i++)
            matrix_read_index_array[read_num + prefix_padding_depth + i] = -1;
    }
    else
    {
        for (size_t i = 0; i < matrix_depth; i++)
            matrix_read_index_array[i] = read_hap_array[i].read_index;
    }
}

void cigar_prefix_length(uint32_t *cigartuples, size_t reference_bases, size_t *ref_bases, size_t *query_bases, size_t left_cigar_index, size_t right_cigar_index, size_t consumed, bool reverse)
{

    size_t ref_pos = 0;
    size_t query_pos = 0;
    for (size_t i = left_cigar_index; i < right_cigar_index; i++)
    {
        size_t index = reverse ? left_cigar_index + right_cigar_index - i - 1 : i;
        size_t cigar_op = bam_cigar_op(cigartuples[index]);
        size_t length = bam_cigar_oplen(cigartuples[index]);

        length = i == left_cigar_index ? consumed : length;
        if (length == 0)
            continue;

        if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF)
        {
            query_pos += length;
            ref_pos += length;
            if (ref_pos >= reference_bases)
            {
                *ref_bases = reference_bases;
                *query_bases = query_pos + reference_bases - ref_pos;
                return;
            }
        }
        else if (cigar_op == BAM_CDEL)
        {
            ref_pos += length;
            if (ref_pos >= reference_bases)
            {
                *ref_bases = reference_bases;
                *query_bases = query_pos;
                return;
            }
        }
        else if (cigar_op == BAM_CINS)
        {
            query_pos += length;
        }
        else if (cigar_op == BAM_CREF_SKIP)
        {
            *ref_bases = reference_bases;
            *query_bases = query_pos;
            return;
        }
    }
}

char *get_ref_seq(char *ref_seq, size_t start, size_t end)
{

    size_t seq_size = end - start;
    char *sub_seq = malloc((seq_size + 1));
    strncpy(sub_seq, ref_seq + start, seq_size);
    sub_seq[seq_size] = '\0';
    return sub_seq;
}

size_t get_read_end(uint32_t *cigartuples, size_t n_cigar, size_t read_start)
{

    size_t ref_pos = read_start;
    for (size_t i = 0; i < n_cigar; i++)
    {
        size_t cigar_op = bam_cigar_op(cigartuples[i]);
        size_t length = bam_cigar_oplen(cigartuples[i]);
        if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF || cigar_op == BAM_CDEL || cigar_op == BAM_CREF_SKIP)
        {
            ref_pos += length;
        }
    }
    return ref_pos;
}

char *get_query_seq(uint8_t *seqi, size_t start, size_t end)
{

    size_t seq_size = end - start;
    char *sub_seq = malloc((seq_size + 1));
    for (size_t i = 0; i < seq_size; i++)
    {
        sub_seq[i] = seq_nt16_str[bam_seqi(seqi, start + i)];
    }
    sub_seq[seq_size] = '\0';
    return sub_seq;
}

void update_haplotype_cost(int allele, int phase_set, int genotype, khash_t(KH_INT_COUNTER) * haplotype_cost)
{

    if (allele == 0)
        return;

    if (allele == genotype)
    {
        kh_int_counter_add(haplotype_cost, phase_set, 1);
    }
    else
    {
        kh_int_counter_add(haplotype_cost, phase_set, -1);
    }
}

int realign_read(Variant *variant, Read *read, size_t i, size_t consumed, size_t query_pos, char *reference, size_t ref_start)
{

    uint32_t *cigartuples = read->cigartuples;
    uint8_t *seqi = read->seqi;
    size_t n_cigar = read->n_cigar;
    size_t middle_op = bam_cigar_op(cigartuples[i]);
    size_t middle_length = bam_cigar_oplen(cigartuples[i]);
    size_t left_consumed = consumed > 0 ? consumed : 0;
    size_t right_consumed = consumed < middle_length ? middle_length - consumed : 0;
    size_t left_ref_bases = 0;
    size_t left_query_bases = 0;
    size_t right_ref_bases = 0;
    size_t right_query_bases = 0;
    size_t left_cigar_size = i + 1;
    size_t right_cigar_size = i;

    cigar_prefix_length(cigartuples, overhang, &left_ref_bases, &left_query_bases, 0, left_cigar_size, left_consumed, true);
    cigar_prefix_length(cigartuples, overhang + 1, &right_ref_bases, &right_query_bases, right_cigar_size, n_cigar, right_consumed, false);

    char *query = get_query_seq(seqi, query_pos - left_query_bases, query_pos + right_query_bases);
    char *ref = get_ref_seq(reference, variant->position - left_ref_bases - ref_start, variant->position + right_ref_bases - ref_start);

    size_t alt_length = left_ref_bases + right_ref_bases + 1;
    char *alt = malloc(alt_length);
    strcpy(alt, ref);
    alt[left_ref_bases] = variant->alt_base;

    size_t distance_ref = levenshtein(query, ref);
    size_t distance_alt = levenshtein(query, alt);

    int allele = 0;
    if (distance_ref < distance_alt)
    {
        allele = 1;
    }
    else if (distance_ref > distance_alt)
    {
        allele = 2;
    }

    free(query);
    free(ref);
    free(alt);

    return allele;
}

int haplotag_read(Variants_info *variants_info, Read *read, char *ref_seq, size_t ref_start)
{

    size_t n = variants_info->variant_num;
    size_t query_pos = 0;
    size_t v_position = 0;
    Variant **variants = variants_info->variants;
    uint8_t *seqi = read->seqi;
    uint32_t *cigartuples = read->cigartuples;
    size_t n_cigar = read->n_cigar;
    size_t j = variants_info->variant_current_pos;
    size_t ref_pos = read->read_start;
    khash_t(KH_INT_COUNTER) *haplotype_cost = kh_init(KH_INT_COUNTER);
    int allele = 0;

    while (j < n && variants[j]->position < ref_pos)
        j += 1;

    for (size_t i = 0; i < n_cigar; i++)
    {
        size_t cigar_op = bam_cigar_op(cigartuples[i]);
        size_t length = bam_cigar_oplen(cigartuples[i]);
        if (j < n)
            v_position = variants[j]->position;

        if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF)
        { // XM=
            while (j < n && v_position < ref_pos + length)
            {
                allele = realign_read(variants[j], read, i, v_position - ref_pos, query_pos + v_position - ref_pos, ref_seq, ref_start);
                update_haplotype_cost(allele, variants[j]->phase_set, variants[j]->genotype, haplotype_cost);
                j++;
                if (j < n)
                    v_position = variants[j]->position;
            }
            query_pos += length;
            ref_pos += length;
        }
        else if (cigar_op == BAM_CINS)
        { // I
            if (j < n && v_position == ref_pos)
            {
                allele = realign_read(variants[j], read, i, 0, query_pos, ref_seq, ref_start);
                update_haplotype_cost(allele, variants[j]->phase_set, variants[j]->genotype, haplotype_cost);
                j++;
                if (j < n)
                    v_position = variants[j]->position;
            }
            query_pos += length;
        }
        else if (cigar_op == BAM_CDEL)
        {
            while (j < n && v_position < ref_pos + length)
            {
                allele = realign_read(variants[j], read, i, v_position - ref_pos, query_pos, ref_seq, ref_start);
                update_haplotype_cost(allele, variants[j]->phase_set, variants[j]->genotype, haplotype_cost);
                j++;
                if (j < n)
                    v_position = variants[j]->position;
            }
            ref_pos += length;
        }
        else if (cigar_op == BAM_CREF_SKIP)
        {
            while (j < n && v_position < ref_pos + length)
            {
                j++;
                if (j < n)
                    v_position = variants[j]->position;
            }
            ref_pos += length;
        }
        else if (cigar_op == BAM_CSOFT_CLIP)
        {
            query_pos += length;
        }
    }

    read->read_end = ref_pos;

    size_t counter_size = 0;
    int max_value = 0;
    int min_value = 0;
    for (khiter_t k = kh_begin(haplotype_cost); k != kh_end(haplotype_cost); ++k)
    {
        if (kh_exist(haplotype_cost, k))
        {
            int val = kh_val(haplotype_cost, k);
            max_value = max(max_value, val);
            min_value = min(min_value, val);
            counter_size++;
        }
    }

    kh_int_counter_destroy(haplotype_cost);

    if (counter_size == 0 || (max_value == 0 && min_value == 0))
    {
        return HAP_UNPHASED;
    }
    else if (max_value > abs(min_value))
    {
        return HAP_1;
    }
    else
    {
        return HAP_2;
    }
}

size_t get_overlap_candidate_num(size_t read_start, size_t read_end, size_t candidate_current_index, size_t flanking_candidates_num, size_t *flanking_candidates)
{
    size_t overlap_num = 0;
    for (size_t i = candidate_current_index; i < flanking_candidates_num; i++)
    {
        if (flanking_candidates[i] >= read_start && flanking_candidates[i] < read_end)
            overlap_num++;
        else
            return overlap_num;
    }
    return overlap_num;
}

fa_data calculate_clair3_full_alignment(const char *region, const char *bam_path, const char *fasta_path, Variant **variants, size_t variant_num, size_t *candidates, size_t candidate_num, bool need_haplotagging, \
size_t min_mq, size_t min_bq, size_t matrix_depth, size_t max_indel_length)
{

    int start, end;
    char *chr = xalloc(strlen(region) + 1, sizeof(char), "chr");
    strcpy(chr, region);
    char *reg_chr = (char *)hts_parse_reg(chr, &start, &end);
    if (reg_chr)
        *reg_chr = '\0';

    // open fasta
    faidx_t *fai = fai_load(fasta_path);
    int len = 0;
    char *ref_seq = NULL;

    const size_t offset_can = no_of_positions * matrix_depth * channel_size;
    const size_t offset_row = no_of_positions * channel_size;
    const size_t offset_col = channel_size;

    int ref_start = max(0, start - expand_reference_region); // 0-index
    int ref_end = end + expand_reference_region;
    ref_seq = faidx_fetch_seq(fai, chr, ref_start, ref_end, &len);

    // open bam
    htsFile *hts_file;
    hts_idx_t *idx;
    bam_hdr_t *header;

    hts_file = sam_open(bam_path, "r");
    idx = sam_index_load(hts_file, bam_path);
    header = sam_hdr_read(hts_file);
    const int tid = bam_name2id(header, chr);
    hts_itr_t *iter = sam_itr_queryi(idx, tid, start, end);
    // initialize an alignment
    bam1_t *alignment = bam_init1();

    size_t reads_num = 0;
    size_t variant_current_pos = 0;
    size_t flanking_candidates_num = 0;
    size_t candidate_current_index = 0;
    size_t read_no_overlap_num = 0;
    Pos_info *pos_info = NULL;

    Variants_info variants_info = {
        .variants = variants,
        .variant_num = variant_num,
        .variant_current_pos = variant_current_pos};

    // dict to store all candidates index
    khash_t(KH_INT_COUNTER) *candidates_p = kh_init(KH_INT_COUNTER);
    // dict to store all flanking candidate index
    khash_t(KH_INT_COUNTER) *flanking_candidates_p = kh_init(KH_INT_COUNTER);
    // dict to count all read name
    khash_t(KH_COUNTER) *read_name_set = kh_init(KH_COUNTER);
    // allocate a position alternative information struct for each candidate
    Pos_alt_info *pos_alt_info = calloc(candidate_num, sizeof(Pos_alt_info));
    // a kvec vector to store all read struct
    kvec_t(Read) read_array;
    kv_init(read_array);

    for (size_t i = 0; i < candidate_num; i++)
    {
        size_t candidate = candidates[i];
        // each candidate is a new key
        kh_int_counter_add(candidates_p, candidate, i);
        pos_alt_info[i].ins_counter = kh_init(KH_COUNTER);
        pos_alt_info[i].del_counter = kh_init(KH_INT_COUNTER);
        pos_alt_info[i].depth = 0;
        for (size_t j = 0; j < 4; j++)
            pos_alt_info[i].acgt_count[j] = 0;

        for (size_t j = 0; j < no_of_positions; j++)
        {
            size_t key = candidate - flanking_base_num + j;
            if (kh_int_counter_val(flanking_candidates_p, key) == -1)
            {
                kh_int_counter_add(flanking_candidates_p, key, flanking_candidates_num++);
            }
        }
    }

    size_t flanking_candidates[flanking_candidates_num];
    for (khiter_t k = kh_begin(flanking_candidates_p); k != kh_end(flanking_candidates_p); ++k)
    {
        if (kh_exist(flanking_candidates_p, k))
        {
            size_t key = kh_key(flanking_candidates_p, k);
            int val = kh_val(flanking_candidates_p, k);
            flanking_candidates[val] = key;
        }
    }

    while (sam_itr_next(hts_file, iter, alignment) >= 0)
    {
        int flag = alignment->core.flag;

        if (flag & SAMTOOLS_VIEW_FILTER_FLAG)
            continue;

        if (alignment->core.qual < min_mq)
        {
            continue;
        }

        const char *q_name = bam_get_qname(alignment);
        //skip the duplicated read name
        int ret = 0;
        khiter_t k = kh_put(KH_COUNTER, read_name_set, q_name, &ret);
        if (ret == 1)
        {
            kh_key(read_name_set, k) = strdup(q_name);
            kh_value(read_name_set, k) = 1;
        }
        else if (ret == 0)
        {
            continue;
        }

        bool is_fwd_strand = (flag & 16) == 16;
        int32_t pos = alignment->core.pos;
        uint32_t l_qseq = alignment->core.l_qseq;
        uint32_t *cigartuples = bam_get_cigar(alignment);
        uint8_t *seqi = bam_get_seq(alignment);
        uint8_t *qual = bam_get_qual(alignment);
        size_t n_cigar = alignment->core.n_cigar;

        Read read = {
            .mq = normalize_mq(alignment->core.qual),
            .read_start = pos,
            .cigartuples = cigartuples,
            .seqi = seqi,
            .qual = qual,
            .strand = normalize_strand(is_fwd_strand),
            .n_cigar = n_cigar,
            .l_qseq = l_qseq,
            .pos_info = NULL,
            .haplotype = HAP_UNPHASED};

        while (variant_current_pos < variant_num && variants[variant_current_pos]->position < pos)
            variant_current_pos++;
        variants_info.variant_current_pos = variant_current_pos;

        while (candidate_current_index < flanking_candidates_num && flanking_candidates[candidate_current_index] < pos)
            candidate_current_index++;

        read.read_end = get_read_end(cigartuples, n_cigar, read.read_start);

        // get the overlap candidates number and skip the alignment if no flanking candidate overlapped
        size_t overlap_candidates_num = get_overlap_candidate_num(pos, read.read_end, candidate_current_index, flanking_candidates_num, &flanking_candidates);
        read.overlap_candidates_num = overlap_candidates_num;
        if (read.overlap_candidates_num == 0)
        {
            read_no_overlap_num++;
            continue;
        }

        // haplotag the read following whatshap haplotagging logic
        if (need_haplotagging && alignment->core.qual >= min_haplotag_mq)
        {
            read.haplotype = haplotag_read(&variants_info, &read, ref_seq, ref_start);
        }

        pos_info = calloc(overlap_candidates_num, sizeof(Pos_info));
        for (size_t i = 0; i < overlap_candidates_num; i++)
        {
            pos_info[i].ins_bases = NULL;
            pos_info[i].ins_length = 0;
            pos_info[i].alt_base = 0;
            pos_info[i].del_length = 0;
            pos_info[i].bq = 0;
        }

        // index of current first overlapped flanking candidate
        size_t flanking_start = kh_int_counter_val(flanking_candidates_p, flanking_candidates[candidate_current_index]);
        read.flanking_start = flanking_start;

        // store all overlapped flanking candidates information and put all centered candidate information
        // into pos_alt_info struct
        size_t ref_pos = read.read_start;
        size_t query_pos = 0;
        for (size_t i = 0; i < n_cigar; i++)
        {
            size_t cigar_op = bam_cigar_op(cigartuples[i]);
            size_t length = bam_cigar_oplen(cigartuples[i]);
            if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF)
            {
                for (size_t p = ref_pos; p < ref_pos + length; p++)
                {
                    int flanking_index = kh_int_counter_val(flanking_candidates_p, p);
                    if (flanking_index != -1)
                    {
                        size_t offset = flanking_index - flanking_start;
                        pos_info[offset].alt_base = bam_seqi(seqi, query_pos);
                        pos_info[offset].bq = normalize_bq(qual[query_pos]);

                        int center_pos_index = kh_int_counter_val(candidates_p, p);
                        if (center_pos_index != -1)
                        {
                            char alt_base = seq_nt16_str[pos_info[offset].alt_base];
                            pos_alt_info[center_pos_index].acgt_count[acgt2num[alt_base - 'A']]++;
                            pos_alt_info[center_pos_index].depth++;
                        }
                    }
                    query_pos++;
                }
                ref_pos += length;
            }
            else if (cigar_op == BAM_CDEL)
            {

                int flanking_index = kh_int_counter_val(flanking_candidates_p, ref_pos - 1);
                if (flanking_index != -1)
                {
                    size_t offset = flanking_index - flanking_start;
                    pos_info[offset].del_length = length;
                    int center_pos_index = kh_int_counter_val(candidates_p, ref_pos - 1);
                    if (center_pos_index != -1)
                    {
                        kh_int_counter_add(pos_alt_info[center_pos_index].del_counter, length, 1);
                    }
                }
                for (size_t p = ref_pos; p < ref_pos + length; p++)
                {
                    int flanking_index = kh_int_counter_val(flanking_candidates_p, p);
                    if (flanking_index != -1)
                    {
                        size_t offset = flanking_index - flanking_start;
                        pos_info[offset].alt_base = -1;
                        int center_pos_index = kh_int_counter_val(candidates_p, p);
                        if (center_pos_index != -1)
                        {
                            pos_alt_info[center_pos_index].depth++;
                        }
                    }
                }
                ref_pos += length;
            }
            else if (cigar_op == BAM_CINS)
            {
                int flanking_index = kh_int_counter_val(flanking_candidates_p, ref_pos - 1);
                if (flanking_index != -1)
                {
                    size_t offset = flanking_index - flanking_start;
                    pos_info[offset].ins_bases = calloc(length + 1, sizeof(char));
                    for (size_t ins_idx = 0; ins_idx < length; ins_idx++)
                    {
                        pos_info[offset].ins_bases[ins_idx] = seq_nt16_str[bam_seqi(read.seqi, query_pos + ins_idx)];
                    }
                    pos_info[offset].ins_bases[length] = '\0';
                    pos_info[offset].ins_length = length;

                    int center_pos_index = kh_int_counter_val(candidates_p, ref_pos - 1);
                    if (center_pos_index != -1)
                    {
                        kh_counter_add(pos_alt_info[center_pos_index].ins_counter, pos_info[offset].ins_bases, 1);
                    }
                }
                query_pos += length;
            }
            else if (cigar_op == BAM_CREF_SKIP)
            {
                ref_pos += length;
            }
            else if (cigar_op == BAM_CSOFT_CLIP)
            {
                query_pos += length;
            }
        }

        //update the read array
        read.pos_info = pos_info;
        reads_num++;
        kv_push(Read, read_array, read);
    }

    // allocate memory of the input matrix of all candidates
    int8_t *matrix = calloc(candidate_num * matrix_depth * no_of_positions * channel_size, sizeof(int8_t));

    HAP read_hap_array[reads_num];
    int matrix_read_index_array[matrix_depth];
    Alt_info *alt_info = malloc(matrix_depth * sizeof(Alt_info));

    char **alt_info_p = calloc(candidate_num, sizeof(char*));
    fa_data data = calloc(1, sizeof(_fa_data));

    // loop each candiate and generate full-alignment input matrix
    for (size_t i = 0; i < candidate_num; i++)
    {
        size_t candidate = candidates[i];
        size_t start_pos = candidate - flanking_base_num;
        size_t end_pos = candidate + flanking_base_num + 1;
        size_t candidate_depth = 0;
        for (size_t j = 0; j < matrix_depth; j++)
        {
            alt_info[j].ins_bases = NULL;
            alt_info[j].alt_base = '\0';
            alt_info[j].del_length = 0;
            alt_info[j].has_alt_info = false;
        }

        for (size_t j = 0; j < matrix_depth; j++)
            matrix_read_index_array[j] = -1;

        size_t overlap_read_num = 0;
        for (size_t j = 0; j < reads_num; j++)
        {
            Read read = kv_A(read_array, j);
            if (read.read_start >= end_pos)
                break;
            if (read.read_end <= start_pos)
                continue;
            read_hap_array[overlap_read_num].read_index = j;
            read_hap_array[overlap_read_num++].haplotype = read.haplotype;
        }

        sort_read_name_by_haplotype(&read_hap_array, &matrix_read_index_array, matrix_depth, overlap_read_num);

        // loop each overlapped read of a candidate
        for (size_t d = 0; d < matrix_depth; d++)
        {
            int read_index = matrix_read_index_array[d];
            if (read_index == -1)
                continue;
            Read read = kv_A(read_array, read_index);
            int8_t hap_v = normalize_hap(read.haplotype);
            int8_t strand_v = read.strand;
            int8_t mq_v = read.mq;

            // loop all flanking position of a read
            for (size_t p = 0; p < no_of_positions; p++)
            {
                size_t cp = p + start_pos;
                size_t flanking_index = kh_int_counter_val(flanking_candidates_p, cp);
                int32_t offset = flanking_index - read.flanking_start;
                bool is_center_pos = p == flanking_base_num;

                if (read.pos_info[offset].alt_base < 0)
                {
                    if (is_center_pos)
                        candidate_depth++;
                    continue;
                }

                if (offset < 0 || offset >= read.overlap_candidates_num)
                    continue;

                int8_t alt_v = 0;
                char ref_base = upper_base(ref_seq[cp - ref_start]);
                int8_t ref_v = num2countbase_fa[ref_base - 'A'];
                int8_t bq_v = read.pos_info[offset].bq;

                if (is_center_pos)
                    candidate_depth++;
                size_t alt_int = read.pos_info[offset].alt_base;
                char alt_base = seq_nt16_str[read.pos_info[offset].alt_base];
                if (read.pos_info[offset].ins_length > 0)
                {
                    size_t ins_length = read.pos_info[offset].ins_length;
                    char *ins_bases = read.pos_info[offset].ins_bases;
                    int8_t ins_v = 0;
                    size_t max_ins_length = ins_length < no_of_positions - p ? ins_length : no_of_positions - p;

                    for (size_t ins_idx = 0; ins_idx < ins_length; ins_idx++)
                    {
                        char ins_alt_base = ins_bases[ins_idx];
                        if (ins_idx < max_ins_length && p < no_of_positions - 1)
                        {
                            ins_v = num2countbase_fa[ins_alt_base - 'A'];
                            matrix[i * offset_can + d * offset_row + (ins_idx + p) * offset_col + 6] = ins_v;
                        }
                    }
                    if (is_center_pos)
                    {
                        alt_info[d].alt_base = alt_base;
                        alt_info[d].ins_bases = ins_bases;
                        alt_info[d].has_alt_info = true;
                    }
                    alt_v = num2countbase_fa['I' - 'A'];
                }
                else if (read.pos_info[offset].del_length > 0)
                {
                    if (is_center_pos)
                    {
                        alt_info[d].del_length = read.pos_info[offset].del_length;
                        alt_info[d].has_alt_info = true;
                    }
                    alt_v = num2countbase_fa['D' - 'A'];
                }
                else if (ref_base - alt_base != 0)
                {
                    if (is_center_pos)
                    {
                        alt_info[d].alt_base = alt_base;
                        alt_info[d].has_alt_info = true;
                    }
                    alt_v = num2countbase_fa[alt_base - 'A'];
                }

                // update the matrix
                matrix[i * offset_can + d * offset_row + p * offset_col + 0] = ref_v;
                matrix[i * offset_can + d * offset_row + p * offset_col + 1] = alt_v;
                matrix[i * offset_can + d * offset_row + p * offset_col + 2] = strand_v;
                matrix[i * offset_can + d * offset_row + p * offset_col + 3] = mq_v;
                matrix[i * offset_can + d * offset_row + p * offset_col + 4] = bq_v;
                matrix[i * offset_can + d * offset_row + p * offset_col + 7] = hap_v;
            }
        }

        // finish the candidate proportion channel;
        candidate_depth = pos_alt_info[i].depth;
        for (size_t j = 0; j < matrix_depth; j++)
        {
            int8_t af_v = 0;
            if (alt_info[j].has_alt_info == false)
                continue;
            if (alt_info[j].ins_bases != NULL)
            {
                size_t count = kh_counter_val(pos_alt_info[i].ins_counter, alt_info[j].ins_bases);
                if (count > 0)
                    af_v = normalize_af(count / (float)candidate_depth);
            }
            else if (alt_info[j].del_length > 0)
            {
                size_t count = kh_int_counter_val(pos_alt_info[i].del_counter, alt_info[j].del_length);
                if (count > 0)
                    af_v = normalize_af(count / (float)candidate_depth);
            }
            else if (alt_info[j].alt_base != '\0')
            {
                size_t offset = alt_info[j].alt_base - 'A';
                af_v = normalize_af(pos_alt_info[i].acgt_count[acgt2num[offset]] / (float)candidate_depth);
            }

            if (af_v > 0)
            {
                for (size_t p = 0; p < no_of_positions; p++)
                {
                    if (matrix[i * offset_can + j * offset_row + p * offset_col + 0] != 0)
                        matrix[i * offset_can + j * offset_row + p * offset_col + 5] = af_v;
                }
            }
        }

        // store the alternative information into string
        size_t max_alt_length = 64;
        char *alt_info_str = calloc(max_alt_length, sizeof(char));
        char center_ref_base = upper_base(ref_seq[candidate - ref_start]);

        sprintf(alt_info_str, "%i-%i-%c-", candidate + 1, candidate_depth, center_ref_base);
        for (size_t j = 0; j < 4; j++)
        {
            if (j != acgt2num[center_ref_base - 'A'] && pos_alt_info[i].acgt_count[j] > 0)
                sprintf(alt_info_str + strlen(alt_info_str), "X%c %i ", ACGT[j], pos_alt_info[i].acgt_count[j]);
        }
        for (khiter_t k = kh_begin(pos_alt_info[i].ins_counter); k != kh_end(pos_alt_info[i].ins_counter); k++)
        {
            if (kh_exist(pos_alt_info[i].ins_counter, k))
            {
                char *key = kh_key(pos_alt_info[i].ins_counter, k);
                int val = kh_val(pos_alt_info[i].ins_counter, k);
                if (strlen(key) <= max_indel_length)
                {
                    if (strlen(alt_info_str) + strlen(key) + 32 >= max_alt_length)
                    {
                        while (strlen(alt_info_str) + strlen(key) + 32 >= max_alt_length)
                            max_alt_length = max_alt_length << 1;
                        alt_info_str = realloc(alt_info_str, max_alt_length*sizeof(char));
                    }
                    sprintf(alt_info_str + strlen(alt_info_str), "I%c%s %i ", center_ref_base, key, val);
                }
            }
        }

        for (khiter_t k = kh_begin(pos_alt_info[i].del_counter); k != kh_end(pos_alt_info[i].del_counter); k++)
        {
            if (kh_exist(pos_alt_info[i].del_counter, k))
            {
                int key = kh_key(pos_alt_info[i].del_counter, k);
                int val = kh_val(pos_alt_info[i].del_counter, k);
                if (key <= max_indel_length)
                {
                    if (strlen(alt_info_str) + key + 32 >= max_alt_length)
                    {
                        while (strlen(alt_info_str) + key + 32 >= max_alt_length)
                            max_alt_length = max_alt_length << 1;
                        alt_info_str = realloc(alt_info_str, max_alt_length*sizeof(char));
                    }
                    sprintf(alt_info_str + strlen(alt_info_str), "D%.*s %i ", key, ref_seq + candidate - ref_start + 1, val);
                }
            }
        }

        alt_info_p[i] = alt_info_str;

    } // end of candidate loop


    data->matrix = matrix;
    data->all_alt_info = alt_info_p;
    data->candidates_num = candidate_num;

    // free all allocated memory
    for (size_t j = 0; j < reads_num; j++)
    {
        Read read = kv_A(read_array, j);
        for (size_t p = 0; p < read.overlap_candidates_num; p++)
        {
            if (read.pos_info[p].ins_bases != NULL)
                free(read.pos_info[p].ins_bases);
        }
        free(read.pos_info);
    }

    for (size_t j = 0; j < candidate_num; j++)
    {
        kh_counter_destroy(pos_alt_info[j].ins_counter);
        kh_int_counter_destroy(pos_alt_info[j].del_counter);
    }

    free(chr);
    free(pos_alt_info);
    free(alt_info);
    kh_counter_destroy(read_name_set);
    kh_int_counter_destroy(candidates_p);
    kh_int_counter_destroy(flanking_candidates_p);
    kv_destroy(read_array);
    bam_destroy1(alignment);
    hts_itr_destroy(iter);
    fai_destroy(fai);

    return data;
}
