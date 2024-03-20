
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "clair3_full_alignment.h"

int main(int argc, char *argv[]) {
    clock_t begin = clock();

    char *region = argv[1];
    char *bam = argv[2];
    char *fasta = argv[3];

    size_t variant_num = 10;
    Variant **variants = malloc(variant_num * sizeof(Variant));
    size_t candidate_num = 10;
    size_t *candidates = malloc(candidate_num * sizeof(size_t));
    bool need_haplotagging = true;
    size_t min_mq = 10;
    size_t minb_bq = 10;
    size_t matrix_depth = 100;
    size_t max_indel_length = 10;


    fa_data data = calculate_clair3_full_alignment(
        region, bam, fasta,
        variants, variant_num,
        candidates, candidate_num,
        need_haplotagging,
        min_mq, min_bq, matrix_depth, max_indel_length);
    destroy_fa_data(data);

    clock_t end = clock();
    fprintf(stderr, "Total CPU time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
    exit(EXIT_SUCCESS);
}
