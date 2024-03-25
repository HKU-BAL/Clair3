#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "clair3_full_alignment.h"
#include "csv.h"

typedef struct Variants {
    size_t num_variants;
    Variant** variants;
} Variants;


typedef struct Candidates {
    size_t num_candidates;
    size_t* candidates;
} Candidates;


Variants read_variants(char *input) {
    // reads a file laid out as:
    // <num_variants>
    // pos,ref_base,alt_base,genotype,phase_set
    // ...
    char* row;
    char* col;
    CsvHandle handle = CsvOpen(input);
    if (handle == NULL) {
        fprintf(stderr, "Failed to read file\n");
        CsvClose(handle);
        return (Variants){0, NULL};
    }

    // first row is number of subsequent lines
    row = CsvReadNextRow(handle);
    size_t num_variants = atoi(row);
    Variant** variants = malloc(num_variants * sizeof(Variant*));

    int nrows = 0;
    Variant* var;
    while (row = CsvReadNextRow(handle)) {
        var = malloc(sizeof(Variant));
        var->position = atoi(CsvReadNextCol(row, handle));
        memcpy(&(var->ref_base), CsvReadNextCol(row, handle), 1);
        memcpy(&(var->alt_base), CsvReadNextCol(row, handle), 1);
        var->genotype = atoi(CsvReadNextCol(row, handle));
        var->phase_set = atoi(CsvReadNextCol(row, handle));
        variants[nrows] = var;
        nrows++;
        if (nrows == num_variants) break;
    }
    //fprintf(
    //    stderr, "%d, %c, %c, %d, %d\n",
    //    var->position, var->ref_base, var->alt_base, var->genotype, var->phase_set);
    CsvClose(handle);
    return (Variants){num_variants, variants};
}


void free_variants(Variants vars) {
    for (size_t i=0; i<vars.num_variants; ++i) {
        free(vars.variants[i]);
    }
    free(vars.variants);
}


Candidates read_candidates(char* input) {
    // reads a file laid out as:
    // <num_candidates>
    // cand1,cand2,cand3,cand4,...,candN
    CsvHandle handle = CsvOpen(input);
    if (handle == NULL) {
        fprintf(stderr, "Failed to read file\n");
        CsvClose(handle);
        return (Candidates){0, NULL};
    }

    // first row is number of subsequent lines
    char* row = CsvReadNextRow(handle);
    size_t num_candidates = atoi(row);
    size_t* candidates = malloc(num_candidates * sizeof(size_t));

    int ncols = 0;
    row = CsvReadNextRow(handle);
    char* col;
    for (size_t i=0; i<num_candidates; ++i) {
        candidates[i] = atoi(CsvReadNextCol(row, handle));
        //fprintf(stderr, "%zu\n", candidates[i]);
    }
    CsvClose(handle);
    return (Candidates){num_candidates, candidates};
}


void free_candidates(Candidates cand) {
    free(cand.candidates);
}


int main(int argc, char *argv[]) {
    // input data as dumped from Python
    //b'chr15:52446617-58814563'
    //b'chr15:52446617-58814563'
    //b'hg002_hac.pass.bam'
    //b'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
    //<cdata 'struct Variant *[]' owning 466496 bytes>
    //58312
    //<cdata 'size_t[10000]' owning 80000 bytes>
    //10000
    //True
    //5
    //0
    //89
    //50
    clock_t begin = clock();

    char *region = argv[1];             // "chr15:52446617-58814563";
    char *bam = argv[2];                // "hg002_hac.pass.bam";
    char *fasta = argv[3];              // "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna";
    char* variants_file = argv[4];      // fa_test_variants.txt
    char* candidates_file = argv[5];    // fa_test_candidates.txt

    Variants variants = read_variants(variants_file);
    if (variants.num_variants == 0) exit(1);
    Candidates candidates = read_candidates(candidates_file);
    if (candidates.num_candidates == 0) exit(1);

    // could take these as input?
    bool need_haplotagging = true;
    size_t min_mq = 5;
    size_t minb_bq = 0;
    size_t matrix_depth = 89;
    size_t max_indel_length = 50;

    fa_data data = calculate_clair3_full_alignment(
        region, bam, fasta,
        variants.variants, variants.num_variants,
        candidates.candidates, candidates.num_candidates,
        need_haplotagging,
        min_mq, min_bq, matrix_depth, max_indel_length);
    // TODO: do something with the data?
    destroy_fa_data(data);

    free_variants(variants);
    free_candidates(candidates);

    clock_t end = clock();
    fprintf(stderr, "Total CPU time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
    exit(EXIT_SUCCESS);
}
