#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>
#define READ_LIMIT 10

int16_t* compute_dwell_times_from_mv_tag(uint8_t *mv_tag_data, size_t l_qseq, bool is_reverse_strand)
{
    if (mv_tag_data == NULL) {
        return NULL;
    }
    if (mv_tag_data[0] != 'B'){
        fprintf(stderr, "mv tag has unexpected layout (type=%c)\n", mv_tag_data[0]);
        return NULL;
    }
    int n = bam_auxB_len(mv_tag_data);
    if (n <= 0) {
        fprintf(stderr, "mv tag has unexpected length (length=%d)\n", n);
        return NULL;
    }
    printf("n = %d\n", n);
    for (int i = 0; i < n; i++) {
        printf(" %d", bam_auxB2i(mv_tag_data, i));
    }
    printf("\n");
    exit(0);
    return 0;

}

static int unpack_mv_tag(uint8_t *mv_tag_data, int32_t *mv_length, uint8_t **mv_array)
{
    if (mv_tag_data == NULL || mv_length == NULL || mv_array == NULL)
    {
        return -1;
    }

    uint8_t type = mv_tag_data[0];

    if (type == 'C')
    {
        *mv_length = *(int32_t *)(mv_tag_data + 1);
        *mv_array = mv_tag_data + 5;
        return 0;
    }

    if (type == 'B' && (mv_tag_data[1] == 'C' || mv_tag_data[1] == 'c'))
    {
        *mv_length = *(int32_t *)(mv_tag_data + 2);
        *mv_array = mv_tag_data + 7;
        return 0;
    }
    fprintf(stderr, "mv tag has unexpected layout (type=%c)\n", type);
    return -1;
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s <input.bam>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *bam_path = argv[1];
    samFile *bam_file = sam_open(bam_path, "r");
    if (bam_file == NULL)
    {
        fprintf(stderr, "Failed to open BAM file: %s\n", bam_path);
        return EXIT_FAILURE;
    }

    bam_hdr_t *header = sam_hdr_read(bam_file);
    if (header == NULL)
    {
        fprintf(stderr, "Failed to read BAM header: %s\n", bam_path);
        sam_close(bam_file);
        return EXIT_FAILURE;
    }

    bam1_t *record = bam_init1();
    if (record == NULL)
    {
        fprintf(stderr, "Failed to allocate BAM record\n");
        bam_hdr_destroy(header);
        sam_close(bam_file);
        return EXIT_FAILURE;
    }

    int processed = 0;
    while (processed < READ_LIMIT && sam_read1(bam_file, header, record) >= 0)
    {
        uint8_t *mv_tag_data = bam_aux_get(record, "mv");
        if (mv_tag_data == NULL)
        {
            continue;
        }

        int32_t mv_length = 0;
        uint8_t *mv_array = NULL;
        if (unpack_mv_tag(mv_tag_data, &mv_length, &mv_array) != 0 || mv_length <= 0)
        {
            fprintf(stderr, "[%d] %s: mv tag has unexpected layout (type=%c)\n",
                    processed + 1, bam_get_qname(record), mv_tag_data[0]);
            continue;
        }

        bool is_reverse = bam_is_rev(record);
        // printf("[%d] %s strand=%c mv_len=%d\n",
        //        processed + 1,
        //        bam_get_qname(record),
        //        is_reverse ? '-' : '+',
        //        mv_length);

        // printf("    mv:");
        // for (int32_t i = 0; i < mv_length; ++i)
        // {
        //     printf(" %u", mv_array[i]);
        // }
        // printf("\n");

        int16_t *dwell_times = compute_dwell_times_from_mv_tag(
            mv_tag_data, record->core.l_qseq, is_reverse);
        if (dwell_times == NULL)
        {
            printf("    compute_dwell_times_from_mv_tag -> NULL\n");
        }
        else
        {
            int preview = record->core.l_qseq < 10 ? record->core.l_qseq : 10;
            printf("    dwell (first %d bp):", preview);
            for (int i = 0; i < preview; ++i)
            {
                printf(" %d", dwell_times[i]);
            }
            if (record->core.l_qseq > preview)
            {
                printf(" ...");
            }
            printf("\n");
        }
        free(dwell_times);

        processed++;
    }

    if (processed == 0)
    {
        fprintf(stderr, "No reads with mv tag found within the first %d alignments\n", READ_LIMIT);
    }

    bam_destroy1(record);
    bam_hdr_destroy(header);
    sam_close(bam_file);

    return EXIT_SUCCESS;
}

