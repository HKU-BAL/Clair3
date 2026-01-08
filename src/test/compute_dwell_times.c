int16_t* compute_dwell_times_from_mv_tag(uint8_t *mv_tag_data, size_t l_qseq, bool is_reverse_strand)
{
    if (mv_tag_data == NULL) {
        return NULL;
    }

    // Parse mv tag - it's stored as array type 'C' (uint8_t)
    uint8_t type = mv_tag_data[0];
    if (type != 'C') {  // Check if it's uint8_t array
        return NULL;
    }

    // Get array length (stored as int32_t after type byte)
    int32_t mv_length = *(int32_t*)(mv_tag_data + 1);
    uint8_t *mv_array = mv_tag_data + 5;  // Skip type byte + 4-byte length

    if (mv_length <= 0) {
        return NULL;
    }

    // Allocate dwell times array
    int16_t *dwell_times = calloc(l_qseq, sizeof(int16_t));
    if (dwell_times == NULL) {
        return NULL;
    }

    // First mv value should be 1
    if (mv_array[0] != 1) {
        free(dwell_times);
        return NULL;
    }

    int16_t dwell_time = 1;
    size_t qpos = 0;

    // Process mv array to compute dwell times
    for (int32_t i = 1; i < mv_length && qpos < l_qseq; i++) {
        if (mv_array[i] == 1) {
            dwell_times[qpos] = dwell_time;
            dwell_time = 1;
            qpos++;
        } else {
            dwell_time++;
        }
    }

    // Handle last position
    if (qpos < l_qseq) {
        dwell_times[qpos] = dwell_time;
    }

    // Reverse if on reverse strand
    if (is_reverse_strand && qpos > 0) {
        for (size_t i = 0; i < qpos / 2; i++) {
            int16_t tmp = dwell_times[i];
            dwell_times[i] = dwell_times[qpos - 1 - i];
            dwell_times[qpos - 1 - i] = tmp;
        }
    }

    return dwell_times;
}