#ifndef _MEDAKA_KHCOUNTER_H
#define _MEDAKA_KHCOUNTER_H

#include "khash.h"

typedef struct kh_counter_stats_t {
    size_t sum;
    size_t max;
} kh_counter_stats_t;

KHASH_MAP_INIT_STR(KH_COUNTER, int)
KHASH_MAP_INIT_INT(KH_INT_COUNTER, int)

// create a counter
static inline khash_t(KH_COUNTER) *kh_counter_init() {
    khash_t(KH_COUNTER) *h = kh_init(KH_COUNTER);
    return h;
}

static inline khash_t(KH_INT_COUNTER) *kh_int_counter_init() {
    khash_t(KH_INT_COUNTER) *h = kh_init(KH_INT_COUNTER);
    return h;
}

// Get a value from a counter 
int kh_counter_val(khash_t(KH_COUNTER) *hash, char *key);

// Clean up a counter
void kh_counter_destroy(khash_t(KH_COUNTER) *hash);

// Increment a counter by one
size_t kh_counter_increment(khash_t(KH_COUNTER) *hash, char *key);

size_t kh_counter_sub(khash_t(KH_COUNTER) *hash, char *key, int val);

// Increment a counter by a given amount
size_t kh_counter_add(khash_t(KH_COUNTER) *hash, char *key, int val);

// Retrieve statistics on counter
kh_counter_stats_t kh_counter_stats(khash_t(KH_COUNTER) *hash);

// Print contents of a counter
void kh_counter_print(khash_t(KH_COUNTER) *hash);

// similar to the kh_counter, except that the key is integer
int kh_int_counter_val(khash_t(KH_INT_COUNTER) *hash, int key);

size_t kh_int_counter_add(khash_t(KH_INT_COUNTER) *hash, int key, int val);

void kh_int_counter_destroy(khash_t(KH_INT_COUNTER) *hash);


#endif
