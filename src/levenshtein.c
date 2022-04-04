// `levenshtein.c` - levenshtein
// MIT licensed.
// Copyright (c) 2015 Titus Wormer <tituswormer@gmail.com>

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "levenshtein.h"

// Returns a size_t, depicting the difference between `a` and `b`.
// See <https://en.wikipedia.org/wiki/Levenshtein_distance> for more information.
size_t
levenshtein_n(const char *a, const size_t length, const char *b, const size_t bLength) {
  // Shortcut optimizations / degenerate cases.
  if (a == b) {
    return 0;
  }

  if (length == 0) {
    return bLength;
  }

  if (bLength == 0) {
    return length;
  }

  size_t *cache = calloc(length, sizeof(size_t));
  size_t index = 0;
  size_t bIndex = 0;
  size_t distance;
  size_t bDistance;
  size_t result;
  char code;

  // initialize the vector.
  while (index < length) {
    cache[index] = index + 1;
    index++;
  }

  // Loop.
  while (bIndex < bLength) {
    code = b[bIndex];
    result = distance = bIndex++;
    index = SIZE_MAX;

    while (++index < length) {
      bDistance = code == a[index] ? distance : distance + 1;
      distance = cache[index];

      cache[index] = result = distance > result
        ? bDistance > result
          ? result + 1
          : bDistance
        : bDistance > distance
          ? distance + 1
          : bDistance;
    }
  }

  free(cache);

  return result;
}

size_t
levenshtein(const char *a, const char *b) {
  const size_t length = strlen(a);
  const size_t bLength = strlen(b);

  return levenshtein_n(a, length, b, bLength);
}
