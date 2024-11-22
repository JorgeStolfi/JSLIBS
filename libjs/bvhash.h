/* Mapping arbitrary byte strings into 64-bit hashes */
/* Last edited on 2024-11-15 19:11:39 by stolfi */

#ifndef bvhash_H
#define bvhash_H

#include <stdint.h>
#include <stdlib.h>

uint64_t bvhash_bytes(void *x, size_t sz);
  /* Returns an integer hash of the {sz} bytes starting at {*x}.
    Suitable for hash tables; not cryptographically secure. */
  
#endif
