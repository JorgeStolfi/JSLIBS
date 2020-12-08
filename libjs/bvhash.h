/* Mapping arbitrary byte strings into 64-bit hashes */
/* Last edited on 2015-10-09 15:31:43 by stolfilocal */

#ifndef bvhash_H
#define bvhash_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>

uint64_t bvhash_bytes(void *x, size_t sz);
  /* Returns an integer hash of the {sz} bytes starting at {*x}.
    Suitable for hash tables; not cryptographically secure. */
  
#endif
