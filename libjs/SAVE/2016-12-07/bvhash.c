/* See {bvhash.h}. */
/* Last edited on 2015-10-09 15:31:35 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>

#include <bvhash.h>

typedef unsigned char ubyte;
    
uint64_t bvhash_bytes(void *x, size_t sz)
  { 
    /* FNV-1a algorithm: */
    ubyte *p = (ubyte *)x;
    uint64_t hash = 14695981039346656037LU;
    size_t i;
    for(i = 0; i < sz; i++)
      { hash ^= (uint64_t)(*p);
        hash *= 1099511628211LU;
        p++;
      }
    return hash;
  }

