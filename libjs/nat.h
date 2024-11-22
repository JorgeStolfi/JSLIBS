#ifndef nat_H
#define nat_H

#include <stdint.h>

#include <vec.h>

/* Another name for "uint32_t" */
/* Last edited on 2024-11-15 19:14:51 by stolfi */

typedef uint32_t nat_t;

vec_typedef(nat_vec_t,nat_vec,nat_t);
  /* Vectors of integers ({nat_t}). */

#endif
