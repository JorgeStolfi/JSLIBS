/* intmerge - in-place merging of integer lists. */
/* Last edited on 2024-11-15 19:13:24 by stolfi */

#ifndef intmerge_H
#define intmerge_H

#include <stdint.h>
  
void imrg_merge(int32_t *a, int32_t *b, int32_t *c, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
/* 
  Merges two consecutive blocks of an array of integers, from {*a}
  through {*(b-1)} and from {*b} through {*(c-1)}, leaving the result
  in the same place.
  
  The elements in each block and in the result are assumed to be
  ordered so that {sgn*cmp(*p,*q) <= 0} for every element addresses
  {p,q} in the same block with {p <= q}. The values returned by {cmp}
  should be consistent with an order relation (anti-symmetric,
  reflexive, and transitive). The {sgn} parameter is a client
  convenience feature, and is usually {+1} or {-1}.  */

#endif
