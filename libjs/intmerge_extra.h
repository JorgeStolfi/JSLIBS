/* intmerge_extra - alternative algorithms for merging of integer lists. */
/* Last edited on 2024-11-15 19:13:29 by stolfi */

#ifndef intmerge_extra_H
#define intmerge_extra_H
  
#include <stdint.h>
 
void imrg_merge_pivot(int32_t *a, int32_t *b, int32_t *c, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
/* 
  Alternative version of {imrg_merge}: chooses the middle element of the
  larger block as a pivot, and locates it in the smaller array by binary
  search.  Swaps blocks by decomposition into cycles.  */

void imrg_merge_symsplit(int32_t *a, int32_t *b, int32_t *c, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
/* 
  Alternative version of {imrg_merge}, uses mirror-symmetric binary
  search to find the best split points (corresponding to a median cut
  of the merged file).  */

#endif
