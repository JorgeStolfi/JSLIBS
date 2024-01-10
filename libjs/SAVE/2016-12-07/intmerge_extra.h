/* intmerge_extra - alternative algorithms for merging of integer lists. */
/* Last edited on 2004-11-02 17:48:19 by stolfi */

#ifndef intmerge_extra_H
#define intmerge_extra_H
  
void imrg_merge_pivot(int *a, int *b, int *c, int cmp(int x, int y), int sgn);
/* 
  Alternative version of {imrg_merge}: chooses the middle element of the
  larger block as a pivot, and locates it in the smaller array by binary
  search.  Swaps blocks by decomposition into cycles.  */

void imrg_merge_symsplit(int *a, int *b, int *c, int cmp(int x, int y), int sgn);
/* 
  Alternative version of {imrg_merge}, uses mirror-symmetric binary
  search to find the best split points (corresponding to a median cut
  of the merged file).  */

#endif
