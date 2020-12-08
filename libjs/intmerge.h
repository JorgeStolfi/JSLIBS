/* intmerge - in-place merging of integer lists. */
/* Last edited on 2004-10-31 16:53:02 by stolfi */

#ifndef intmerge_H
#define intmerge_H
  
void imrg_merge(int *a, int *b, int *c, int cmp(int x, int y), int sgn);
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
