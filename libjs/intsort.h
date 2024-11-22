/* intsort - various algorithms for sorting integer vectors. */
/* Last edited on 2024-11-16 12:11:31 by stolfi */

#ifndef intsort_H
#define intsort_H

#include <stdint.h>

/* 
  These routines sort an array of integers {h[0..n-1]} according
  to a linear order with equalities specified by the parameters
  {cmp} and {sgn}.  
    
  More precisely, the integers {h[0..n-1]} are rearranged so that
  {sgn*cmp(h[i],h[j]) <= 0} for every valid index pair {i,j} with {i
  <= j}. The values returned by {cmp} should be consistent with an
  order relation (anti-symmetric, reflexive, and transitive). The
  {sgn} parameter is a client convenience feature, and is usually
  {+1} or {-1}. 
  
  NUMBER OF INVERSIONS
  
  In the comments, {K} stands for the number of inversions in the
  original list, namely index pairs {i,j} such that {i<j} but
  {sgn*cmp(h[i],h[j]) > 0}. Thus {K = 0} for a list that is
  already ordered, {K = n(n-1)/2} for a list in the reverse order, and
  {K = n(n-1)/4} (expected) for a randomly ordered list. */
  
void isrt_mergesort(int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
  /* Sorts {h[0..n-1]} by in-place merge-sort. Takes {O(n log(n))} 
    time in the worst case. */

void isrt_heapsort(int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
  /* Sorts {h[0..n-1]} by heap-sort.  Takes {O(n log(n))} time
    in the worst case. */

void isrt_binssort(int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
  /* Sorts {h[0..n-1]} by insertion sort with binary search. Makes
    {O(n*log(n))} comparisons but {O(K)} moves, where {K} is the number
    of inversions in the input order; which can be {O(n^2)} in the
    worst case. */

void isrt_inssort(int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
  /* Sorts {h[0..n-1]} by classic insertion sort. Makes
    {O(n+K)} comparisons and {O(K)} moves, where {K} is the number
    of inversions in the input order; which can be {O(n^2)} in the
    worst case. */

#endif
