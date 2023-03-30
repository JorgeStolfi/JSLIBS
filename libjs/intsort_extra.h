/* intsort_extra - medley of experimetnal sorting algorithms. */
/* Last edited on 2023-03-18 11:24:08 by stolfi */

#ifndef intsort_extra_H
#define intsort_extra_H

#define _GNU_SOURCE
#include <stdint.h>

/* 
  CLASSIC HEAP-SORT */

void isrt_heapsort_classic(int32_t *h, int32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
/* Classic textbook version of heapsort: insert the new element at root 
  and bubble it down. */

/*
  VACANCY-SINKING HEAPSORT */

void isrt_heapsort_vacsink(int32_t *h, int32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
/* Improved heap-sort: sink root vacancy to a leaf,
  swap with new element, and bubble it up. */  

/* 
  QUICKSORT, PIVOT AT MIDDLE */

/* Threshold and sorter for oldquick-sort->binsertion-sort switch: */
#define isrt_quicksort_middle_SMALL 5
#define isrt_quicksort_middle_SMALLSORT isrt_inssort

void isrt_quicksort_middle(int32_t *h, int32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
/* Quicksort pivoted on the middle element. For small lists, switches to
  {isrt_quicksort_middle_SMALLSORT}. */  

/*
  QUICKSORT, MEDIAN-OF-THREE PIVOT */

/* Threshold for newquick-sort->binsertion-sort switch: */
#define isrt_quicksort_median3_SMALL 13
#define isrt_quicksort_median3_SMALLSORT isrt_inssort

void isrt_quicksort_median3(int32_t *h, int32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
/* Quicksort with median-of-three pivot. For small lists, switches to
  {isrt_quicksort_middle_SMALLSORT}. */


/* 
  MERGE-SORT, OPTIMIZED */

/* Threshold for merge-sort->binsertion-sort switch: */
#define isrt_mergesort_pivot_SMALL 1
#define isrt_mergesort_pivot_SMALLSORT isrt_inssort

void isrt_mergesort_pivot(int32_t *h, int32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn); 
/* Experimental version of merge-sort. */


/* 
  MERGE-SORT, ALTERNATIVE */

/* Threshold for merge-sort->binsertion-sort switch: */
#define isrt_mergesort_symsplit_SMALL 1
#define isrt_mergesort_symsplit_SMALLSORT isrt_inssort

void isrt_mergesort_symsplit(int32_t *h, int32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn);
/* Another variant of merge-sort. */

#endif
