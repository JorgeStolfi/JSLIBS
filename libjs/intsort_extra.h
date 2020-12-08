/* intsort_extra - medley of experimetnal sorting algorithms. */
/* Last edited on 2004-11-02 14:13:33 by stolfi */

#ifndef intsort_extra_H
#define intsort_extra_H

/* 
  CLASSIC HEAP-SORT */

void isrt_heapsort_classic(int *h, int n, int cmp(int x, int y), int sgn);
/* Classic textbook version of heapsort: insert the new element at root 
  and bubble it down. */

/*
  VACANCY-SINKING HEAPSORT */

void isrt_heapsort_vacsink(int *h, int n, int cmp(int x, int y), int sgn);
/* Improved heap-sort: sink root vacancy to a leaf,
  swap with new element, and bubble it up. */  

/* 
  QUICKSORT, PIVOT AT MIDDLE */

/* Threshold and sorter for oldquick-sort->binsertion-sort switch: */
#define isrt_quicksort_middle_SMALL 5
#define isrt_quicksort_middle_SMALLSORT isrt_inssort

void isrt_quicksort_middle(int *h, int n, int cmp(int x, int y), int sgn);
/* Quicksort pivoted on the middle element. For small lists, switches to
  {isrt_quicksort_middle_SMALLSORT}. */  

/*
  QUICKSORT, MEDIAN-OF-THREE PIVOT */

/* Threshold for newquick-sort->binsertion-sort switch: */
#define isrt_quicksort_median3_SMALL 13
#define isrt_quicksort_median3_SMALLSORT isrt_inssort

void isrt_quicksort_median3(int *h, int n, int cmp(int x, int y), int sgn);
/* Quicksort with median-of-three pivot. For small lists, switches to
  {isrt_quicksort_middle_SMALLSORT}. */


/* 
  MERGE-SORT, OPTIMIZED */

/* Threshold for merge-sort->binsertion-sort switch: */
#define isrt_mergesort_pivot_SMALL 1
#define isrt_mergesort_pivot_SMALLSORT isrt_inssort

void isrt_mergesort_pivot(int *h, int n, int cmp(int x, int y), int sgn); 
/* Experimental version of merge-sort. */


/* 
  MERGE-SORT, ALTERNATIVE */

/* Threshold for merge-sort->binsertion-sort switch: */
#define isrt_mergesort_symsplit_SMALL 1
#define isrt_mergesort_symsplit_SMALLSORT isrt_inssort

void isrt_mergesort_symsplit(int *h, int n, int cmp(int x, int y), int sgn);
/* Another variant of merge-sort. */

#endif
