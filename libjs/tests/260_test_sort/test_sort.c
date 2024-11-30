#define PROG_NAME "test_sort"
#define PROG_DESC "test of sorting routines"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-22 20:43:47 by stolfi */
/* Created on 2004-11-02 (or earlier) by J. Stolfi, UNICAMP */

#define test_sort_COPYRIGHT \
  "Copyright © 2004  by the State University of Campinas (UNICAMP)"

#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include <affirm.h>
#include <vec.h>
#include <intmerge.h>
#include <intmerge_extra.h>
#include <intheap.h>
#include <intsort.h>
#include <intsort_extra.h>

/* Which tests to run: */
#define TEST_ihp_push_pop             TRUE

#define TEST_imrg_merge               TRUE
#define TEST_imrg_merge_pivot         TRUE
#define TEST_imrg_merge_symsplit      TRUE

#define TEST_isrt_inssort             TRUE
#define TEST_isrt_binssort            TRUE
#define TEST_isrt_heapsort            TRUE
#define TEST_isrt_heapsort_classic    TRUE
#define TEST_isrt_heapsort_vacsink    TRUE
#define TEST_isrt_mergesort           TRUE
#define TEST_isrt_mergesort_pivot     TRUE
#define TEST_isrt_mergesort_symsplit  TRUE
#define TEST_isrt_quicksort_median3   TRUE
#define TEST_isrt_quicksort_middle    TRUE

/* Max number of elements in test array: */
#define NMAX (1<<24)

/*
  TYPES AND PROTOTYPES

  In the following, {cmp} and {sgn} define an order for integers:
  an integer vector {h} is ordered iff {sgn*cmp(h[i-1],h[i]) <= 0} for all {i}.

  For the tests to be effective, the result of {cmp(a,b)} should bear
  no relation to {a < b}.  */

typedef int32_t int32_t_cmp_t(int32_t a, int32_t b); 
  /* A signed comparison predicate for integers (indices, etc.). */

typedef void tsr_sort_t(int32_t *h, uint32_t n, int32_t_cmp_t cmp, int32_t sgn);
  /* A procedure that sorts {h[0..n-1]} */

typedef void tsr_merge_t(int32_t *a, int32_t *b, int32_t *c, int32_t_cmp_t cmp, int32_t sgn);
  /* A procedure that, given two consecutive sorted blocks of elements
    in {h}, from {*a} through {*(b-1)} and from {b} through {*(cm-1)},
    interleaves them in order. */

void tsr_test_sort(int32_t *h, uint32_t n, tsr_sort_t sort, char *name, int32_t_cmp_t cmp);
/* Tests the sorting routine {sort} on the array {h[0..n-1]}, with comparison
  procedure {cmp} and both directions ({sgn = +1} and {sgn = -1}).
  The {name} should be the name of the sorting routine.

  Before the test, the array {h} is initialized with {h[i]=i}. */

void tsr_test_merge(int32_t *h, uint32_t n, tsr_merge_t merge, char *name, int32_t_cmp_t cmp);
/* Tests the merge routine {merge} on the array {h[0..n-1]}, with comparison
  procedure {cmp} and both directions ({sgn = +1} and {sgn = -1}).
  The {name} should be the name of the sorting routine.

  Before applying {merge}, the array {h} is initialized with {h[i]=i},
  divided into two unequal blocks, and each block is sorted with
  binsertion sort. */

void tsr_check_order(int32_t *h, uint32_t n, int32_t_cmp_t cmp, int32_t sgn, bool_t print);
/* Checks whether the order of {h[0..n-1]} is consistent with {cmp}
  and {sgn}. If {print = TRUE}, also prints the elements {h[i]} of {h}
  and their "hidden" values {c[h[i]]}. */

void tsr_parse_options (int32_t argc, char **argv, uint32_t *np);

/* THE TEST ORDER */

double_vec_t c;

int32_t tsr_cmp(int32_t i, int32_t j);
/* Compares {i} and {j} in a contrived order.
  (Actually, compares {c[i]} with {c[j]}.) */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    uint32_t n;
    tsr_parse_options(argc, argv, &n);
    int32_t h[n];

    /* Prepare table of actual values to be sorted: */
    c = double_vec_new(n);
    for (uint32_t i = 0; i < n; i++)
      { /* Insert about 10% equal values, if n >= 10: */
        c.e[i] = (i < (9*n + 9)/10 ? (1.0 + sin(22.2*i*i))/2 : c.e[n-i-1]);
      }

    fprintf(stderr, "Original order:\n");
    for (uint32_t i = 0;  i < n; i++) { h[i] = (int32_t)i; }
    for (uint32_t i = 0;  i < n; i++)
      { fprintf(stderr, "%5d %5d %6.3f\n", i, h[i], c.e[h[i]]); }

    /* Tests insertion sort routines first (used by other tests): */

    if (TEST_isrt_inssort)
      { tsr_test_sort(h, n, isrt_inssort, "isrt_inssort", tsr_cmp); }
    if (TEST_isrt_binssort)
      { tsr_test_sort(h, n, isrt_binssort, "isrt_binssort", tsr_cmp); }

    /* Tests quicksort algorithms: */

    if (TEST_isrt_quicksort_middle)
      { tsr_test_sort(h, n, isrt_quicksort_middle, "isrt_quicksort_middle", tsr_cmp); }
    if (TEST_isrt_quicksort_median3)
      { tsr_test_sort(h, n, isrt_quicksort_median3, "isrt_quicksort_median3", tsr_cmp); }

    /* Tests heap push/pop routines (used by some heapsorts): */
    if (TEST_ihp_push_pop)
      { fprintf(stderr, "Heapified:\n");
        uint32_t k = 0;
        while (k < n) { ihp_heap_insert(h, &k, h[k], tsr_cmp, +1); }
        for (uint32_t i = 0;  i < n; i++)
          { fprintf(stderr, "%5d %5d %6.3f\n", i, h[i], c.e[h[i]]); }

        fprintf(stderr, "De-heapified:\n");
        /* Pop them out in reverse order and store them at the end: */
        k = n;
        while (k > 0) { h[k-1] = ihp_heap_pop(h, &k, tsr_cmp, +1); }
        for (uint32_t i = 0;  i < n; i++)
          { fprintf(stderr, "%5d %5d %6.3f\n", i, h[i], c.e[h[i]]); }
      }

    /* Tests heapsort algorithms: */
    if (TEST_isrt_heapsort)
      { tsr_test_sort(h, n, isrt_heapsort, "isrt_heapsort", tsr_cmp); }
    if (TEST_isrt_heapsort_vacsink)
      { tsr_test_sort(h, n, isrt_heapsort_vacsink, "isrt_heapsort_vacsink", tsr_cmp); }
    if (TEST_isrt_heapsort_classic)
      { tsr_test_sort(h, n, isrt_heapsort_classic, "isrt_heapsort_classic", tsr_cmp); }

    /* Tests in-place merge routine (used by some mergesorts): */
    if (TEST_imrg_merge)
      { tsr_test_merge(h, n, imrg_merge, "imrg_merge", tsr_cmp); }
    if (TEST_imrg_merge_pivot)
      { tsr_test_merge(h, n, imrg_merge_pivot, "imrg_merge_pivot", tsr_cmp); }
    if (TEST_imrg_merge_symsplit)
      { tsr_test_merge(h, n, imrg_merge_symsplit, "imrg_merge_symsplit", tsr_cmp); }


    /* Tests in-place mergesorts: */
    if (TEST_isrt_mergesort)
      { tsr_test_sort(h, n, isrt_mergesort, "isrt_mergesort", tsr_cmp); }
    if (TEST_isrt_mergesort_pivot)
      { tsr_test_sort(h, n, isrt_mergesort, "isrt_mergesort_pivot", tsr_cmp); }
    if (TEST_isrt_mergesort_symsplit)
      { tsr_test_sort(h, n, isrt_mergesort_symsplit, "isrt_mergesort_symsplit", tsr_cmp); }

    return(0);
  }

int32_t tsr_cmp(int32_t i, int32_t j)
  { affirm((i >= 0) && (i < c.ne), "bad i");
    affirm((j >= 0) && (j < c.ne), "bad j");
    double ci = c.e[i], cj = c.e[j];
    if (ci < cj)
      { return -1; }
    else if (ci > cj)
      { return +1; }
    else
      { return 0; }
  }

/* TEST PROCEDURES */

void tsr_test_sort(int32_t *h, uint32_t n, tsr_sort_t sort, char *name, int32_t_cmp_t cmp)
  {
    fprintf(stderr, "============================================================\n");
    fprintf(stderr, "Sorting with %s:\n", name);
    
    bool_t print = (n <= 100);

    int32_t sgn;
    for (sgn = +1; sgn >= -1; sgn -= 2)
      { 
        fprintf(stderr, "Order = %d (%s):\n", sgn, (sgn > 0 ? "increasing" : "decreasing"));

        /* Restore {h[0..n-1]} to the original increasing order (random by {cmp}): */
        for (uint32_t i = 0;  i < n; i++) { h[i] = (int32_t)i; }

        /* Sort them by {sgn*cmp}: */
        for (uint32_t i = 0;  i < n; i++) { h[i] = (int32_t)i; }
        sort(h, n, cmp, sgn);
        
        /* Check and print: */
        tsr_check_order(h, n, cmp, sgn, print);
        
        fprintf(stderr, "\n");
     }

  }

void tsr_test_merge(int32_t *h, uint32_t n, tsr_merge_t merge, char *name, int32_t_cmp_t cmp)
  { 
    fprintf(stderr, "============================================================\n");
    fprintf(stderr, "Merging with %s:\n", name);
    
    uint32_t m = 3*n/5; /* Size of first block. */
    bool_t print = (n <= 100);

    int32_t sgn;
    for (sgn = +1; sgn >= -1; sgn -= 2)
      { 
        fprintf(stderr, "Order = %d (%s):\n", sgn, (sgn > 0 ? "increasing" : "decreasing"));
        
        /* Restore {h[0..n-1]} to the original increasing order (random by {cmp}): */
        for (uint32_t i = 0;  i < n; i++) { h[i] = (int32_t)i; }
        
        /* Sort the two blocks: */
        fprintf(stderr, "Sorting blocks [0..%d], [%d..%d]:\n", m-1, m, n-1);
        isrt_binssort(h, m, cmp, sgn);
        tsr_check_order(h, m, cmp, sgn, print);
        
        isrt_binssort(h+m, n-m, cmp, sgn);
        tsr_check_order(h+m, n-m, cmp, sgn, print);

        /* Merge the two blocks: */
        fprintf(stderr, "Merging blocks with %s:\n", name);
        merge(h, h+m, h+n, cmp, sgn);
        
        /* Check and print: */
        tsr_check_order(h, n, cmp, sgn, print);
        fprintf(stderr, "\n");
      }
  }

void tsr_check_order(int32_t *h, uint32_t n, int32_t_cmp_t cmp, int32_t sgn, bool_t print)
  {
    bool_t printed;
            
    auto void prt(uint32_t i);
    void prt(uint32_t i)
      { fprintf(stderr, "%5d %5d %6.3f", i, h[i], c.e[h[i]]);
        printed = TRUE;
      } 
      
    bool_t unstable = FALSE;
    bool_t buggy = FALSE;
    
    for (uint32_t i = 0; i < n; i++)
      { printed = FALSE;
        if (print) { prt(i); }
        if (i > 0)
          { int32_t b = sgn*cmp(h[i-1],h[i]);
            if (b > 0)
              { if (! printed) { prt(i); }
                fprintf(stderr, " ** out of order");
                buggy = TRUE;
              }
            else if ((b == 0) && (h[i] < h[i-1]))
              { if (! printed) { prt(i); }
                fprintf(stderr, " !! equal swapped");
                unstable = TRUE;
              }
          }
        if (printed) { fprintf(stderr, "\n"); }
      }
    affirm(! buggy, "** buggy sort");
    if (unstable)
      { fprintf(stderr, "!! unstable sort\n"); }
  }

void tsr_parse_options (int32_t argc, char **argv, uint32_t *np)
  { int64_t nn;
    char *p = NULL;
    
    if (argc != 2)
      { fprintf(stderr, "usage: %s NELEMS\n", PROG_NAME); exit(1); }
    nn = strtol(argv[1], &p, 10);
    if (((*p) != '\000') || (nn < 0) || (nn > NMAX))
      { fprintf(stderr, "invalid n = %ld\n", nn); exit(1); }

    (*np) = (uint32_t)nn;
  }
