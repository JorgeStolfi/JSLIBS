#ifndef wt_median_H
#define wt_median_H

/* Weight tables for filtering digital signals */
/* Last edited on 2023-11-21 18:11:46 by stolfi */

#define wt_median_H_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <argparser.h>

/* WEIGHT TABLES WITHOUT ALLOCATION

  The procedures in this section create filter weight tables of 
  arbitrary length {n>0}, even or odd, stored into client-given arrays.  */

#define wt_median_WSUM_MAX (INT32_MAX/2 - 1)
  /* Maximum allowed weight sum. */

double wt_median_sorted
  ( int32_t ns,     /* Count of samples. */
    double xs[],    /* The samples are {xs[0..ns-1]}. */
    int32_t ws[],   /* The respective weights are {ws[0..ns-1]}. */
    bool_t interp   /* Should interpolate the median? */
  );
  /* Computes the weighted median of {ns} samples {xs[0..ns-1]} with
    weights {ws[0..ns-1]}. The numbers must be all distinct and finite,
    and must be sorted in strictly increasing order. 
    
    The weights {ws[0..ns-1]} must be striclty positive and their sum
    must be less than {INT32_MAX/2}.
    
    If {ns} is zero, the result is {NAN}. Otherwise, the procedure tries
    to return a value {xm} that minimizes the absolute value of {F(xm)},
    the sum of weights {ws[i]} of samples with {xs[i]<xm} minus the
    weights {ws[j]} of samples with {xs[j]>xm}.
    
    More precisely, if there is a sample {xs[i]} in the set {X} with
    {F(xs[i])==0}, the procedure returns that sample. This is
    the case, in particular, if {ns} is 1. The
    {interp} flag is ignored in this case.
    
    If there is no {i} such that {F(xs[i])==0}, then there will
    be two consecutive indices {ia,ib} in {0..ns-1} such that
    {F(xs[ia])<0} and {F(xs[ib])>0}. In that case, if {interp} is false,
    the procedure returns either {xs[ia]} or {xs[ib]}, depending on
    which of the two {F} values is closer to zero; breaking ties in
    favor of the even index.  If {interp} is true, the
    procedure will pretend that {F} is continuous piecewise linear
    (affine), and estimates its zero as a value intermediate 
    between {xs[ia]} and {xs[ib]} by linear (affine) interpolation. */

double wt_median_unsorted
  ( int32_t n,     /* Count of samples. */
    double x[],     /* The samples are {x[0..n-1]}. */
    int32_t w[],    /* The respective weights are {w[0..n-1]}. */
    bool_t interp,  /* Should interpolate the median? */
    int32_t *ns_P,  /* (OUT) Number of distinct sample values with nonzero weight. */
    double xs[],    /* (WORK) Condensed sample table with {n} entries, or {NULL}. */
    int32_t ws[],   /* (WORK) Condensed weight table with {n} entries, or {NULL}. */
    int32_t kx[]    /* (WORK) Index table with {n} entries, or {NULL}. */
  );
  /* Computes the weighted median of {ns} samples {x[0..n-1]}
    with weights {w[0..n-1]}.  Similar to {wt_median_sorted} but accepts 
    samples in any order and/or repeated, and zero weights.
    
    The weights {w[0..n-1]} must be non-negative, and their sum must be
    less than {INT32_MAX/2}.
    
    If all the weights are zero (in particular, if {n} is zero), the
    result is {NAN}. Otherwise, as in {wt_median_sorted}, the procedure
    tries to return a value {xm} that minimizes the absolute value of
    {F(xm)}, the sum of weights {w[i]} of samples with {x[i]<xm} minus
    the weights {w[j]} of samples with {x[j]>xm}.
    
    More precisely, if there is a sample {x[i]} in the set {X} with
    non-zero weight and {F(x[i])==0}, the procedure returns that sample.
    This is the case, in particular, if {n} is 1. The {interp} flag is
    ignored in this case.
    
    If there is no such {i}, then there must be two indices {ia,ib} in
    {0..n-1} such that {w[ia]} and {w[ib]} are positive, {F(x[ia])<0},
    {F(x[ib])>0}, and {F(x[ia])} and {F(x[ib])} are as close as possible
    to zero with those conditions. In that case, if {interp} is false,
    the procedure returns either {x[ia]} or {x[ib]}, depending on which
    of the two {F} values is closer to zero; breaking ties in favor of
    the even index. If {interp} is true, the procedure will pretend that
    {F} is continuous piecewise linear (affine), and estimates its zero
    as a value intermediate between {x[ia]} and {x[ib]} by linear
    (affine) interpolation.
    
    The procedure determines the number {ns} of /distinct/ sample values
    in {x[0..n-1]} that have positive weight, and returns that number in
    {*ns_P} if {ns_P} is not {NULL}.
    
    If the work arrays {ns,ws,kx} are not {NULL}, they must have at
    least {n} elements each. They are used by the procedure to compute
    the median. On exit, {xs[0..ns-1]} will be the distinct sample
    values with positive weight, {ws[i]} will be the total weight of
    sample value {xs[i]}, and {kx[0..n-1]} will be the indices of
    samples {x[0..n-1]} (including those with zero weight) in order of
    increasing value.  If any or all of these parameters is
    {NULL}, a work arrays is allocated internally and freed 
    before return.  */

int32_t wt_median_gather_samples
  ( int32_t nx,     /* Count of samples. */
    double x[],     /* The samples are {x[0..nx-1]}. */
    int32_t ix,     /* Lowest sample index in window. */
    int32_t nw,     /* Window width. */
    int32_t w[],    /* The window weights are {w[0..nw-1]}. */
    int32_t kx[],   /* (IN) The sorted window indices are {kx[0..nw-1]}. */
    double xs[],    /* (OUT) the rearranged and condensed samples are {xs[0..*ns-1]}. */
    int32_t ws[]    /* (OUT) the corresponding condensed weights are {ws[0..*ns-1]}. */
  );
  /* Given an array {x[0..nx-1]} of samples, and array {w[0..nw-1]} of weights, 
    and an array {kx[0..nw-1]}  of indices in {0..nx-1}, extracts the {nw} 
    samples {x[kx[0..nw-1]]}, together with the respective weights, 
    sorts them,  condenses repeated values, and stores the resulting {ns} samples and
    their condensed weights in {xs[0..ns-1]} and {ws[0..ns-1]}.
    
    The procedure assumes that the indices in {kx[0..nw-1]} are
    in the range {ix..ix+nw-1}, and that the {w[k]} is the weight of
    sample {x[ix+k]}.
    
    Samples and weights must be finite, and the samples must be
    non-negative. Samples with zero weight are omitted from the output.
    
    The procedure assumes that the array {kx[0..nw-1]} is sorted by
    sample value; that is, that the value of {x[kx[k]]} is
    non-decreasing as {k} ranges in {0..nw-1}.
    
    The condensed sample count {ns} is returned as the result. It will
    be zero if all weights are zero or {nw} is zero; otherwise it will
    be at least 1. */
    
void wt_median_index_set_sort(int32_t nx, double x[], int32_t nw, int32_t kx[], int32_t np);
void wt_median_index_set_quick_sort(int32_t nx, double x[], int32_t nw, int32_t kx[]);
void wt_median_index_set_insertion_sort(int32_t nx, double x[], int32_t nw, int32_t kx[]);
  /* These procedures sort the indices {kx[0..nw-1]}, which must be a subset of {0..nx-1}, so that 
    the elements {x[kx[i]]} are in non-decreasing order when {i} ranges from 0 to {nw-1}.
    
    The procedure {wt_median_index_set_quick_sort} uses the C {qsort}
    routine, while {wt_median_index_set_insertion_sort} uses insertion
    sort. The latter may be faster if the list {kx[0..nw-1]} is already
    almost sorted.

    The version {wt_median_index_set_sort} assumes that at most {np} elements of {kx[0..nw-1]}
    are out of order, and chooses between the quick and insertion sort based on that hint. */

#endif
