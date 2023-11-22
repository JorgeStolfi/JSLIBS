#ifndef wt_window_median_H
#define wt_window_median_H

/* Running median filter tools */
/* Last edited on 2023-11-21 09:47:03 by stolfi */

#define wt_window_median_H_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <argparser.h>

#define wt_window_median_WSUM_MAX (INT32_MAX/2 - 1)
  /* Maximum allowed weight sum. */
    
double wt_window_median
  ( int32_t nx,    /* Total number of samples. */ 
    double x[],    /* Samples are {x[0..nx-1]}. */
    int32_t ix,    /* First sample in window. */
    int32_t nw,    /* Window width. */
    int32_t w[],   /* Weights of samples in window are {w[0..nw-1]}. */
    bool_t interp, /* True interpolates median between sample values, false returns nearest sample. */
    int32_t nk,    /* (IN) number of sample indices in {kx}. */
    int32_t kx[],  /* (IN) {kx[0..nk-1]} some sample indices; (OUT) sorted window sample indices. */
    int32_t *ns_P, /* (OUT) Num of distinct sample values with nonzero weight. */
    double xs[],   /* (OUT) condensed window samples are {xs[0..ns-1]}. */
    int32_t ws[]   /* (OUT) respective condensed weights are {ws[0..ns-1]}. */
  );
  /* Finds the weighted median of samples {x[ix..jx]}, where {jx=ix+nw-1},
    with weights {w[0..nw-1]}. 
    
    Those samples and weights are sorted and equal values are condensed
    as described under {wt_median_gather_samples}. The sorted and
    condensed samples and weights are stored into {xs[0..ns-1]} and
    {ws[0..ns-1]}, respectively, and the count {ns} is returned in
    {*ns_P}.  These vectors must have at least {nw} elements each.
    
    The {interp} parameter has the same meaning as in {wt_window_median}. Note
    that samples with zero weight are ignored for all purposes, so that
    the result is either a sample value with nonzero weight or is
    interpolated between the two closest values with nonzero weight that
    bracket the median.
    
    The parameters {nk} and {kx} are hints that can be used to speed up
    the computation.  On entry, {kx[0..nk-1]} should be a list
    of distinct consecutive indices in the range {0..nx-1}. The
    procedure updates that list with {wt_window_median_index_set_update}, so
    that is becomes a permutation of {ix..jx}. It is OK if {nk} is
    zero or greater than {nw}, but the vector {kx} must have at least
    {nw} elements. The procedure then sorts {kx[0..nw-1]} with
    {wt_median_index_set_sort}, so that {x[kx[k]]} is strictly
    non-decreasing as {k} ranges from 0 to {nw-1}. The vector {kx} is
    then used for {wt_window_median_gather_samples}. When computing multiple
    medians of  multiple overlapping windows, the use of {kx} between
    successive windows may reduce the sorting time substantially. */

int32_t wt_window_median_index_set_update
  ( int32_t nx,     /* Count of samples. */
    int32_t ix,     /* Index of first sample in window. */
    int32_t nw,     /* Window width. */
    int32_t nk,     /* (IN) Count of indices in {kx}. */
    int32_t kx[]    /* (IN/OUT) Sorted indices are {kx[0..nk-1]}. */
  );
  /* This procedure can be used to speed up the computation of medians
    (specifically, of {wt_median_gather_samples}) of many overlapping
    windows.
  
    The procedeure assumes that, on entry, {kx[0..nk-1]} is an arbitrary
    list of {nk} distinct indices from the set {0..nx-1}. On exit,
    {kx[0..nw-1]} will contain the {nw} consecutive indices
    {ix..jx}, where {jx=ix+nw-1}, which must be a sub-range of {0..nx-1}.
    The input set size {nk} may be bigger or smaller than the window
    width {nw}, but the vector {kx} must have at least {nw} elements.
    
    Indices in {kx[0..nk-1]} which are not in the range
    {ix..jx} are deleted, and those which are in that range are
    copied into consecutive positions, preserving their order.  Indices
    in {ix..jx} which were not in {kx[0..nk-1]}, if any, are then
    appended at the end of {kx}.  The procedure returns the
    number of {np} of indices that were added, which is in {0..nw} 
    (but may be more than {nw-nk}).
    
    Thus, if the input index list {kx[0..nk-1]} was sorted by some
    criterion, on output the indices {kx[0..nw-np-1]} will stil be
    sorted, while the indices {kx[nw-np..nw-1]} will need to be sorted
    and inserted among them. */

#endif
