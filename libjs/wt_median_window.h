#ifndef wt_median_window_H
#define wt_median_window_H

/* Running median filter tools */
/* Last edited on 2024-11-22 03:24:11 by stolfi */

#define wt_median_window_H_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <argparser.h>

#define wt_median_window_WSUM_MAX (INT32_MAX/2 - 1)
  /* Maximum allowed weight sum. */
    
double wt_median_window
  ( uint32_t nx,    /* Total number of samples. */ 
    double x[],     /* Samples are {x[0..nx-1]}. */
    int32_t ix,     /* First sample in window. */
    uint32_t nw,    /* Count of samples in window. */
    uint64_t w[],   /* Weights of samples in window are {w[0..nw-1]}. */
    bool_t interp,  /* True interpolates median, false returns nearest sample. */
    uint32_t nk,    /* (IN) number of sample indices in {kx}. */
    uint32_t kx[],  /* (IN/OUT) prev sample indices; (OUT) {ix..jx} sorted by weight. */
    uint32_t *ns_P, /* (OUT) Num of distinct sample values with nonzero weight. */
    double xs[],    /* (OUT) condensed window samples are {xs[0..ns-1]}. */
    uint64_t ws[]   /* (OUT) respective condensed weights are {ws[0..ns-1]}. */
  );
  /* Finds the weighted median of samples {x[ix..jx]}, where {jx=ix+nw-1},
    with weights {w[0..nw-1]}. 
    
    Those samples and weights are sorted, then equal values are condensed
    as described under {wt_median_gather_samples}. The sorted and
    condensed samples and weights are stored into {xs[0..ns-1]} and
    {ws[0..ns-1]}, respectively, and the count {ns} is returned in
    {*ns_P}.  These vectors must have at least {nw} elements each.
    
    The {interp} parameter has the same meaning as in {wt_median_window}. Note
    that samples with zero weight are ignored for all purposes, so that
    the result is either a sample value with nonzero weight or is
    interpolated between the two closest values with nonzero weight that
    bracket the median.
    
    The parameters {nk} and {kx} are hints that can be used to speed up
    the computation. On entry, {nk} should be at most {nx} and
    {kx[0..nk-1]} should be a list of distinct CONSECUTIVE indices
    {ix0..jx0} in the range {0..nx-1}, in any order. The procedure
    updates that list with {wt_median_window_index_set_update}, so that
    it becomes a list {kx[0..nw-1]} that is a permutation of {ix..jx}.
    Note that the new list length {nw} may be greater than {nk}, so the
    vector {kx} must have at least {nw} elements. The procedure then
    sorts {kx[0..nw-1]} with {wt_median_index_set_sort}, so that
    {x[kx[k]]} is strictly non-decreasing as {k} ranges from 0 to
    {nw-1}. The vector {kx} is then used for
    {wt_median_window_gather_samples}. When computing multiple medians
    of multiple overlapping windows, the use of {kx} between successive
    windows may reduce the sorting time substantially. */

uint32_t wt_median_window_index_set_update
  ( uint32_t nx,     /* Total count of samples. */
    uint32_t nk,     /* (IN/OUT) Count of indices in {kx}. */
    uint32_t kx[],   /* (IN/OUT) Sorted indices are {kx[0..nk-1]}. */
    int32_t ix,      /* Index of first sample in new window. */
    uint32_t nw      /* Count of samples in new window. */
  );
  /* This procedure can be used to speed up the computation of medians
    (specifically, of {wt_median_gather_samples}) of many overlapping
    windows.  It updates the set of sample indices of a previous 
    window to those of a new window, while preserving the order of
    the indices that are in both windows.
  
    Specifically, on entry {nk} must be at most {nx}, and {kx[0..nk-1]}
    must be a list of {nk} distinct CONSECUTIVE indices {ix0..ij0} from
    the set {0..nx-1}, in any order. The parameter {ix} must be an index
    {0..nx-1}, and {nw} must be at most {nx}.
    
    On exit, {kx} will contain the {nw} consecutive indices {ix..jx} of
    the samples in the new window, where {jx=ix+nw-1}, in some order.
    This too must be a sub-range of {0..nx-1}. The new list size {nw}
    may be bigger or smaller than the input size {nk}, thus the vector
    {kx} must have at least {nw} elements.
    
    Indices in the input list {ix0..jx0} which are not in the new
    range {ix..jx} are deleted, and those which are in that range are
    copied into consecutive positions, starting at {kx[0]}, preserving
    their order. Indices in {ix..jx} which were NOT in {ix0..jx0},
    if any, are then stored into {kx} after all the preserved ones. The
    procedure returns the number of {nkept} of indices that were preserved,
    which is in {0..nk1}.
    
    Thus, if the input index list {kx[0..nk0-1]} was sorted by some
    criterion, on output the indices {kx[0..nkept-1]} will stil be
    sorted, while the indices {kx[nkept..nk1-1]} will need to be sorted
    and inserted among them. Preserving the ordering of the indices that
    are not removed is the whole point of this procedure.
    
    In particular, if {nk} is zero or the old window indices {ix0..jx0} are
    disjoint from the new ones {ix..jx}, the entries {kx[0..nw-1]} will be
    set {ix..jx}, in increasing order. */

#endif
