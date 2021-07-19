/* Multidimensional sample arrays stored as k-d-trees. */
/* Last edited on 2021-07-16 09:48:40 by jstolfi */

#ifndef kdtom_test_H
#define kdtom_test_H

#define _GNU_SOURCE
#include <stdio.h>

#include <ppv_array.h>
#include <bool.h>

#include <kdtom.h>
#include <kdtom_array.h>
#include <kdtom_const.h>
#include <kdtom_split.h>

ppv_sample_t *kdtom_test_pick_max_samples(int32_t nms);
  /* Returns a vector of {ppv_sample_t} containing various values of
    {maxsmp} to use in tests, that require from 0 to {ppv_MAX_BPS} bits
    per sample in roughly linear progression. */

void kdtom_test_choose_array_size(ppv_dim_t d, ppv_size_t size[]);
  /* Fills {size[0..d-1]} with suitable values, not all equal. */

ppv_array_t *kdtom_test_array_make(ppv_dim_t d, ppv_size_t size[], ppv_sample_t maxsmp);
  /* Creates an array {A} with the given attributes {d,size[0..d-1],maxsmp},
    and fills it with a suitable test pattern. */

typedef ppv_sample_t kdtom_test_get_sample_proc_t(ppv_index_t ix[]);
  /* Type of a procedure that returns a sample value from a index. */

void kdtom_test_get_sample(kdtom_t *T, kdtom_test_get_sample_proc_t *getsmp);
  /* Tests {kdtom_get_sample(T,ix)} for points with varous indices {ix} just inside and 
    just outside the core domain {T.DK}.  Checks whether the index vectors outside {T.DK} yield
    {T.fill}, and those inside {T.DK} yield {getsmp(ix)}. */  

void kdtom_test_translate(kdtom_t *T);
  /* Tests {kdtom_translate} with some arbitrary displacement vector {dx[0..T.d-1]}.
    
    The translated version {S} of {T} should have {S.DK} equal to {T.DK}
    eaccept for the translation by {dx}. The sample values of {S} will
    be tested with {kdtom_test_get_sample}. The test requires {S.V[ix]}
    be the same as {T.V[ix-dx]} if {ix} is in {S.DK}, and equal to
    {T.fill} otherwise. */

void kdtom_test_clip_core(kdtom_t *T, ppv_axis_t ax, ppv_index_t ixlo_ax, ppv_size_t size_ax);
  /* Tests {kdtom_clip_core(T,ixlo,size)} with a box that 
    is slightly bigger than the core domain {T.DK}, except that
    along axis {ax} is spans indices {ixlo_ax .. ixlo_ax+size_ax-1}.
    
    Each clipped version {S} of {T} will be tested with
    {kdtom_test_get_sample}. The test requires {S.V[ix]} be the same as
    {T.V[ix]} if {ix} is in {S.DK}, and equal to {T.fill} otherwise. */

int32_t ktdom_test_remove_dup_indices(int32_t n, ppv_index_t ix[]);
  /* Sorts the array {ix[0..n-1]} in increasing order, removing duplicates. 
    Returns the number {nu} of unique indices, and those indices in {ix[0..nu-1]}. */

typedef void kdtom_test_range_proc_t(ppv_index_t ixlo, ppv_size_t size);
  /* Type of a procedure that processes a range of indices {ixlo..ixlo+size-1}. */

void kdtom_test_enum_ranges
  ( int32_t nlo, 
    ppv_index_t ixlo_hot[], 
    int32_t nhi, 
    ppv_index_t ixhi_hot[], 
    kdtom_test_range_proc_t process
  );
  /* Enumerates a bunch of index ranges {ixlo..ixhi} and calls
    {process(ixlo,size)} on them, where {size = ixhi-ixlo+1}.
    
    The low index {ixlo} will be taken from the list
    {ixlo_hot[0..nlo-1]} and the high index {ixhi} will be taken from
    the list {ixhi_hot[0..nhi-1]}, in all combinations. The {process}
    routine will be called at most once with an empty range. */

void kdtom_test_enum_ranges_single(ppv_index_t ixlo, ppv_size_t size, kdtom_test_range_proc_t process);
  /* Enumerates a bunch of index ranges {ixlo_r..ixhi_r} and calls 
    {process(ixlo_r,size_r)} on them, where {size_r = ixhi_r - ixlo_r + 1}.
    
    The low and high indices of those ranges will be just below {ixlo}, at {ixlo},
    somewhere between {ixlo} and {ixhi = ixlo+size-1}, at {ixhi}, and just above {ixhi}.
    Exactly one of the ranges will be empty. */

void kdtom_test_enum_ranges_split
  ( ppv_index_t ixlo,
    ppv_size_t size, 
    ppv_size_t size0,
    kdtom_test_range_proc_t process
  );
  /* Like {kdtom_test_enum_ranges}, but the interesting values that may be tried include
    also just below {ixhi1 = ixlo + size0} and at {ixhi1}. */

void kdtom_test_paint_bullseye(ppv_array_t *A, double ctr[], double R);
  /* Paints a suitable test pattern into {A}. Namely, two thin spherical
    shells where the original values of {A} are preserved,
    constant values outside, between, and inside the shells,
    and a different value at the center. The pattern will 
    be centered at the point {ctr}, and the outer radius of the
    outer shell will be {R}. */
    
void kdtom_test_check_tree(kdtom_t *T);
  /* Prints data about the tree {T}. */

void kdtom_test_plot(kdtom_t *T);
  /* Plots the array {T}.  Only if {T.d} is 2. */


#endif
