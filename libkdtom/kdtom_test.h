/* Multidimensional sample arrays stored as k-d-trees. */
/* Last edited on 2021-07-13 01:07:06 by jstolfi */

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

void kdtom_test_show_box(FILE *wr, ppv_dim_t d, char *pref, ppv_index_t ixlo[], ppv_size_t size[], char *suff);
  /* Writes the elements {ixlo[0..d-1]} and {size[0..d-1]} as "0(10) 10(100) ..." 
    preceded by {pref} and followed by {suff}, if these strings are not {NULL}. */

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

void kdtom_test_clip(kdtom_t *T, int32_t loclip, int32_t hiclip);
  /* Tests {kdtom_clip(T,ixlo,size)} with various clip boxes.
    
    For each axis {ax}, as {loclip} ranges in {-2..+2}, the low index of
    the clip box will be just below {T.DK}, at the low end of {T.DK},
    somewhere inside {T.DK}, at the upper end of {T.DK}, and just above that end.
    
    The high index of the clip box along axis {ax} will be similarly
    selected as {hiclip} ranges in {-2..+2}.  Note that some combinations
    of {loclip} and {hiclip} may result in an empty box.
    
    Along axes other than {ax}, the clipbox will strictly contain 
    the core domain {T.DK}.
    
    Each clipped version {S} of {T} will be tested with
    {kdtom_test_get_sample}. The test requires {S.V[ix]} be the same as
    {T.V[ix]} if {ix} is in {S.DK}, and equal to {T.fill} otherwise. */

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
