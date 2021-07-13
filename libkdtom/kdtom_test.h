/* Multidimensional sample arrays stored as k-d-trees. */
/* Last edited on 2021-07-12 22:27:46 by jstolfi */

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

void kdtom_test_choose_array_size(ppv_dim_t d, ppv_size_t sz[]);
  /* Fills {sz[0..d-1]} with suitable values, not all equal. */

void kdtom_test_show_box(FILE *wr, ppv_dim_t d, char *pref, ppv_index_t ixlo[], ppv_size_t size[], char *suff);
  /* Writes the elements {ixlo[0..d-1]} and {size[0..d-1]} as "0(10) 10(100) ..." 
    preceded by {pref} and followed by {suff}, if these strings are not {NULL}. */

ppv_array_t *kdtom_test_make_array(ppv_dim_t d, ppv_size_t sz[], ppv_sample_t maxsmp);
  /* Creates an array {A} with the given attributes {d,sz[0..d-1],maxsmp},
    and fills it with a suitable test pattern. */

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
