/* Compact drawings of asexual descent trees. */
#ifndef drtree_compact_H
#define drtree_compact_H
/* Last edited on 2023-06-16 19:40:21 by stolfi */

#define drtree_H_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdint.h>

#include <drtree.h>

/* See {drtree_plot.h} for explanation of the plot grid. */

void drtree_compact_arrange
  ( int32_t ni, 
    drtree_node_t dt[], 
    int32_t tMin, 
    int32_t tMax, 
    int32_t rdr[],
    int32_t *ncols_P,
    int32_t *nrows_P
  );
  /* Assigns plot grid rows to all individuals described in
    {dt[0..ni-1]}, trying to minimize the number of rows used, without
    regard for planarity.
    
    If {q=dt[iq]} is a null node, sets {rdr[iq]} to {-1}. Otherwise
    requires that the life span of {q} be contained in {tMin..tMax}.
    Returns in {*ncols_P} the number of plot grid columns (which is
    {tMax-tMin+1}) and in {*nrows_P} the number of rows used (at most
    the number of non-null nodes, but hopefully much less). */

#endif

