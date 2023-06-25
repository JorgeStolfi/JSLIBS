/* Planar drawings of asexual descent trees. */
#ifndef drtree_planar_H
#define drtree_planar_H

/* Last edited on 2023-06-17 06:27:09 by stolfi */

#define drtree_planar_H_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdint.h>

#include <drtree.h>

void drtree_planar_arrange
  ( int32_t ni, 
    drtree_node_t dt[], 
    int32_t tMin, 
    int32_t tMax,
    int32_t rdr[],
    int32_t *ncols_P,
    int32_t *nrows_P
  );
  /* Assigns rows {rdr[0..ni-1]} to all individuals {dt[0..ni-1]}, using
    the corresponding visual information in {vf[0..ni-1]},
    trying to minimize {nrows} while keeping the layout planar. If a
    node {dt[iq]} is not visible, its assigned row is set to {-1}.
    
    Returns in {*ncols_P} and {*nrows_P} the number of columns and rows
    used.

    A layout is said to be /planar/ if, for any individual {q} with a
    parent {p} and visible birth, there is no other individual occupying
    any cell on the column {j = q.tbr-tMin} between {rdr(q)} and
    {rdr(p)}, except possibly for other chidlren of {p} with same time
    of birth as {q}. Then one can draw a line inside that column
    connecting the trace of {p} to the first cell of {q} and those other
    children, without crossing any other traces. */

#endif

