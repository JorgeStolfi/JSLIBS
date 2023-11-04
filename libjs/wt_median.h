#ifndef wt_median_H
#define wt_median_H

/* Weight tables for filtering digital signals */
/* Last edited on 2023-11-03 19:34:57 by stolfi */

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

double wt_median(int32_t nx, double x[], int32_t ix, int32_t nw, double w[], bool_t interp, int32_t *nk_P, int32_t kx[]);
  /* Computes the weighted median of {nw} consecutive elements of {x[0..nx-1]}
    with window weights {w[0..nw-1]}.  
    
    The window width {nw} must be odd, {nw = 2*hw+1}. The window weights
    must be positive and must add to 1. The elements to be processed are
    {X = x[ix-hw..ix+hw]}. These elements must exist; that is, we must have
    {ix-hw >= 0} and {ix+hw < nx}.
    
    Each element {x[ix+k]} is assumed to have weight {w[hw+k]}, for {k}
    in {-hw..+hw]}. Roughly, the procedure returns a value {xm} such
    that the sum {Slo} of the weights of elements with {x[ix+k] <= xm}
    and the the total weight {Shi} of the elements with {x[ix+k] >= xm}
    are as close as possible.
    
    More precisely, if there is an element {xm} in the set {X}
    with that property, the procedure returns that element.
    The {interp} flag is ignored in this case.
   
    If there is no such element, then there must be two elements {xa}
    and {xb} in {X} such that {xa < xb}, there is no other element {xc}
    in {X} with {xa < xc < xb}, the sum {Sb} of weights of elements
    {x[ix+k] <= xb} is at more than 0.5, and the sum {Sa} of all weights
    of elements {x[ix+k] >= xa} is more than 0.5. In that case, if
    {interp} is false, the procedure returns either {xa} or {xb},
    depending on which choice will make the difference {Slo-Shi} closer
    to zero, with a suitable tie-breaking criterion. If {interp} is
    true, the procedure will return some weighted average of {xa} and
    {xb}.
    
    If {kx} is not {NULL}, it should point to a vector of at least {nw}
    elements that will be used to speed up the search. In that case, on
    entry {*nk_P} must be an integer {nk} in {0..nx}, and {kx} must
    point to a vector with at least {max(nw,nk)} elements. Elements
    {kx[0..nk-1]} must be a set of {nk} consecutive indices in
    {0..nx-1}. It is okay if {nk} is zero. On exit, {*nk_P} will have been
    set to {nw}, and {kx[0..nw-1]} will be the indices {ix-hw..ix+hw} of the
    window elements, sorted so that {x[kx[0..nw-1]]} are in increasing
    order. */

#endif
