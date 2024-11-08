/* Test tools for {hr2_pmap_opt.h} test programs  */
/* Last edited on 2024-11-04 06:39:36 by stolfi */

#ifndef hr2_pmap_opt_test_tools_H
#define hr2_pmap_opt_test_tools_H

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <jsrandom.h>
#include <affirm.h>

#include <rn.h>
#include <r3x3.h>
#include <r3.h>
#include <r2x2.h>
#include <r2.h>

#include <hr2.h>
#include <hr2_pmap.h>

void hr2_pmap_opt_test_choose_r2_point_pairs
  ( hr2_pmap_type_t type,
    sign_t sgn,
    bool_t tight,
    bool_t ident,
    double pert,
    bool_t verbose,
    int32_t *np_P,
    r2_t **p1_P,
    r2_t **p2_P,
    double **w_P,
    hr2_pmap_t *M_P
  );
  /* Chooses a suitable number {np} of point pairs {p1[i],p2[i]} and
    weights {w[i]}, with {i} in {0..np-1}, for testing the quadratic
    optimization with maps of the given {type} and {sgn}.
    
    If {tight}, {np} will be the minimum number {nr} of points needed to
    define such a map, and that map will be exact.
    
    The procedure will choose two maps {M1,M2} of the given {type} and
    {sgn}, then generate each pair of points {p1[i],p2[i]} by mapping
    the same random point {p0[i]} through these two maps, and adding to
    both random perturbations with max size {pert}. Thus, for all {i},
    {p1[i]} and {p2[i]} will be approximately related by the nominal map
    {M=M1^{-1} M2}.
    
    If {ident} is true, the nominal map {M} will be the identity if
    {sgn} is {+1}, or the {XY} swap map if {sgn} is {-1}. Otherwise {M}
    will be a random map of the given {type} and {sgn}.
    
    Returns {np,p1,p2,w} in {*np_P,*p1_P,*p2_P,*w_P}.  Sometimes,
    the weight table {w} will be {NULL}.
    
    Also returns in {*M_P} the nominal projective map {M=M1^{-1} M2}. */

void hr2_pmap_opt_test_print_map(char *name, hr2_pmap_t *M, double f2M);
  /* Prints the map {M} (direct and inverse) to {stderr},
    prefixed by a line"  {name} = ".  if {f2M} is not {NAN}
    prints it too. */

#endif
