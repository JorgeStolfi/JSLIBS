#ifndef cpk_build_H
#define cpk_build_H

/* Building geographic incompatibility graphs. */
/* Last edited on 2024-12-31 14:12:54 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <bool.h>
#include <r2.h>
#include <vec.h>

#include <cpk_basic.h>
#include <cpk_main.h>
#include <hxg_canvas.h>

/* THE STATIOW PLACEMENT GRAPH */

void cpk_build_graph
  ( cpk_domain_t *C,     /* Geographic data. */
    cpk_policy_t *P,     /* Location planning parameters. */
    r2_vec_t *V,         /* OUT: List of candidate points. */
    double_vec_t *W,     /* OUT: Demand coverage of each candidate auction point. */
    ui2_vec_t *E,        /* OUT: Proximity (= incompatibility) edges. */
    bool_t verbose       /* TRUE prints diagnostic messages. */
  );
  /* Builds the graph {V,E} that describes the candidates for auction
    centers and the pairs of candidates that are incompatible because
    they are too close.  Also computes a `desirability' weight {W[v]}
    for each vertex {v}. 
    
    The vertices {V} are a finite subset of all the points of
    {\Real^2} that are valid according to the domain description {C}.
    The points are taken from a regular hexagonal grid selected by the
    procedure. If {P->for_demand} is TRUE, the set {V} has only the
    {C}-valid points that lie within distance {P->rAuc} from the
    declared demand points {P->DI[0..]}. Otherwise {V} has a dense
    grid of points, covering all the {C}-valid region of the plane.
    
    In any case, the edges {E} of the graph are all pairs {i,j} such
    that {i<j} and {dist(V[i],V[j]) < P->dMin}.
    
    The vectors {V}, {W}, and {E} are allocated by the procedure. */

#endif
