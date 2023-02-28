#ifndef delaunay_H
#define delaunay_H

/* Delaunay triangulation. */
/* Last edited on 2023-02-20 06:11:02 by stolfi */

/* The divide-and-conquer algorithm for computing the Delaunay
  triangulation of a set of points (sites) on the plane. 
  For details, see 

    "Primitives for the Manipulation of General Subdivisions 
    and the Computation of Voronoi Diagrams"
    L. Guibas, J. Stolfi, ACM Transactions on Graphics, April 1985

  Implemented by Jim Roth (DEC CADM Advanced Group) on May 1986.
  Adapted by J. Stolfi in April 1993.
  See the copyright notice at the end of this file.
*/

#define _GNU_SOURCE
#include <stdint.h>

#include <quad.h>
#include <bool.h>
#include <sign.h>
#include <i3.h>
#include <hi2.h>
#include <r2.h>

#define delaunay_MAX_BITS_COORD (14)
  /* Number of bits used by a homogeneous site coordinate, excluding the sign. */

#define delaunay_MAX_COORD ((1 << delaunay_MAX_BITS_COORD)-1)
  /* Max absolute value of any homogeneous site coordinate. */

typedef struct delaunay_site_t
  { hi2_point_t pt;  /* Homogeneous site coordinates (cooordinate 0 is the weight). */
    int32_t index;   /* Site index. */
  } delaunay_site_t;

/* The Delaunay triangulation: */

quad_arc_t delaunay_build(delaunay_site_t st[], int32_t nsites);
  /* Builds the Delaunay trinagulation of {site[0..nsites-1]}.
    Returns an arc {e} on the perimeter of the triangulation,
    oriented so that {LEFT(e)} is the outer (unbounded) face. */

/* Quad-edge data pointers: */

#define ORG(e) ((delaunay_site_t *) quad_odata(e))
#define DST(e) ((delaunay_site_t *) quad_ddata(e))

#define SET_ORG(e,v) quad_set_odata(e, (void *)(v))
#define SET_DST(e,v) quad_set_ddata(e, (void *)(v))

/* Topology tools: */

quad_arc_t delaunay_make_edge(void);
  /* Creates an isolated non-loop edge with undefined {ORG} and {DST} fields. */

void delaunay_destroy_edge(quad_arc_t e);
  /* Disconnects the edge {e} from anything it may be connected to, and reclaims
    its storage. */

quad_arc_t delaunay_connect(quad_arc_t a, quad_arc_t b);
  /* Adds a new edge {e} connecting sites {DST(a)} to {ORG(b)},
    in such a way that {quad_onext(a) == e} and {quad_onext(b) == sym(e)}.
    If {a} and {b} are in distinct connected components, these will
    be joined and the two sides of {e} will be the same face.
    If {a} and {b} are in the same component, then their left faces
    must be the same, and the new edge {e} will split that face in two. */

/* Geometric tools: */

r2_t delaunay_r2_from_hi2(hi2_point_t *ph);
  /* Converts a point from homogeneous {int32_t} integer coordinates {h}
    to {double} Cartesian cooordinates. The result may have floating-point
    roundoff error.  Retuns {(INF,INF)} with arbitrary signs
    if the weight {h.c.c[0]} is zero. */

hi2_point_t delaunay_hi2_from_r2(r2_t *pc);
  /* Converts a point from {double} Cartesian cooordinates {p} to
    homogeneous {int32_t} integer coordinates whose absolute
    value does not exceed {delaunay_MAX_COORD}.  The maximum roundoff
    error is about {0.5/delaunay_MAX_COORD} in each coordinate. */

sign_t delaunay_orient(delaunay_site_t *a, delaunay_site_t *b, delaunay_site_t *c);
  /* Orientation of the three sites {a,b,c}: positive means CCW,
    negative means CW, zero means they are collinear. */
   
bool_t delaunay_leftof (delaunay_site_t *s, quad_arc_t e);
  /* TRUE iff the site {s} lies to the left of the line {ORG(e)-->DEST(e)}. */
    
bool_t delaunay_rightof (delaunay_site_t *s, quad_arc_t e);
  /* TRUE iff the site {s} lies to the right of the line {ORG(e)-->DEST(e)}. */
    
bool_t delaunay_incircle (delaunay_site_t *a, delaunay_site_t *b, delaunay_site_t *c, delaunay_site_t *d);
  /* TRUE iff the site {d} lies inside the circle {a,b,c}. */

#endif

/*
** Copyright notice:
**
** Copyright 1996 Institute of Computing, Unicamp.
**
** Permission to use this software for any purpose is hereby granted,
** provided that any substantial copy or mechanically derived version
** of this file that is made available to other parties is accompanied
** by this copyright notice in full, and is distributed under these same
** terms. 
**
** DISCLAIMER: This software is provided "as is" with no explicit or
** implicit warranty of any kind.  Neither the authors nor their
** employers can be held responsible for any losses or damages
** that might be attributed to its use.
**
** End of copyright notice.
*/
