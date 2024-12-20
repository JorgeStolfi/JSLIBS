#ifndef delaunay_debug_H
#define delaunay_debug_H

/* Debugging tools for {delaunay.c}. */
/* Last edited on 2024-12-05 10:25:03 by stolfi */

#include <stdint.h>
 
#include <quad.h>
#include <bool.h>
#include <sign.h>

#include <delaunay.h>

#define DEBUGGING_QUAD FALSE
  /* Define as TRUE to check consistency of the quad-edge data structure. */

#define DEBUGGING_DELAUNAY FALSE
  /* Define as TRUE to check topology and geometry of the Delaunay diagram. */

/* Face data pointers (only for debugging): */

typedef struct delaunay_face_t { } delaunay_face_t;

#define LEFT(e) ((delaunay_face_t *) quad_ldata(e))
#define RITE(e) ((delaunay_face_t *) quad_rdata(e))

#define SET_LEFT(e,F) quad_set_ldata(e, (void *)(F))
#define SET_RITE(e,F) quad_set_rdata(e, (void *)(F))

void deldebug_print_site (char *msg, delaunay_site_t *a);
  /* Prints site {a} to {stderr} labeled with {msg}. */

void deldebug_print_edge (char *msg, quad_arc_t e);
  /* Prints arc {e} to {stderr} labeled with {msg}. */

void deldebug_check_quad_edge(quad_arc_t e);
  /* Performs consistency checks on the fields of {quad-edge(e)}. */

void deldebug_check_quad_all_edges(void);
  /* Applies {delaunay_check_quad_edge} on all edges created so far. */

quad_arc_t deldebug_make_edge(void);
  /* Same as {quad_make_edge} but also saves the edge for checking. */

void deldebug_destroy_edge(quad_arc_t e);
  /* Dymmy {quad_destroy_edge} that does not actually free {e}. */

void deldebug_check_delaunay_edge_property(quad_arc_t e);
  /* Checks the Delaunay diagonal property for edge {e}, assuming that
    it is an internal edge. */

void deldebug_check_interior_face(quad_arc_t e);
  /* Checks whether the face {LEFT(e)} has three or more sides and
    is convex. If it has more than 3 sides, requires that all sites
    be cocircular. */

void deldebug_set_all_left_faces(quad_arc_t e, delaunay_face_t *F);
  /* Sets {F} as the {LEFT} pointers of all arcs with the same left face as {e}. */

void deldebug_check_left_triangle(quad_arc_t b, int32_t depth);
  /* Checks the triangle at the left of edge {e}, assuming that
    it is an internal triangle. If {depth > 0}, also 
    looks at the two other edges {e,f} of the triangle, oriented 
    clockwise. If {e} is an internal edge, tests its Delaunay 
    property and then calls recursively {delaunay_check_left_triangle(e,depth-1)}.
    Ditto for {f}. */

bool_t deldebug_has_left_triangle(quad_arc_t e);
  /* Returns TRUE iff {LEFT(e)} is a counterclockwise triangle. */

#endif
