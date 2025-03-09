#ifndef delaunay_grid_H
#define delaunay_grid_H

/* Grids from delaunay meshes. */
/* Last edited on 2025-02-16 19:32:42 by stolfi */

#include <stdint.h>

#include <quad.h>
#include <bool.h>
#include <sign.h>
#include <i3.h>
#include <hi2.h>
#include <r2.h>
#include <delaunay.h>

quad_arc_t delaunay_grid_build(int32_t *ns_P, delaunay_site_t **site_P);
  /* On input, {*site_P} must be a vector {site[0..ns-1]}, where {ns} is 
    the input value of {*ns_P}.  The procedure first builds a Delaunay triangulation
    of those sites.  Then it looks for pairs of adjacent triangles that are 
    close to a rectangle, and replaces each pair by four triangles 
    with an extra vertex at the midpoint of their common edge.
    Appends those new vertices at the end of the {site} array.
    reallocating it as needed.  Returns its new address in {*site_P},
    and the new count in {*ns_P}. */

#endif
