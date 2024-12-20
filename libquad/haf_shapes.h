#ifndef haf_shapes_H
#define haf_shapes_H

/* Procedures that build half-edge structures for various simple maps. */
/* Last edited on 2024-12-05 10:39:06 by stolfi */

#define haf_shapes_H_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP)\n\n" jslibs_copyright
  
#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>

#include <haf.h>

haf_arc_t haf_shapes_torus(void);
  /* A map on the torus, with two edges and one face. */

haf_arc_t haf_shapes_bitorus(void);
  /* A map on the bitorus. */

haf_arc_t haf_shapes_tritorus(void);
  /* A map on the tritorus, with tetrahedral symmetry. */

haf_arc_t haf_shapes_star(int32_t n);
  /* A map on the sphere with one face, one vertex, and
    {n} edges radiating from that vertex. */

haf_arc_t haf_shapes_ring(int32_t n);
  /* A map on the sphere with two faces separated by
    a ring with {n} edges and {n} vertices. */

haf_arc_t haf_shapes_pyramid(int32_t n);
  /* A map on the sphere with one base face with {n} sides,
    adjacent to {n} triangles that share a common vertex. */

haf_arc_t haf_shapes_orange(int32_t n);
  /* Builds a map on the sphere with two vertices (`poles') connected
    by {n} edges (`meridians'), that divide the sphere
    into {n} wedge-like faces. Namely, the dual of {haf_shapes_ring}. */

haf_arc_t haf_shapes_tetra(void);
  /* The map of a tetrahedron. */
  
haf_arc_t haf_shapes_cube (void);
  /* The map of a cube. */

void haf_shapes_buld_tower(int32_t m, int32_t h, haf_arc_t a);
  /* Builds a cylindrical tower on the face {haf_right(a)}, 
    which must have {m} edges. The tower will have {h} stages
    and a roof; each stage will be a ring with {m} square faces.
    The roof of the tower will be a single face of {m} edges. */

#endif
