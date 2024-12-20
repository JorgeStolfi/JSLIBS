#ifndef oct_shapes_H
#define oct_shapes_H

/* Procedures that build oct-edge structures for various simple maps. */
/* Last edited on 2024-12-05 10:39:40 by stolfi */

#define oct_shapes_H_copyright \
  "Copyright © 1996, 2006 State University of Campinas (UNICAMP)\n\n" jslibs_copyright
  
#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>

#include <oct.h>

void buld_tower(uint m, uint h, oct_arc_t a);
  /* Builds a cylindrical tower on the face {oct_right(a)}, 
    which must have {m} edges. The tower will have {h} stages
    and a roof; each stage will be a ring with {m} square faces.
    The roof of the tower will be a single face of {m} edges. */

oct_arc_t make_ring (int32_t n);
  /* Builds a map on the sphere with two faces separated by
    a ring with {n} edges and {n} vertices. */

oct_arc_t make_orange(uint n);
  /* Builds a map on the sphere with two vertices (`poles') connected
    by {n} edges (`meridians'), that divide the sphere
    into {n} wedge-like faces. Namely, the dual of {make_ring}. */

oct_arc_t make_tetra(void);
  /* The map of a tetrahedron. */
  
oct_arc_t make_cube (void);
  /* The map of a cube. */
  
oct_arc_t make_stick(void);
  /* A map on the sphere with a single edge, two vertices,
    and a single face. Namely, the dual of {make_edge}. */

oct_arc_t make_sausage (uint len);
  /* A sausage-like map made by taking a strip of {2 × len} squares,
   sewing its long sides together, and closing the ends with two digons. */
  
oct_arc_t make_fork(uint np, uint len);
  /* A star-shaped object consisting of {np} `prongs' radiating
    from a `hub'.  The hub is the result of {make_orange(np)},
    and each prong is similar to the result of {make_sausage(len)},
    minus one of the end-caps. Each face of the hub is deleted and 
    its border is sewed to the open end of the prong. */

oct_arc_t make_projetive (void);
  /* A map on the projetive plane, obtained by taking a map with a 
    single non-look edge and splicing the ends of that edge, with
    a twist. */
  
oct_arc_t make_torus(void);
  /* A map on the torus, with two edges and one face. */

oct_arc_t make_bitorus(void);
  /* A map on the bitorus. */

oct_arc_t make_tritorus(void);
  /* A map on the tritorus, with tetrahedral symmetry. */

oct_arc_t make_klein(void);
  /* A map on the Klein bottle, with two edges and one face. */

oct_arc_t make_klein2(void);
  /* A map on the Klein bottle, with four edges. */

oct_arc_t make_klein3(void);
  /* A map on the Klein bottle, with six edges. */

#endif
