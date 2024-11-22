/* r3_hedron.h --- building special polyhedra */
/* Last edited on 2024-11-20 12:55:48 by stolfi */

#ifndef r3_hedron_H
#define r3_hedron_H

#define _GNU_SOURCE
#include <stdint.h>

#include <r3.h>

void r3_hedron_tetra_vertices(double R, uint32_t n, r3_t r[]);
  /* Sets {r[0..n-1]} to the 4 corners of a regular tetrahedron with 
    center at the origin and circumradius {R}.  The mid-edge radius 
    will be {R/sqrt(3)}, the inradius will be {R/3}, and the side
    will be {R*sqrt(8/3)}. Requires {n==4}. */
  
void r3_hedron_octa_vertices(double R, uint32_t n, r3_t r[]);
  /* Sets {r[0..n-1]} to the 6 corners of a regular octahedron with 
    center at the origin and circumradius {R}.  The mid-edge radius
    will be {R/sqrt(2)}, the inradius will be {R/sqrt(3)}, and the 
    side will be {R*sqrt(2)}. Requires {n==6}. */
  
void r3_hedron_hexa_vertices(double R, uint32_t n, r3_t r[]);
  /* Sets {r[0..n-1]} to the 8 corners of a regular hexahedron (cube) with 
    center at the origin and circumradius {R}.  The mid-edge radius
    will be {R/sqrt(2)}, the inradius will be {R/sqrt(3)},
    and the side will be {2*R/sqrt(3)}. Requires {n==8}. */
  
void r3_hedron_icosa_vertices(double R, uint32_t n, r3_t r[]);
  /* Sets {r[0..n-1]} to the 12 corners of a regular icosahedron with 
    center at the origin and circumradius {R}.  The mid-edge radius
    will be {???}, the inradius will be {???}, and the side
    will be {R*sqrt(2 - 2/sqrt(5))}. Requires {n==12}. */
  
void r3_hedron_dodeca_vertices(double R, uint32_t n, r3_t r[]);
  /* Sets {r[0..n-1]} to the 20 corners of a regular dodecahedron with 
    center at the origin and circumradius {R}.  The mid-edge radius
    will be {???}, the inradius will be {???}, and the side
    will be {R*(sqrt(5)-1)/sqrt(3)}. Requires {n==20}. */

void r3_hedron_cylinder(double H, double R, uint32_t m, uint32_t n, double skew, r3_t r[]);
  /* Stores into {r[0..nr*nh-1]} the vertices of a cylindrical mesh centered
    at the origin whose axis is the {Z} axis, inscribed in a vertical cylinder
    with radius {R} and half-height {H}. The mesh has {m} horizontal rings
    each with {n} equally spaced vertices.  Each ring is twisted counterclockwise
    relative to the one below by {skew} times the angular spacing of its vertices.
    
    In particular, {m=1} gives a regular {n}-gon on the {Z=0} plane;
    {m=2,skew=0} gives the right {n}-prism with regular bases;
    and {m=2,skew=0.5} gives the {n}-antiprism with regular bases. */

#endif
