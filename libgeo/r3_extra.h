/* r3_extra.h --- additional operations on points and vectors of R^3 */
/* Last edited on 2014-01-12 12:43:14 by stolfilocal */

#ifndef r3_extra_H
#define r3_extra_H

#define _GNU_SOURCE

#include <r3.h>

double r3_pick_ortho (r3_t *u, r3_t *r);
  /* Sets {r} to an arbitrary vector orthogonal to {u},
    with the same L-infinity norm as {u} (which is the result of the call).
    In particular, if {u} is {(0,0,0)}, sets {r} to {(0,0,0)} and returns 0.
    Note that {r} is not a continuous function of {u}.  */

void r3_tetrahedron_vertices(double R, int n, r3_t r[]);
  /* Sets {r[0..n-1]} to the 4 corners of a regular tetrahedron with 
    center at the origin and circumradius {R}.  The mid-edge radius 
    will be {R/sqrt(3)}, the inradius will be {R/3}, and the side
    will be {R*sqrt(8/3)}. Requires {n==4}. */
  
void r3_octahedron_vertices(double R, int n, r3_t r[]);
  /* Sets {r[0..n-1]} to the 6 corners of a regular octahedron with 
    center at the origin and circumradius {R}.  The mid-edge radius
    will be {R/sqrt(2)}, the inradius will be {R/sqrt(3)}, and the 
    side will be {R*sqrt(2)}. Requires {n==6}. */
  
void r3_hexahedron_vertices(double R, int n, r3_t r[]);
  /* Sets {r[0..n-1]} to the 8 corners of a regular hexahedron (cube) with 
    center at the origin and circumradius {R}.  The mid-edge radius
    will be {R/sqrt(2)}, the inradius will be {R/sqrt(3)},
    and the side will be {2*R/sqrt(3)}. Requires {n==8}. */
  
void r3_icosahedron_vertices(double R, int n, r3_t r[]);
  /* Sets {r[0..n-1]} to the 12 corners of a regular icosahedron with 
    center at the origin and circumradius {R}.  The mid-edge radius
    will be {???}, the inradius will be {???}, and the side
    will be {R*sqrt(2 - 2/sqrt(5))}. Requires {n==12}. */
  
void r3_dodecahedron_vertices(double R, int n, r3_t r[]);
  /* Sets {r[0..n-1]} to the 20 corners of a regular dodecahedron with 
    center at the origin and circumradius {R}.  The mid-edge radius
    will be {???}, the inradius will be {???}, and the side
    will be {R*(sqrt(5)-1)/sqrt(3)}. Requires {n==20}. */

void r3_cylindrical_grid(double H, double R, int m, int n, double skew, r3_t r[]);
  /* Stores into {r[0..nr*nh-1]} the vertices of a cylindrical mesh centered
    at the origin whose axis is the {Z} axis, inscribed in a vertical cylinder
    with radius {R} and half-height {H}. The mesh has {m} horizontal rings
    each with {n} equally spaced vertices.  Each ring is twisted counterclockwise
    relative to the one below by {skew} times the angular spacing of its vertices.
    
    In particular, {m=1} gives a regular {n}-gon on the {Z=0} plane;
    {m=2,skew=0} gives the right {n}-prism with regular bases;
    and {m=2,skew=0.5} gives the {n}-antiprism with regular bases. */

#endif
