/* Last edited on 2025-04-27 11:30:40 by stolfi */
/* Geometry of a curved line that represents a graph edge in plots. */ 

#ifndef psgr_path_H
#define psgr_path_H

/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <i2.h>
#include <r2.h>
#include <bool.h>

typedef struct psgr_path_t
  { uint32_t n;     /* Number of internal vertices. */
    r2_t* v;        /* Coordinates of internal nodes (px). */
    bool_t reverse; /* TRUE if the nodes are to be taken in reversed order. */
  } psgr_path_t;
  /* A value {P} of this type is a descriptor of a directed polygonal
    line on the plane.
  
    The initial and final vertices {pa} and {pb} of the polygonal line
    are defined externally.  The other vertices (the /nodes/ of the path)
    are {P.v[0 .. .n-1]}. The coordinates of these nodes are in units of
    fractional pixels, and need to be scaled by the pixel size before
    being used for plotting.
    
    If {n=0} the implied polygonal line will be the straight line
    segment {S} from {pa} to {pb}.
  
    If {reverse} is false, the path nodes should be taken in the order
    {.v[0]} to {.v[n-1]}. If {reverse} is true, they should be taken in
    the reverse order, from {.v[n-1]} to {.v[0]}. Note that in this case
    the initial and final vertices {pa,pb} should be swapped
    too. */

#define psgr_path_NULL ((psgr_path_t){ .n = 0, .v = NULL, .reverse = FALSE })
  /* A null value for a {psgr_path_t} variable. */

#define psgr_path_MAX_NODES (4*4096)
  /* Should be far more than enough for most cases. */

psgr_path_t psgr_path_copy(psgr_path_t P);
  /* Returns a copy of path {P}, consisting of a copy of the node list 
    with the same {.reverse} bit. */

bool_t psgr_path_equal(psgr_path_t P0, psgr_path_t P1);
  /* Returns true iff paths {P0} and {P1} are equal: same node count {.n},
    same {.reverse} bit, and the coordinates of corresponding nodes
    are all equal. */

void psgr_path_free(psgr_path_t P);
  /* Reclaims the storage used by {P} (but not the descriptor {P} itself). */

psgr_path_t psgr_path_reverse(psgr_path_t P);
  /* Returns the path {r} which is the reverse of path {P}.  
    NOTE: {r} and {P} will share the same node list (that is, {r.v == P.v}),
    so at most one of them should be reclaimed. */

r2_t psgr_path_center(i2_t *ipo, psgr_path_t P, i2_t *ipd);
  /* Returns the nominal midpoint {ctr} of the path {P}.

    More precisely, if {P.n} is zero, {ctr} is the midpoint of {ipo} and {ipd};
    else, if {P.n} is odd, {ctr} is the middle node of {P.v};
    else {ctr} is the midpoint of the two middle nodes of {P.v}. */

r2_t psgr_path_central_dir(i2_t *ipo, psgr_path_t P, i2_t *ipd);
  /* Returns the a vector {dir} parallel to 
    the mean direction of the path around its nominal midpoint.

    More precisely, if {p.n} is 0 or 1 {dir} is {ipd-ipo};
    else, if {P.n} is odd, {dir} is the difference between the 
    nodes after and before the middle one; else, {ctr} the
    difference between the two middle nodes. */

r2_t psgr_path_starting_dir(i2_t *ipo, psgr_path_t P, i2_t *ipd);
  /* Returns the initial direction of the path {P}, assuming that 
    it starts at {ipo} and ends at {ipd}.   
    
    Specifically, if {P.n} is zero, returns the direction from {ipo} to {ipd}.
    Othwewise, if {P.reverse} is false, returns the direction from {ipo} to {P.v[0]}.
    Otherwise returns the direction from {ipo} to {P.v[P.n-1]}. */

psgr_path_t psgr_path_concatenate(psgr_path_t P0, r2_t *pm, psgr_path_t P1);
  /* Returns a path that is the concatenation of {P0} and {P1}, with 
    an extra node {pm} inserted between them.
    
    Typically, {P0} and {P1} describe the plot polygonals for two arcs, one
    that ends at a vertex with pixel coordinates {pm} and another
    that starts at that vertex.  However, {pm} may be any point in the domain
    {[0 _ NX-1] × [0 _ NY-1]} of the underlying image, even with fractional
    coordinates. */
 
psgr_path_t psgr_path_in_star_wedge(i2_t *ip1, i2_t *ipm, i2_t *ip2);
  /* Computes a {psgr_path_t} that is suitable for a new edge
    that will result from a star-cycle transformation.
    
    Specifically, given the pixel indices {ipm} of the star's central
    vertex and {ip1,ip2} of the destination vertices of those arcs,
    creates a path from {ip1} to {ip2} that follows close to the
    polygonal {ip1--ipm--ip2} but displaced towards the barycenter of
    those three vertices. */

psgr_path_t psgr_path_throw(i2_t *ip0, uint32_t n, i2_t *ip1);
  /* Generates a random path {P} that is supposed to start at {ip0}
    and end at {ip1}.  The path will have {n} internal vertices {P.v[0..n-1]},
    and {P.reverse} will be false. 
    
    If {n} is not zero, the points {ip0}, {P.v[0..n-1]}, and {ip1} will be {n+2}
    equally spaced points along the {S}, some of them with random
    displacements perpendicular to {S}. Specifically, if {n>=2}, points
    {P.v[0]} and {P.v[n-1]} (hence the first and last segments of the
    polygonal line) will lie on {S}. If {n=1}, the vertex {P.v[0]} will
    be displaced laterally too. */

void psgr_path_write(FILE* wr, psgr_path_t P);
  /* Writes to {wr} a decription of the path {P}, in a format compatible with
    {psgr_read_path}. */

psgr_path_t psgr_path_read(FILE* rd, uint32_t NX, uint32_t NY);
  /* Reads from {rd} a description of a {psgr_path_t}. It should be an integer {n},
    a reversal bit {rev} (0 or 1), and {n} pairs "({X},{Y})" where {X}
    and {Y} are {double} coordinates in fractional pixels, in the
    ranges {[0 _ NX-1]} and {[0 _ NY-1]} . */

#endif
