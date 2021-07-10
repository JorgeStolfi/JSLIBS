/* Marching octahedra algorithm -- generating STL of isosurfaces. */
/* Last edited on 2021-07-07 22:56:11 by jstolfi */

#ifndef mroc_stl_H
#define mroc_stl_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <ppv_types.h>
#include <r3.h>
#include <i3.h>

#include <mroc.h>

void mroc_stl_classify_tetra_vertices(double f[], int32_t *niP, int32_t ki[], int32_t *noP, int32_t ko[]);
  /* Scans the values {f[0..3]} counting the numbers {ni,no} of corners
    of the tetrahedron that are respectively inside and outside the object.  Stores
    the indices of those corners in {ki[0..ni-1]} and {ko[0..no-1]}, and their 
    counts in {*niP,*noP}. */
    
void mroc_stl_write_surface_in_tetra
  ( FILE *wr, 
    r3_t p[], 
    double f[], 
    float eps,
    int32_t *ntP,
    int32_t *neP
  );
  /* Assumes {f[0..3]} are the occupancy values at the 
    four corners {p[0..3]} of a tetrahedron.
    If the tetrahedron straddles the 0.5 isosurface,
    writes the relevant part of that isosurface (one or two triangles)
    into {wr}. 
    
    The triangle corners are rounded to EVEN integer multiples
    of {eps}, and the triangle is discarded if two vertices coincide
    as a result.
    
    Adds to {*ntP} the number of triangles written out,
    and to {*neP} the number of triangles discarded. */

void mroc_stl_write_i3_triangle
  ( FILE *wr, 
    i3_t *p0, 
    i3_t *p1, 
    i3_t *p2,
    float eps,
    int32_t *ntP,
    int32_t *neP
  );
  /* Checks whether the three quantized points {p0,p1,p2}
    are all distinct.  If so, writes to {wr} the triangle with corners {p0,p1,p2},
    after multiplying all coordinates by {eps}, and increments {*ntP}.
    Otherwise writes nothing and increments {*neP}. */

void mroc_stl_write_r3_triangle
  ( FILE *wr, 
    r3_t *p0, 
    r3_t *p1, 
    r3_t *p2
  );
  /* Writes the triangle {p0,p1,p2} to [wr}. */
 
r3_t mroc_stl_unround_point(i3_t *p, float eps);
  /* Converts the quantized point {*p} to a non-quantized point,
    with coordinates in mm, by multiplying it by {eps}. */

r3_t mroc_stl_compute_normal(r3_t *p0, r3_t *p1, r3_t *p2);
  /* Computes the normal of the triangle with vertices {p0,p1,p2}. */
 
i3_t mroc_stl_edge_crossing(r3_t *p0, double f0, r3_t *p1, double f1, float eps);
  /* Assumes that {f0,f1} are the function values at points {p0,p1},
    one of them less than 0.5, the other one greater than 0.5. 
    Computes the point {q} along that segment where the function
    is exactly 0.5, assuming affine interpolation. Returns
    {q} with coordinates expressed as EVEN integer multiples of {eps},
    with proper rounding. */

/* DEBUGGING */

#define t2s_TETRA_SHRINK (0.8)
  /* Tetrahedron shrinking factor for {t2s_write_tetra_faces}. */

void mroc_stl_write_tetra_faces
  ( FILE *wr, 
    r3_t p[], 
    double f[], 
    int32_t *ntP
  );
  /* Assumes {f[0..3]} are the occupancy values at the 
    four corners {p[0..3]} of a tetrahedron.
    If the tetrahedron straddles the 0.5 isosurface,
    writes into {wr} the faces of the tetrahedron, 
    after shrinking it by some factor. 
    
    Adds to {*ntP} the number of triangles created. */

#endif
