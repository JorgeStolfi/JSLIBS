/* Handling of mesh slice polygons.  */
/* Last edited on 2024-10-05 07:39:28 by stolfi  */

#ifndef tosl_slice_H
#define tosl_slice_H

#define _GNU_SOURCE
#include <stdint.h>

#include <tosl.h>

typedef struct tosl_slice_t
  { int32_t NV;               /* Total vertex count of slice. */
    tosl_coord_t Z;       /* Z-coordinate of slice. */
    tosl_arc_id_t *iarc;  /* Arcs that determine the vertices of the slice. */
    int32_t NV_max;           /* Max vertex count of slice. */
  } tosl_slice_t;
  /* A {tosl_slice_t} record {S} stores a slice of the mesh by a
    horizontal plane. The slice consists of zero or more polygonal
    contours, with a total of {S.NV} vertices.

    For each {jv} in {0..S.NV-1}, vertex {jv} of the slice is the
    intersection of the plane at height {S.Z} with the arc {a} of the
    mesh whose {haf_arc_id} is {S.iarc[jv]}. The arc {a} must cross the
    slicing plane. The arc {a} should be directed downwards if that is
    the last vertex of a contour, and upwards otherwise. 
    
    The {iarc} vector has space for {NV_max} vertices. */

tosl_slice_t *tosl_slice_new(int32_t NV_max, tosl_coord_t Zp);
  /* Creates a new {tosl_slice_t} record {S} with space for {NV_max} polygon vertices.
    Sets {S.Z} to {Zp} and initalizes {S.NV} to 0. */

void tosl_slice_free(tosl_slice_t *S);
  /* Releases the storage of {S}, including internal tabels and the record {*S} itself. */

#endif
