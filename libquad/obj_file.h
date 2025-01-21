#ifndef obj_file_H
#define obj_file_H

/* Common definitions for OBJ file reading and writing. */ 
/* Last edited on 2025-01-09 23:24:44 by stolfi */

#define obj_file_H_copyright \
  "Copyright (C) 2024 Jorge Stolfi, UNICAMP.\n\n" jslibs_copyright

#include <stdio.h>
#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <r3.h>

vec_typedef(obj_file_face_vec_t,obj_file_face_vec,int32_vec_t);
  /* An extensible vector of extensible vectors of integers. */ 
  
typedef struct obj_file_data_t
  { /* Coordinate tables: */
    r3_vec_t V;              /* Vertices. */
    r3_vec_t T;              /* Texpoints. */
    r3_vec_t N;              /* Normals. */
    string_vec_t VL;         /* Vertex labels. */
    /* Face corner data: */
    obj_file_face_vec_t FV;  /* Vertex indices. */
    obj_file_face_vec_t FT;  /* Texpoint indices. */
    obj_file_face_vec_t FN;  /* Normal indices. */
  } obj_file_data_t;
  /* A set of tables containing the main information from 
    an OBJ file.
        
    The components {V}, {T}, and {N} are lists of points or vectors in {\RR^3}.
    
      {V.e[0..NV-1]} vertex coordinates (from "v" lines).
      {T.e[0..NT-1]} texpoint coordinates (from "vt" lines).
      {N.e[0..NN-1]} normal vectors (from "vn" lines).
      
    where {NV = V.ne}, {NT = T.ne}, and {NN = N.ne}.
      
    A /textpoint/ is a texture mapping point, a point in the domain of
    some 2D or 3D texture image.
    
    The vertex label {VL.e[k]} is the '#'-comments on the "v" line
    that defined the vertex {V.e[k]}, or {NULL} there was no comment
    there.  Each non-{NULL} label is separately allocated on the heap.
    
    The parameters {FV}, {FT}, and {FN} are tables that specify the
    attributes of each corner of each face. The number of faces {NF} is
    assumed to be {FV.ne} which should be equal to {FT.ne} and {FN.ne}.
    
    For each {i} in {0..NF-1}, element {FV.e[i]} is the list of the
    indices into {V.e} of the vertices on the border of face {i}.
    Likewise, {FT.e[i]} is the list of indices into {T.e} of the
    texpoints at those corners, and {FN.e[i]} is the list of indices
    into {N.e} of their normals.
    
    More specifically, the number of corners of face {i} is
    {FV.e[i].ne}, which must be equal to {FT.e[i].ne} and {FN.e[i].ne}.
    For each {j} in {0..FV.e[i].ne-1},
    
      {FV.e[i].e[j]} is the index into {V.e} of the coordinates of corner {j} of face {i}.
      {FT.e[i].e[j]} is the index into {T.e} of the texpoint of corner {j} of face {i}.
      {FN.e[i].e[j]} is the index into {N.e} of the normal of corner {j} of face {i}.
      
    The vertex index {FV.e[i].e[j]} must be in the range {0..NV-1}.
    The texpoint index {FT.e[i].e[j]} must be either in the range
    {0..NT-1}, or {-1} if no texpoint was specified for that corner.
    Similarly the normal index {FN.e[i].e[j]} must be either in the
    range {0..n-1}, or {-1} if no normal was specified for that
    corner.
    
    Note that these indices start at zero. In the OBJ file the indices
    start at 1, so vertex {V.e[i]} will be vertex number {i+1} in the
    file; and similary for the texpoint and normal indices. */

obj_file_data_t *obj_file_data_new(void);
  /* Creates a new data structure to represent the data in an OBJ file.
    
    The {V,T,F,VL} tables will be initially allocated with some arbitrary
    size, but will contain garbage. 
    
    The corner lists {FV}, {FT}, {FN} will be allocated with some
    arbitrary size, but their entries will be invalid {int32_vec_t}
    descriptors.
    
    The client must extend those tables as needed, initialize the
    individual face lists, fill them with the proper data, and trim them
    to the correct size. */
    
void obj_file_data_free(obj_file_data_t *D);
  /* Frees all space used by {D} and its tables, including the label
    strings. It assumes that all the tables have been initialized and
    trimmed, so any pointers in them are {NULL} point to existing heap
    areas. */

#define obj_file_data_prec_normal 7
  /* Number of decimal digits to use when writing coordinates of normal vectors. */
  
#define obj_file_data_tol_normal 0.0000002
  /* Tolerance to use when comparing normals after output and input. */

bool_t obj_file_data_compare(obj_file_data_t *D1, obj_file_data_t *D2, double tol, bool_t verbose);
  /* Compares the models {D1} and {D2}.  Returns {FALSE} iff the indices differ
    or the coordinates of vertices and texpoints differ by more than {tol},
    of the normal coordinates differ by more than {obj_file_data_tol_normal}.
    
    If {verbose} is true, prints messages explaining the differences. */

#endif
