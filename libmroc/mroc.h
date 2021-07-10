/* Marching octahedra algorithm. */
/* Last edited on 2021-07-07 17:56:36 by jstolfi */

#ifndef mroc_H
#define mroc_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <ppv_types.h>
#include <r3.h>

/* 
  THE MARCHING OCTAHEDRA ALGORITHM
  
  The marching octahedra algorithm is a systematic way of exploring a
  real-valued trivariate function {F} that is sampled on a regular grid
  of points, and interpolated linearly to the rest of some region {B} of
  space.
  
  The region of interest {B} is assumed to be a three-dimensional box
  divided into a regular orthogonal grid of /cells/ or /voxels/. The
  sampling points are assumed to be the centers and the corners of those
  cells. (The function can be sampled at only one of these two sets of
  points: the vertex values can be estimated as the average of the 8
  adjacent cell center values, or vice-versa.)
  
  The algorithm is a variant of the marching cubes algorithm. It tiles
  the region of interest with a mesh of octahedra, each defined by the
  centers of two adjacent cells of the cubical voxel mesh and the four
  corners of their shared square face. Then its splits each octahedron
  into 4 tetrahedra with a shared edge connecting the two voxel centers.
  
  The result is a mesh of tetrahedra that covers the box {B}. Each
  tetrahedron has two corners that are voxel centers and two corners
  that are voxel corners. The function {F} is assumed to be affine
  ("linear") in each tetrahedron, i.e. obtained by interpolation of the
  four corner samples.
  
  The marching octahedra algorthm simply sweeps this mesh in a
  systematic order: roughly, layer by layer in the {Z} direction, then
  row by row in the {Y} direction inside each layer, and then voxel by
  voxel in the {X} direction inside each row. The
  processing of each tetrahedron is specified by user-defined procedural
  parameter.
  
  Note that each tetrahedron spans two adjacent voxels, and each voxel
  is overlapped by 24 tetrahedra. Therefore the mesh has approximately
  12 tetrahedra for each voxel (apart from boundary adjustments).
  
  More precisely, there is one octahedron for each voxel face in the
  voxel mesh, including faces that are on the boundary of {B}. Thus, if
  the voxel mesh has {NX}, {NY}, and {NZ} cells along the {X}, {Y}, and
  {Z} directions, the number of octahedra is {(NX+1)*NY*NZ +
  NX*(NY+1)*NZ + NX*NY*(NZ+1)}, that is, {3*NX*NY*NZ + NY*NZ + NX*NZ +
  NX*NY}; and the number of tetrahedra is four times that.
  
  The voxel faces that are on the boundary of {B} define octahedra that
  lie partly outside {B}. In order to process those octahedra, the
  algorithm requires extra samples at 'dummy' voxels tha are outside {B}
  but adjacent to it. The client of the algorithm must provide suitable
  'padding' function values for the centers of those dummy voxels. 
  
  The algorithm is defined here as two procedures {mroc_2D_inter} and
  {mroc_2D_intra}, that together process only one layer of the voxel
  array. The client must iterate them to sweep the whole region {B},
  and provide suitable padding values for the first and last layer. */

typedef void mroc_tetra_proc_t(r3_t *p[], double f[]); 
  /* Type of a procedure that processes one tetrahedron
    with corners {*(p[0..3])} and corresponding function 
    values {f[0..3]}.
    
    Points {p[0]} and {p[1]} are voxel centers (with integer-plus-half
    coordinates), while {p[2]} and {p[3]} are voxel corners (with
    integer coordinates). */

void mroc_2D_inter
  ( int32_t iz,
    ppv_size_t NCX,
    ppv_size_t NCY, 
    double *antcL, 
    double *poscL, 
    double *midvL,
    mroc_tetra_proc_t tetra_proc  /* Tetrahedron processing function. */
  );
  /* Processes the tetrahedra that span the consecutive voxel layers
    {iz-1} and {iz}, assumed to have {NCY} rows of {NCX} voxels.
    
    The cell-centered function sample values are supposed to be stored in {antcL}
    and {poscL}, both with {NCY+2} rows of {NCX+2} elements (including
    dummy padding voxels along each border, which are ignored). Assumes
    that {midvL} is an array with {NCY+1} rows of {NCX+1} elements
    containing the function sample values at the vertices between those two
    voxel layers. */

void mroc_2D_intra
  ( int32_t iz,
    ppv_size_t NCX,
    ppv_size_t NCY, 
    double *thiscL, 
    double *botvL, 
    double *topvL,
    mroc_tetra_proc_t tetra_proc  /* Tetrahedron processing function. */
  );
  /* Processes the tetrahedra that span adjacent voxels 
    inside layer {iz}, whose cell-centered function sample values are stored in
    the array {thiscL} with {NCY+2} rows of {NCX+2} elements
    (including padding border elements, which are used).  Assumes that {botvL} 
    and {topvL} are arrays with {NCY+1} rows of {NCX+1} elements,
    containing the function sample values at the 
    vertices respectively below and above that voxel layer. */
  
void mroc_process_octahedron
  ( r3_t pu[], 
    double fu[], 
    r3_t pv[], 
    double fv[],
    mroc_tetra_proc_t tetra_proc  /* Tetrahedron processing function. */
  );
  /* Assumes {pu,fu} have 4 elements and {pv,fv} have 2 elements,
    and each element of {fu,fv} is the function value at the 
    corresponding point {pu,pv}.
    
    Processes the 4 tetrahedra whose corners are {pv[0],pv[1]},
    and each 2 succesive elements of {pu[0..3]}.  */

#endif
