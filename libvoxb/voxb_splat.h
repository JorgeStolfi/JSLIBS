/* voxb_splat.h --- voxel-based modeling of 3D objects (binary tomogram version) */
/* Last edited on 2021-06-14 20:55:38 by jstolfi */

#ifndef voxb_splat_H
#define voxb_splat_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <r3_motion.h>
#include <ppv_array.h>

/* VOXEL-BASED SOLID MODELING 

  The procedures in this interface generate models of 3D objects by
  "painting" or "splatting", them into a binary tomogram -- a 3D array of
  binary voxels.
  
  The array must have non-trivial size along axes 0, 1, and 2, which are
  identified with the coordinate axes {Z}, {Y} and {X}, in that order.
  That is, the index 0 selects an horizontal layer of the array, index 1
  selects a row in that layer parallel to the {X} axis, and index 2
  selects a voxel in that row.
  
  This module represents shapes as union of voxels. For the purpose of
  this module, each voxel is imagined as a cube of unit size with
  integer corner coordinates. If an array has {N=a.size[j]} voxels along
  some axis {j}, the voxels are supposed to span the coordinate interval
  {[0 _ N]} along that axis. The index {i=ix[j]} of a voxel along that
  axis ranges from 0 to {N-1}; the voxel spans the range {[i _ i+1]} on
  that axis, and its center is at {i+0.5}.
  
  Voxel values in a binary tomogram {a} (/samples/) are booleans, either 0 (false)
  meaning 'exterior', or 1 (true) meaning 'interior'.  */

/* SINGLE VOXEL SPLATTING */

typedef enum {
    voxb_op_OR,  /* a = a | v    (union of shapes). */
    voxb_op_AND, /* a = a & v    (intersecion). */
    voxb_op_SUB, /* a = a & (~v) (subtraction). */
    voxb_op_XOR, /* a = a XOR v  (symm. difference). */
    voxb_op_LIM  /* Limit of enum values. */
  } voxb_op_t; 
  /* Specifies the operation to be performed by {voxb_splat_voxel} and the like. */
 
void voxb_splat_voxel(ppv_array_desc_t *A, int32_t kx, int32_t ky, int32_t kz, bool_t val, voxb_op_t op);
 /* Modifies the stored value {oldv} of voxel {kx} of row {ky} of layer {kz} of {A} 
   with the value {val}. Specifically, replaces the old value {ak} by {ak op val}.
   
   Note that {kx,ky,kz} are indices {2,1,0} of {A}, in that order. */

/* SPLATTING OBJECTS FROM POINT PREDICATES */
  
void voxb_splat_object
  ( ppv_array_desc_t *A,
    r3_pred_t *obj,
    r3_motion_state_t *S,
    double maxR,
    voxb_op_t op
  );
  /* Splats into the voxel array {A} the primitive object defined by the
    occupancy function {obj}, modified by the matrix {S.M} and
    translated by {S.p}. The matrix {S.M} must be an isometry (a
    rotation or a reflection).
    
    The predicate {obj(&p)} receives a point of the array's domain,
    translated by minus {S.p} and transformed by the investe of {S.M}.
    The coordinates are in the geometric order {(X,Y,Z)}, reversed from
    the voxel index order. The predicate should return {TRUE} if the
    given point {p} is inside the object {obj}, {FALSE} if outside.
    
    Assumes that the modified and translated object has zero occupancy
    at any point that differs more than {maxR} units from {S.p} along
    any axes. */

void voxb_splat_object_multi
  ( ppv_array_desc_t *A,
    r3_pred_t *obj,
    int32_t ns,
    r3_motion_state_t S[],
    double maxR,
    voxb_op_t op
  );
  /* Splats into the voxel array {A} the object {obj} modified by the
    matrces {S[k].M} and translated by {S[k].p}, for {k} in {0..ns-1}.
    Assumes that the modified object extends at most {maxR} units from
    its reference point, for any {k}. */

/* HACKS */

extern bool_t voxb_splat_debug;
  /* Set to true to see all voxels being splatted. */

#endif

