/* voxm_splat.h --- voxel-based modeling of antialiased 3D objects */
/* Last edited on 2021-06-22 13:48:13 by jstolfi */

#ifndef voxm_splat_H
#define voxm_splat_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <r3_motion.h>
#include <ppv_array.h>

/* VOXEL-BASED SOLID MODELING 

  The procedures in this interface generate models of 3D objects by
  "painting" or "splatting", them into a tomogram -- a 3D array of
  voxels.
  
  The array must have non-trivial size along axes 0, 1, and 2, which are
  identified with the coordinate axes {Z}, {Y} and {X}, in that order.
  That is, the index 0 selects an horizontal layer of the array, index 1
  selects a row in that layer parallel to the {X} axis, and index 2
  selects a voxel in that row.
  
  For the purpose of this module, each voxel is imagined as a cube of
  unit size with integer corner coordinates. If an array has {N=a.size[j]}
  voxels along some axis {j}, the voxels are supposed to span
  the coordinate interval {[0 _ N]} along that axis. The index {i=ix[j]}
  of a voxel along that axis ranges from 0 to {N-1}; the voxel spans the
  range {[i _ i+1]} on that axis, and its center is at {i+0.5}.
  
  Voxel values in a tomogram {a} (/samples/) are unsigned integers from
  0 to {MAXVAL = 2^a.bps}. They are interpreted as fractional
  (/occupancy values/) between 0.0 and 1.0; where voxel value 0 means
  occupancy 0.0, {MAXVAL} means 1.0, and other voxel values interpolate
  linearly between them.
  
  Futhermore, the occupancy values 0.0 and 1.0 are assumed to mean that
  the center of the voxel is well outside or well inside the object,
  respectively. Intermetiate values are used in a /fuzzy layer/ near the
  object boundary to improve the smoothness of its implied surface,
  which is asumed to be the locus of points with occupancy 0.5.
  
  Splatting overlapping objects will produce the union of their definite
  interiors (voxels with value {MAXVAL}) and the intersection of their
  definite exteriors (voxels with value 0). Intermediate values will
  also be correct in voxels of the result where the original contents or
  the splatted item were 0.  However, in voxels where both had
  intermediate values, the resulting intermediate value will generally
  be incorrect.   */

/* SINGLE VOXEL SPLATTING */
 
void voxm_splat_voxel(ppv_array_t *A, int32_t kx, int32_t ky, int32_t kz, double val, bool_t sub);
 /* Modifies the stored value {oldv} of voxel {kx} of row {ky} of layer {kz} of {A} 
   with the value {val}, after quantizing {val} to an unsigned integer {newv}.
   
   Specifically, if {sub} is false, stores into the voxel the maximum of
   {oldv} and {newv}. If {sub} is true, stores instead the minimum of
   {oldv} and {MAXVAL-newv}.
   
   Note that {kx,ky,kz} are indices {2,1,0} of {A}, in that order. */

/* SPLATTING OBJECTS FROM OCCUPANCY FUNCTIONS */

typedef double voxm_splat_obfun_t(r3_t *p);
  /* Type of the occupancy function of a fuzzy object, that 
    returns 0 if {p} is well outside the object, 1 if it
    is well inside it, and fractional values in the fuzzy layer. */
   
void voxm_splat_object
  ( ppv_array_t *A,
    r3_double_func_t *obj,
    r3_motion_state_t *S,
    double maxR,
    bool_t sub
  );
  /* Splats into the voxel array {A} the object defined by the occupancy function
    {obj}, modified by the matrix {S.M} and translated by {S.p}.
    The matrix {S.M} must be an isometry (a rotation or a reflection).
    
    The function {obj} should return 0 if {p} is well outside the
    object, 1 if it is well inside it, and fractional values in the
    fuzzy layer.
    
    Assumes that the modified and translated object has zero occupancy
    at any point that differs more than {maxR} units from {S.p} along
    any axes. Note that this range must completely include any fuzzy
    layer. */

void voxm_splat_object_multi
  ( ppv_array_t *A,
    r3_double_func_t *obj,
    int32_t ns,
    r3_motion_state_t S[],
    double maxR,
    bool_t sub
  );
  /* Splats into the voxel array {A} the object {obj} modified by the
    matrces {S[k].M} and translated by {S[k].p}, for {k} in {0..ns-1}.
    Assumes that the modified object extends at most {maxR} units from
    its reference point, for any {k}. */

#endif

