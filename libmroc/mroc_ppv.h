/* Applies the marching octahedra algorithm to a 3D sample array (a {ppv_array_t}). */
/* Last edited on 2021-07-08 14:35:05 by jstolfi */

#ifndef mroc_ppv_H
#define mroc_ppv_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <ppv_types.h>
#include <ppv_array.h>
#include <r3.h>

#include <mroc.h>

/* !!! Provide option for toroidal array domain. !!! */

void mroc_ppv_process_array
  ( ppv_array_t *A,
    double pad,
    mroc_tetra_proc_t tetra_proc
  );
  /* Applies the marching octahedra algorithm to the 
    voxel array {A}, calling {tetra_proc} on each tetrahedron of the 
    algorithm's mesh.
    
    The array {A} must be non-empty. Axes 0, 1, and 2 are assumed to be
    {Z}, {Y}, and {X}, respectively. The array must be trivial ({size =
    1}) along other axes.
    
    Each voxel is assumed to be a cubical mesh cell with unit side. If
    the array has sizes {NX,NY,NZ} along axes {2,1,0}, then its voxels
    are assumed to span the box {B = [0 _ NX] × [0 _ NY] × [0 _ NZ]} of
    {\RR^3}.
    
    Each sample {smp} of {A} is converted to a value between 0 and 1
    with {mroc_ppv_floatize(smp,A->maxsmp)}, and is assigned to the
    corresponding voxel center. The values at voxel corners are obtained
    by averaging the center values of the 8 surrounding voxels.
    
    In order to provide vertex values along the boundary of {B}, and
    voxel center values for octahedra that span that boundary, the array
    {A} is implicitly surrounded by a layer of voxels with center values
    set to the {pad} parameter. However, tetrahedra that are entirely
    inside this padding layer are not visited. */

double mroc_ppv_floatize(ppv_sample_t smp, ppv_sample_t maxsmp);
  /* Converts a sample {smp} from integer to {double} by the affine map
    that takes 0 and {maxsmp} to 0 and 1, respecively. Fails if {smp >
    maxsmp}. Returns {NAN} if {maxsmp} is zero. The result is never {0.5}. */

#endif
