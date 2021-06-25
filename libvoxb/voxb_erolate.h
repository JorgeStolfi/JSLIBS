/* voxb_erolate.h --- erosion and dilation of PPV arrays. */
/* Last edited on 2021-06-22 13:46:43 by jstolfi */

#ifndef voxb_erolate_H
#define voxb_erolate_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <ppv_array.h>
#include <ppv_brush.h>

/* OPERATIONS WITH BRUSHES */

void voxb_erolate_with_brush(ppv_array_t *A, ppv_brush_t *b, bool_t erode);
 /* Erodes or dilates the shape stored in the binary PPV array {A} 
   with the brush (structuring element) {b}.
   
   Specifically, if {erode} is true, sets {A[ix]} to '0' if there is any
   voxel '0' in the neighborhood that is the brush translated by {ix}.
   This operation erodes the set of '1' voxels by expanding every '0'
   voxel into a set of '0' voxels whose indices are the central
   reflection of the indices in {b}.
   
   Conversely, {erode} is false, sets {A[ix]} to '1' if there is any voxel '1' in the
   neighbrhood.  This operation dilates the set of '1' voxels by expanding every '1'
   voxel into a set of '1' voxels with those same indices. 
   
   This procedure runs in time proportional to {N*} where {N} and {M}
   are the voxel counts of {A} and {b},respectively. */

#endif

