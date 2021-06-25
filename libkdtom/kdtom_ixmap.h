/* An inrernal k-d-tree node that permutes, reverses, and shifts indices. */
/* Last edited on 2021-06-24 01:30:30 by jstolfi */

#ifndef kdtom_ixmap_H
#define kdtom_ixmap_H

#define _GNU_SOURCE
#include <stdio.h>

#include <kdtom.h>
#include <ppv_array.h>
#include <bool.h>

typedef struct kdtom_ixmap_t 
  { kdtom_t head;                /* Common {kdtom_t} fields. */
    ppv_axis_t *prax;            /* Index remapping vector. */
    bool_t *rvix;                /* Index reversing flags. */
    ppv_index_t *ixlo;           /* Index shift. */
    ppv_sample_t fill;           /* Surround fill value. */
    kdtom_t *sub;                /* Sub-block. */
  } kdtom_ixmap_t;
  /* A descriptor for a block {vt} of voxels in the {d}-dimensional
    infinite voxel grid, represetnted as another k-d-tree node {sub}
    with some remapping of the indices. 
    ??? Complete ??? */
 
kdtom_ixmap_t *kdtom_ixmap_make
  ( ppv_dim_t d,
    ppv_size_t size[],
    ppv_axis_t prax[],
    bool_t rvix[],
    ppv_index_t ixlo[],
    ppv_sample_t fill, 
    kdtom_t *sub
  );
  /* Creates an index-mapping node {T} with sub-node {sub}, ???. */
   
ppv_sample_t kdtom_ixmap_get_sample(kdtom_ixmap_t *T, ppv_index_t ix[]);
  /* Obtains the sample {T.v[ix]} if {ix} is a valid index vector for
    {T}, otherwise bombs out. */

#endif
