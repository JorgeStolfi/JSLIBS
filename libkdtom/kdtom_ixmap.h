/* An inrernal k-d-tree node that permutes, reverses, and shifts indices. */
/* Last edited on 2021-06-27 10:16:11 by jstolfi */

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

size_t kdtom_ixmap_node_bytesize(ppv_dim_t d);
  /* Size in bytes of a {kdtom_ixmap_t} node {T}, including the {T.sub}
    pointer and the vectors {T.head.size}, {T.prax}, and {T.ixlo} but
    NOT the nodes in the sub-tree pointed to {T.sub}. */

size_t kdtom_ixmap_bytesize(kdtom_ixmap_t *T, bool_t total);
  /* If {total} is false, returns the size in bytes used by the
    node record {*T}. The result includes the storage used by the {T.head}
    header and the vectors {T.head.size,T.prax,T.rvix,T.ixlo}
    and the {T.sub} pointer; but NOT the other nodes accessed 
    through that pointer.
    
    If {total} is true, the result also includes the total size in
    bytes of all nodes that can be reached through {T.sub},
    as returned by {kdtom_bytesize(T.sub,TRUE)}. .*/

#endif
