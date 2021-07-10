/* An inrernal k-d-tree node that shifts and pads a sub-tree. */
/* Last edited on 2021-07-08 15:51:22 by jstolfi */

#ifndef kdtom_trans_H
#define kdtom_trans_H

#define _GNU_SOURCE
#include <stdio.h>

#include <kdtom.h>
#include <ppv_array.h>
#include <bool.h>

typedef struct kdtom_trans_t 
  { kdtom_t h;                /* Common {kdtom_t} fields. */
    ppv_index_t *ixlo;           /* Index shift. */
    ppv_sample_t fill;           /* Surround fill value. */
    kdtom_t *sub;                /* Sub-tree. */
  } kdtom_trans_t;
  /* A descriptor {T} for a block {T.V} of voxels in the {d}-dimensional
    infinite voxel grid, represetnted as another k-d-tree {sub}
    but shifted and padded.
    
    Namely, a sample T.V[ix]} of {T}, as usual, has index vector {ix[0..d-1]} where,
    for each {k}, {ix[k]} is in the range {0..tsz[k]-1]} and {tsz = T.h.size}.
    The sample is undefined for other index vectors.
    
    If, for every {k}, {ix[k]} is in the range {lo[k]..lo[k]+usz[k]-1}}, where {lo = T.ixlo]}
    and {usz = T.sub.h.size}, then that sample is the sample of {T.sub}
    with index vector {jx}, where {jx[k] = ix[k] - lo[k]}.  Otherwise, that sample
    has the value {T.fill}. 
    
    Note that any parts of the domain of {T.sub} that, shfted by {ixlo},
    fall outside the domain of {T}, are inacessible through {T}, and thus
    implicitly clipped to {T}'s domain. */
 
kdtom_trans_t *kdtom_trans_make
  ( ppv_dim_t d,
    ppv_size_t size[],
    ppv_index_t ixlo[],
    ppv_sample_t fill, 
    kdtom_t *sub
  );
  /* Creates a shifting and padding node {T} with {d} axes and
    dimensions specified by {size[0..d-1]}, with the parameters {fill} and
    {ixlo[0..d-1]}, and the sub-tree {sub}. */
   
ppv_sample_t kdtom_trans_get_sample(kdtom_trans_t *T, ppv_index_t ix[]);
  /* Obtains the sample {T.V[ix]} if {ix} is a valid index vector for
    {T}, otherwise bombs out. */

size_t kdtom_trans_node_bytesize(ppv_dim_t d);
  /* Size in bytes of a {kdtom_trans_t} node {T}, including the {T.sub}
    pointer and the vectors {T.h.size} and {T.ixlo} but
    NOT the nodes in the sub-tree pointed to {T.sub}. */

size_t kdtom_trans_bytesize(kdtom_trans_t *T, bool_t total);
  /* If {total} is false, returns the size in bytes used by the
    node record {*T}. The result includes the storage used by the {T.h}
    header and the vectors {T.h.size,T.ixlo}
    and the {T.sub} pointer; but NOT the other nodes accessed 
    through that pointer.
    
    If {total} is true, the result also includes the total size in
    bytes of all nodes that can be reached through {T.sub},
    as returned by {kdtom_bytesize(T.sub,TRUE)}. .*/

#endif
