/* An inrernal k-d-tree node that divides the domain betwee two nodes. */
/* Last edited on 2021-06-25 06:37:47 by jstolfi */

#ifndef kdtom_split_H
#define kdtom_split_H

#define _GNU_SOURCE
#include <stdio.h>

#include <kdtom.h>
#include <ppv_array.h>
#include <bool.h>

typedef struct kdtom_split_t 
  { kdtom_t head;       /* General node parameters. */
    ppv_axis_t ax;           /* Coordinate axis perpendicular to cut. */
    kdtom_t *sub[2];    /* Representation of the two half-blocks. */
  } kdtom_split_t;
  /* A record {T} of this type describes a block of voxels as the
    justaposition of two blocks {T.sub[0]} and {T.sub[1]} that share a
    face perpendicular to axis {T.ax}.
    
    The two blocks must have {sz[k]} voxels along each axis {k}, where {sz}
    is {T.head.size}; except when {k} is {T.ax}.  Along that axis,
    block {T.sub[0]} must have a certain size {sz0} in {0..sz[ax]}, and {T.sub[1]} 
    must have size {T.head.sz[ax]-sz0}.
    
    The voxel {T.v[ix]} is {T.sub[0].v[ix]} if {ix[ax]} is less than 
    {sz0 = sub0.size[ax]}, and {T.sub[1].v[jx]} otherwise, where
    {jx] is {ix} with {sz0} subtracted from {ix[ax]}.
    
    The field {head} must be the first field in the record. */

kdtom_split_t *kdtom_split_make(ppv_axis_t ax, kdtom_t *sub0, kdtom_t *sub1);
  /* Creates a split node {T} with sub-nodes {sub0} and {sub1}, the
    latter being displaced along axis {ax} so that it shares with {sub0}
    a face perpendicular to axis {ax}. Their sizes along each axis {k}
    must be equal, and will be the sizes of {T}; except along axis {ax}.
    where the sizes can be arbitrary, and their sum will be the size of
    {T}. */

ppv_sample_t kdtom_split_get_sample(kdtom_split_t *T, ppv_index_t ix[]);
  /* Obtains the sample {T.v[ix]} if {ix} is a valid index vector for
    {T}, otherwise bombs out. The index vector {ix} is be modified by
    the procedure, but is restored before returning. */

size_t kdtom_split_node_size(ppv_dim_t d);
  /* Size in bytes of a {kdtom_split_t} node {T}, including the {T.head.size}
    vector and the {T.sub[0..1]} pointers, but NOT the nodes in those
    sub-trees. */

#endif
