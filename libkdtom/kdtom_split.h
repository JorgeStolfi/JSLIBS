/* An inrernal k-d-tree node that divides the domain betwee two nodes. */
/* Last edited on 2021-07-11 18:47:05 by jstolfi */

#ifndef kdtom_split_H
#define kdtom_split_H

#define _GNU_SOURCE
#include <stdio.h>

#include <kdtom.h>
#include <ppv_array.h>
#include <bool.h>

typedef struct kdtom_split_t 
  { kdtom_t h;          /* General node parameters. */
    ppv_axis_t ax;      /* Coordinate axis perpendicular to cut. */
    ppv_index_t ixcut;  /* Index value of the cut. */
    kdtom_t *sub[2];    /* Representation of the two half-blocks. */
  } kdtom_split_t;
  /* A record {T} of this type describes an infinite grid {T.V} of
    voxels whose core consists of the the splice of two grids {T0.V} and
    {T1.V}, separated by a plane perpendicular to axis {ax} between the
    voxels with indices {ixcut-1} and {ixcut} on that axis. 
    
    The core {T.K} of {T} is defined by {T.h.ixlo} and {T.h.size}, as in
    any {kdtom_t} variant; and voxels outside that region are all equal
    to {T.fill}. Note that {T.fill}, {T.sub[0].fill}, and
    {T.sub[1].fill} may be all different.

    The axis counts {T.sub[0].d} and {T.sub[1].d} must be equal, and
    that will be the count {T.h.d}. Same for the {.maxsmp} fields.
    
    The axis {T.ax} must be in {0..d-1} where {d = T.h.d}, and the value
    of {T.ixcut} must be such that {T.ixcut - T.d.ixlo[ax]} is in
    {1..T.d.size[ax]-1}. Note that the core {T.K} must not be empty, and
    {t.h.size[ax]} must be at least 2.
    
    The core domain {T.DK} is split into two non-empty sub-domains:
    {T.DK0}, with with indices up to {T.ixcut-1} along axis {ax}, and
    {T.DK1} with indices from {T.ixcut} and higher. The value {T.K[ix]}
    of a voxel in {T}'s core is {T0.V[ix]} or {T1.V[ix]}, respectively.

    Any element {T0.V[ix]} with {ix} outside {T.DK0} is inaccessible
    through {T}. Thus normally the core domain {T0.DK} of {T0} is
    contained in the sub-domain {T.DK0}. Conversely, if {T0.DK} contains
    the whole of {T.DK0}, the fill value of {T0} is irrelevant. The same
    is true of {T1.V[ix]} relative to {T.DK1}.
     
    The grid {T0.V} is the grid {T.sub[0].V} implicitly translated so
    that its indices are relative to the low corner {T.ix0} of {T.DK0}.
    Similarly {T1.V} is {T.sub[1].V} implicitly translated by the low
    corner {T.ix1} of {T.DK1}. That is, a voxel {T.K[ix]} inside the
    core domain {T.DK} is either {T.sub[0].V[ix - T.ix0]} or
    {T.sub[1].V[ix - T.ix1]}, depending on whether {ix} is in {T.DK0} or
    {T.DK1}. The vector {T.ix0} is just {T.ixlo}, and {T.ix1} is the
    same except that {T.ix1[ax]} is T.ixcut}.
    
    The {.d} and {.maxval} fields of {T.sub[0]} and {T.sub[1]} must be equal to 
    those of {T.d}.  The {.fill} values may be all different.
   
    The header field {h} must be the first field in the record. */

kdtom_t *kdtom_split_make
  ( ppv_sample_t fill,
    ppv_index_t ixlo[], 
    ppv_size_t size[], 
    ppv_axis_t ax, 
    ppv_index_t ixcut,
    kdtom_t *sub0,
    kdtom_t *sub1
  );
  /* Tries to create a {kdtom_split_t} node {T} with sub-nodes {T.sub[0]
    = sub0} and {T.sub[1] = sub1}, with the specified values of
    {T.h.ixlo}, {T.h.size}, {T.h.fill}, {T.ax}, and {T.ixcut}. However,
    for certain values of the parameters, may return a simpler {kdtom_t}
    node that is equivalent to that.
    
    The normal case is when no element of {size} is zero, {size[ax]} is
    at least 2, and {ixcut-ixlo[ax]} is in {1..size[ax]-1}. Then in
    principle {T} will be a {kdtom_split_t} node as described above. In
    particular, the grids {T0.V} and {T1.V} will be {sub0.V} and
    {sub1.V} implicitly translated by the respective low index vectors
    {T.ix0} and {T.ix1} of the sub-domains {T.DK0} and {T.DK1}, and then
    clipped to those sub-domains.
    
    Special cases
    
    If any {size[k]} is zero, the core {T.K} will be empty, and {T.size}
    will be set to all zeros. Then the result {T} will be a
    {kdtom_const_t} node with the given {fill} value.
    
    If {ixcut-ixlo[ax]} is not positive, so that {T.DK0} is empty, then
    {T.K} will be just {T1.V} clipped to {T.DK}. Similarly, if
    {ixcut-ixlo[ax]} is {size[ax]} or more, so that {T.DK1} is empty,
    {T.K} will be just {T0.V} clipped to {T.DK}. In both cases the value
    outside {T.DK} will be the given {fill} value.
    
    More generally, if the grid {T.V} can be described by a single
    {kdtom_t} node with savings of space, it the replacement will take
    place. */

kdtom_split_t *kdtom_split_do_make
  ( ppv_sample_t fill,
    ppv_size_t size[], 
    ppv_axis_t ax, 
    kdtom_t *T0,
    ppv_size_t sz0,
    kdtom_t *T1,
    ppv_size_t sz1
  );
  /* Similar to {kdtom_split_make}, but assumes that the two sub-nodes cannot be simplified 
    or joned into a single node.  Therefore, always returns a {kdtom_split_t} node.
    
    Assumes that {T0} and {T1} are the sub-nodes, already clipped to the corresponding
    sub-domains ({T.DK0} and {T.DK1}) and each translated so that the corner of its
    sub-domain is at the origin.  Assumes that,along axis {ax}, the indices
    in the core {T0.DK} are contained in {0..sz0-1}, and those of{T1.DK}
    are contained in  0..sz1-1}. 
    
    The split node {T} created will have {T.h.ixlo} all zeros
    and {T.ixcut = sz0}. */

ppv_sample_t kdtom_split_get_core_sample(kdtom_split_t *T, ppv_index_t dx[]);
  /* Returns the sample {T.V[T.ixlo + dx]}. IMPORTANT: Assumes that this
    index is insde the core domain {T.DK}; that is, {dx} is in
    {0..T.size[k]-1}. The outcome is undefined otherwise. */

kdtom_t *kdtom_split_clip(kdtom_split_t *T, ppv_index_t ixlo[], ppv_size_t size[]);
  /* Returns a {kdtom_t} structure {S} that describes the same voxel
    grid as {T}, except that its core {S.K} is the box {B} with
    low corner {ixlo[0..T.d-1]} and size {size[0..T.d-1]}.
    
    IMPORTANT: The procedures assumes that {B} is not empty and is
    contained in {T.DK}.
    
    The parameters {S.maxsmp} and {S.fill} will be copied from {T}.
    
    The result {S} will be a newly allocated node, even it describes the
    same grid as {T}. It may not have the same kind as {T}. Note that
    {T} and some of its descendants may not be reachable from {S}. */

size_t kdtom_split_node_bytesize(ppv_dim_t d);
  /* Size in bytes of a {kdtom_split_t} node {T}, including the {T.h.ixlo} and {T.h.size}
    vector and the {T.sub[0..1]} pointers, but NOT the nodes in those
    sub-trees. */

size_t kdtom_split_bytesize(kdtom_split_t *T, bool_t total);
  /* If {total} is false, returns the size in bytes used by the node
    record {*T}. The result includes the storage used by the {T.h}
    header and the {T.h.ixlo} and {T.h.size} vectors, and the pointers
    {T.sub[0..1]}; but NOT the other nodes accessed through those two
    pointers.
    
    If {total} is true, the result also includes the total size in
    bytes of all nodes that can be reached through {T.sub[0..1]},
    as returned by {kdtom_bytesize(T.sub[i],TRUE)}. .*/

#endif
