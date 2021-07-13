/* An internal k-d-tree node that has a constant value over the whole domain. */
/* Last edited on 2021-07-13 01:14:11 by jstolfi */

#ifndef kdtom_const_H
#define kdtom_const_H

#define _GNU_SOURCE
#include <stdio.h>

#include <kdtom.h>
#include <ppv_array.h>
#include <bool.h>

typedef struct kdtom_const_t 
  { kdtom_t h;          /* General node parameters. */
    ppv_sample_t smp;   /* Sample value. */
  } kdtom_const_t;
  /* A record {T} of this type describes an infinite grid {T.V} of voxels whose
    values are all equal to {T.smp} within the core region {T.K}, and to {T.h.fill}
    outside it.  
    
    The core value {T.smp} is equal to {T.fill} if and only if the core is empty.
    
    The header field {h} must be the first field in the record. */

ppv_sample_t kdtom_const_get_core_sample(kdtom_const_t *T, ppv_index_t dx[]);
  /* Returns the sample {T.V[T.ixlo + dx]}. IMPORTANT: Assumes that this
    index is insde the core domain {T.DK}; that is, {dx} is in
    {0..T.size[k]-1}. The outcome is undefined otherwise. */

kdtom_const_t *kdtom_const_make
  ( ppv_dim_t d,
    ppv_sample_t maxsmp, 
    ppv_sample_t fill,
    ppv_index_t ixlo[],
    ppv_size_t size[], 
    ppv_sample_t smp
  );
  /* Creates a constant node {T} with the given parameters {T.h.d,
    T.h.maxsmp, T.h.fill, T.smp}. The number of axes {d} should not be zero.
    
    If {smp==fill} or {size == NULL} or {size[k] == 0} for any axis {k},
    the core will be empty ({T.size[k]==0} for all {k}) otherwise
    {T.size} will be copied from {size}. If {ixlo} is {NULL}, {T.ixlo} will be set
    to all zeros. */

kdtom_const_t *kdtom_const_clone(kdtom_const_t *T);
  /* Returns a copy of the node {T}. Copies all the fixed fields
    as well as the internally allocated vectors {T.ixlo,T.size}. */

kdtom_const_t *kdtom_const_clip(kdtom_const_t *T, ppv_index_t ixlo[], ppv_size_t size[]);
  /* Returns a {kdtom_t} structure {S} that describes the same voxel
    grid as {T}, except that its core {S.K} is the box {B} with
    low corner {ixlo[0..T.d-1]} and size {size[0..T.d-1]}. 
    
    IMPORTANT: The procedures assumes that {B} is not empty and is
    contained in {T.DK}.
    
    The parameters {S.maxsmp} and {S.fill} will be copied from {T}.
    
    The result {S} will be a newly allocated node, even it describes the
    same grid as {T}. It may not have the same kind as {T}. Note that
    {T} and some of its descendants may not be reachable from {S}. */

bool_t kdtom_const_is_all_fill(kdtom_const_t *T, ppv_sample_t fill);
  /* True iff {T.V} is everywhere equal to {fill}. */

kdtom_const_t *kdtom_const_join_nodes
  ( ppv_size_t size[],
    ppv_axis_t ax,
    kdtom_const_t *T0, 
    ppv_size_t sz0, 
    kdtom_const_t *T1, 
    ppv_size_t sz1
  );
  /* If nodes {T0} and {T1} can be joined in a single {kdtom_const_t}
    node, creates and returns that node. Otherwise returns {NULL}.
    
    Assumes that {T0} and {T1} are clipped so that the indice in their
    cores, along the axis {ax}, are contained in {0..sz0-1} and
    {0..sz1-1}, respectively; and in {0..size[k]-1} along any other axis {k}.
    
    The joined node {T}, if any, will have all zeros in {T.ixlo}, and
    will be such that {T.V[ix]} will be {T0.V[ix]} if {ix[ax] < sz0},
    and {T1.V[jx]} otherwise; were {jx[k]} is {ix[k]} for all {k},
    except that {jx[ax] = ix[ax]-sz0}. */

size_t kdtom_const_node_bytesize(ppv_dim_t d);
  /* Size in bytes of a {kdtom_const_t} node {T}, including the {T.h.size}
    vector. The number of axes {d} should not be zero. */

size_t kdtom_const_bytesize(kdtom_const_t *T);
  /* Returns the size in bytes used by the record {*T}. The result
    includes the storage used by the {T.h} header and its
    {T.h.size} vector, as well as any internal type-specific
    fields. */

#endif
