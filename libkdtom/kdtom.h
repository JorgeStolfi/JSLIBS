/* Multidimensional sample arrays stored as k-d-trees. */
/* Last edited on 2021-07-12 22:01:27 by jstolfi */

#ifndef kdtom_H
#define kdtom_H

#define _GNU_SOURCE
#include <stdio.h>

#include <ppv_array.h>
#include <bool.h>

typedef enum
  { kdtom_kind_CONST,  /* All voxels are constant. */
    kdtom_kind_ARRAY,  /* The voxels are represented by an array of samples. */
    kdtom_kind_SPLIT,  /* The voxels are split in twosub-blocks. */
  } kdtom_kind_t;
  /* A code that identifies the kind of a {kdtom_t}. */

#define kdtom_kind_FIRST kdtom_kind_CONST
#define kdtom_kind_LAST  kdtom_kind_SPLIT

typedef struct kdtom_t 
  { kdtom_kind_t kind;    /* Kind of node. */
    ppv_dim_t d;          /* Dimensinality (number of axes), positive. */
    ppv_sample_t maxsmp;  /* Max sample value. */
    ppv_sample_t fill;    /* Surround sample value. */
    ppv_index_t *ixlo;    /* Low index of core. */
    ppv_size_t *size;     /* Count of voxels along each axis. */
  } kdtom_t;
  /* A {kdtom_t} {T} describes an infinite multi-dimensional grid of
    /voxels/ {T.V}, each holding a /sample/ (unsigned integer) value.
    The grid consists of a finite /core/ grid {T.K}, described by a
    k-d-tree structure, surrounded by {T.fill} samples on all sides.
    
    Compared to a simple {ppv_array_t}, a {kdtom_t} is much more
    economical of space if the array has large regions with uniform
    sample values and relatively simple boundaries. Unlike a
    {ppv_array_t}, it is read-only -- the voxel values cannot be
    changed.
    
    Axes, indices, and sample values
    
    The voxel grid {T.V} has a positive number {T.d} axes (dimensions).
    Each voxel {T.v[ix]} is identified by an /index vector/ {ix[0..d-1]}
    of signed integers.
    
    Each sample is a {ppv_sample_t} value (an unsigned integer) in the
    range {0..T.maxsmp}. In particular,if {T.maxsmp} is zero, all samples have
    value zero.
    
    Core and fill
    
    The core {T.K} has {T.size[k]} voxels along each axis {k}, and its
    lowest index along that axis is {T.ixlo[k]}. Thus, {T.V[ix]} is in
    the core if {ix[k]-ixlo[k]} is in {0..T.size[k]-1} for each axis
    {k}.  The set of index vectors of the core voxels is the /core
    domain/, denoted {T.DK}.
    
    If {ix} is in the core domain {T.DK}, the value of {T.V[ix]} is
    obtained from the k-d-tree structure. Otherwise the value of
    {T.V[ix]} is the fixed sample {T.fill}, which must be in
    {0..T.maxsmp}.
    
    If {T.size[ax]} is zero for some axis {ax}, the core
    grid {T.K} is /empty/ -- has no voxels. In that case, all elements
    of {T.ixlo} and {T.size} will be zero, and all samples {T.V} are
    equal to {T.fill}.
    
    Allocation
    
    The type {kdtom_t} is a union type; that is, a {kdtom_t} record is
    actually one of several record types that are distinguished by the
    field {kind}.  The record will actually be a {kdtom_const_t} record if
    {kind} is {kdtom_kind_CONST}, a {kdtom_array_t} record if {kind} is
    {kdtom_kind_ARRAY}, and so on. Thus a {kdtom_t} record should never
    be allocated; one of those specific record types should be allocated
    instead, and its address can then be cast as a {kdtom_t*}.
    
    When any kind of {kdom_t} node record {T} is allocated on the heap,
    the {T.size} and {T.ixlo} vectors are usually allocated within the
    same heap memory area. Then the command {free(T)} will reclaim the
    space used by the {kdtom_t} record AND by those two vectors, plus
    any other fields of the record that are specific to the variant. In
    that case one should NOT call {free(T->size)} or {free(T->ixlo)}, or
    severe tire damage may result.
    
    Modifying the attributes
    
    The fields {T.d} and {T.kind} should never be modified after the
    node is created. The fields {T.size} and {T.maxsmp} can be modified
    only in some cases and for some kinds of nodes, and may require
    changing other nodes in the whole tree.

    The field {T.fill} can be changed to any value not exceeding {T.maxsmp}.
    
    The fields {T.ixlo} of the root of a tree can be changed at will,
    having the effect of translating the whole grid {T.V}.. */
    
bool_t kdtom_is_all_fill(kdtom_t *T, ppv_sample_t fill);
  /* True iff {T} is a {kdtom_const_t} node and {T.V} is everywhere equal to {fill}. */

bool_t kdtom_has_empty_core(kdtom_t  *T);
  /* Returns {TRUE} if and only if {T} has empty core. */

kdtom_t *kdtom_join_nodes
  ( ppv_size_t size[], 
    ppv_axis_t ax,
    kdtom_t *T0, 
    ppv_size_t sz0, 
    kdtom_t *T1, 
    ppv_size_t sz1
  );
  /* If nodes {T0} and {T1} are {kdtom_const_t} nodes that can be joined
    in a single {kdtom_const_t} node, creates and returns that node.
    Otherwise returns {NULL}.
    
    Assumes that {T0} and {T1} are clipped so that the indice in their
    cores, along the axis {ax}, are contained in {0..sz0-1} and
    {0..sz1-1}, respectively; and in {0..size[k]-1} along any other axis {k}.
    
    The joined node {T}, if any, will have all zeros in {T.ixlo}, and
    will be such that {T.V[ix]} will be {T0.V[ix]} if {ix[ax]<sz0}, and
    {T1.V[jx]} otherwise; were {jx[k]} is {ix[k]} for all {k}, except
    that {jx[ax] = ix[ax]-sz0}. */

ppv_sample_t kdtom_get_sample(kdtom_t  *T, ppv_index_t ix[]);
  /* Obtains the sample {T.V[ix]}; from the code {T.K} if {ix} is in the core
  domain {T.DK}, otherwise returns the {T.fill} value. */
    
/* RECURSIVE TREE BUILDING AND MANIPULATION */

void kdtom_translate(kdtom_t  *T, ppv_index_t dx[]);
  /* Translates the whole grid {T.V} by the vector {dx[0..T.d-1]}, by
    adding {dx} to the {T.ixlo} vector; unless the core domain {T.DK} is
    empty, in which case the operation is a no-op. */

void kdtom_realloc_array_nodes(kdtom_t *T);
  /* Replaces the storage area of every {kdtom_array_t} node {S} in {T}
    by an independent newly allocated array {B} of voxels, with the same
    parameters {S.d, S.maxsmp, S.size},and copies the samples of {S.K}
    into it. The parameters {B.bps}, {B.bpw}, and {B.step} are
    recomputed, based on {S.maxsmp} and {S.size}, so as to pack the
    elements as tightly as possible.
    
    This procedure can be called after {kdtom_grind_array(A,fill)} to
    allow that {A.el} can be safely reclaimed. */

kdtom_t *kdtom_clip(kdtom_t *T, ppv_index_t ixlo[], ppv_size_t size[]);
  /* Returns a {kdtom_t} structure {S} that describes the same voxel
    grid as {T}, except that its core {S.K} is the core {T.K} clipped
    to the box {B} with low corner {ixlo[0..T.d-1]} and size {size[0..T.d-1]}.
    If {ixlo} is NULL, it is assumed to be the all-zeros vector. 
    Note that {B} may be empty.
    
    The parameters {S.maxsmp} and {S.fill} will be copied from {T}.
    
    The result {S} will be a newly allocated node, even it describes the
    same grid as {T}. It may not have the same kind as {T}. Note that
    {T} and some of its descendants may not be reachable from {S}. */
    
bool_t kdtom_box_is_empty(ppv_dim_t d, ppv_size_t  size[]);
  /* Returns {TRUE} if any one of the components {size[0..d-1]}
    is zero. */

void kdtom_intersect_boxes
  ( ppv_dim_t d,
    ppv_index_t ixlo_A[], 
    ppv_size_t  size_A[],
    ppv_index_t ixlo_B[], 
    ppv_size_t  size_B[],
    ppv_index_t ixlo_R[], 
    ppv_size_t  size_R[]
  );
  /* Takes two {d}-dimensional boxes -- {A}, defined by the low-corner {ixlo_A[0..d-1]}
    and the sizes  {size_A[0..d-1]}, and {B}, defined similarly by {ixlo_B} and {size_B};
    and computes their intersection {R}, defined by {ixlo_R} and {size_R}. 
    
    If the intersection is empy (in particular, if any {size_A[k]} or {size_B[k]} is empty),
    sets {ixlo_R} and {size_R} to alll zeros.*/
    
/* FOR USE BY VARIANT IMPLEMENTATIONS */

void kdtom_node_init
  ( kdtom_t *T, 
    kdtom_kind_t kind,
    ppv_dim_t d, 
    ppv_sample_t maxsmp,
    ppv_sample_t fill,
    ppv_index_t ixlo[], 
    ppv_size_t size[], 
    size_t sztot,
    char **pendP
  );
  /* Assumes that {T} is pointing to a memory area with {sztot} bytes,
    which is large enough to contain a bare {d}-dimensional {kdtom_t}
    node record including the {ixlo} and {size} vectors. Also assumes
    that {*pendP) is the address within that memory area where the those
    vectors can be stored.
    
    The procedure sets the size vector addresses {T.ixlo} and {T.size}
    with {kdtom_alloc_internal_vector} to point into the free area
    specified by {*pendP}, updating {*pendP} to point to the unused
    space still available in that area. It fails if there is not enough
    space from{*pendP} to the end of the record, as specified by
    {sztot}.
    
    The procedure also sets the fields {T->kind, T->d, T->maxsmp, T->fill} to the
    given values, and copies the vectors {ixlo[0..d-1]} and
    {size[0..d-1]} into {T.ixlo} and {T.size}. However, if any element
    {size[k]} is zero, it sets both vectors to all zeros.
    
    All other fields of {T} are left undefined. */
    
char *kdtom_alloc_internal_vector(kdtom_t *T, size_t sztot, ppv_dim_t d, size_t elsz, char **pendP);
  /* Assumes that {T} is pointing to an area with {sztot} bytes and
    {*pendP} is the address within that area. Returns {*pendP} after
    making sure that there is enough space from there to the end of the
    record to hold a vector of {d} elements each with {elsz} bytes. 
    Also updates {*pendP} to point to the end of that vector. */

size_t kdtom_bytesize(kdtom_t *T, bool_t total);
  /* If {total} is false, returns the size in bytes used by the
    record {*T}, taking into account its kind. The result includes the
    storare used by the {T.ixlo} and {T.size} vectors, as well as any internal
    type-specific fields. In case of an array node, however, it does NOT
    include its voxels storage area.

    If {total} is true, returns instead the total size in bytes used by
    all the k-d-tree nodes reached from {T}, includin the node records
    themselves and the nominal size of the storage area of every array
    node (as returned by {kdtom_array_total_bytesize}). */

#endif
