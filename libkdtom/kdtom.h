/* Multidimensional sample arrays stored as k-d-trees. */
/* Last edited on 2021-06-28 00:35:51 by jstolfi */

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
    kdtom_kind_IXMAP   /* The voxels indices are remapped. */
  } kdtom_kind_t;
  /* A code that identifies the kind of a {kdtom_t}. */

typedef struct kdtom_t 
  { kdtom_kind_t kind;  /* Kind of node. */
    ppv_dim_t d;        /* Dimensinality (number of axes). */
    ppv_nbits_t bps;    /* Bits per sample. */
    ppv_size_t *size;   /* Count of voxels along each axis. */
  } kdtom_t;
  /* A {kdtom_t} {T} describes a finite multi-dimensional array of
    samples {T.v}, internally stored as a k-d-tree.
    
    Compared to a simple {ppv_array_t}, a {kdtom_t} is much more
    economical of space if the array has large regions with uniform
    sample values and relatively simple boundaries. Unlike a
    {ppv_array_t}, it is read-only -- the voxel values cannot be
    changed.
    
    The array has {T.d} axes (dimensions). Each sample is a
    {ppv_sample_t} value (an unsigned inteegr) in the range
    {0..2^bps-1}. If {T.bps} is zero, all samples have value zero.
    
    Each element of {T.v} (a /voxel/) is identified by an index vector
    {ix[0..d-1]}. The voxel exists iff {ix[k]} is in {0..T.size[k]-1}
    for each axis {k}; and undefined otherwise.
    
    The type {kdtom_t} is a union type; that is, a {kdtom_t} record is
    actually one of several record types that are distinguished by the
    field {kind}.  The record will actually be a {kdtom_const_t} record if
    {kind} is {kdtom_kind_CONST}, a {kdtom_array_t} record if {kind} is
    {kdtom_kind_ARRAY}, and so on. Thus a {kdtom_t} record should never
    be allocated; one of those specific record types should be allocated
    instead, and its address can then be cast as a {kdtom_t*}.
    
    A {kdtom_t} may have zero dimensions ({T.d = 0}), in which case it
    has a single voxel whose index vector is the empty vector. On the
    other hand, if {T.D > 0} but {T.size[ax]} is zero for some axis
    {ax}, the {kdtom_t} is /empty/ -- there are no voxels. The other
    elements of {T.size} are still significant and may make a difference
    in some operations.
    
    When any kind of {kdom_t} node record {T} is allocated on the heap,
    the {T.size} vector is usually allocated within the same heap memory
    area. Then the command {free(T)} will reclaim the space used by the
    {kdtom_t} record AND by the vector {T->size}, plus any other fields
    of the record that are specific to the variant. In that case one
    should NOT call {free(T->size)}, or serious tire damage may result. */
    
/* SAMPLE EXTRACTION */

ppv_sample_t kdtom_get_sample(kdtom_t  *T, ppv_index_t ix[]);
  /* Obtains the sample {T.v[ix]} if {ix} is a valid index vector, otherwise
    bombs out. */
    
/* RECURSIVE TREE BUILDNG AND MANIPULATION */

kdtom_t *kdtom_grind_array(ppv_array_t *A);
  /* Creates a {kdtom_t} tree from the array {A} by recursively splitting
    and analyzing each part.
    
    The internal nodes of the tree will be {kdtom_split_t} nodes. The
    leaves will be either {kdtom_const_t} or {kdtom_array_t} nodes.
    
    Any substantial part of {A} found in the search whose volxels have
    all the same value is represented by a const node.
    
    Otherwise the part is represented by an array node, with a newly
    allocated sample storage area, and the relevant voxels of {A} are
    copied ino that area. Thus {A} can be wholly reclaimed after the
    procedure returns.
    
    The recursive subdivision will proceed until the part of {A} under
    analysis would take more space if represented as a k-d-tree than if
    it was kept as a single array node. 
    
    If the two sides of a split are const nodes with the same  value,
    they are merged into one. Ditto if they are both array nodes.*/

kdtom_t kdtom_clip(kdtom_t *T, ppv_axis_t ax, ppv_index_t ixlo, ppv_index_t sz);
  /* Returns a {kdtom_t} structure {S} that describes the same voxel array as {T},
    except that along the axis {ax} it is clipped to the index range {ixlo..ixlo+sz-1}.
    This range must be contained in the range {0..T.size[ax]-1}.
    
    That is, voxel {S.v[ix]} is equal to {T.v[jx]} where {jx[k]=ix[k]}
    for each {k} in {0..T.d-1}, except that {jx[ax]=ix[ax]+ixlo};
    provided that {ix[ax]} is in {0..sz-1}, and {jx[k]} is in
    {0..T.size[k]} fr all {k}.
    
    In any case (even if {S} is empty) {S.size[k]} will be {sz} if
    {k=ax}, and {T.size[k]} otherwise.*/
    
/* FOR USE BY VARIANT IMPLEMENTATIONS */

void kdtom_node_init(kdtom_t *T, size_t sztot, ppv_dim_t d, char **pendP);
  /* Assumes that {T} is pointing to a memory area with {sztot} bytes,
    which is large enough to contain a bare {d}-dimensional {kdtom_t}
    node record including the {size} vector. Also assumes that {*pendP)
    is the address within that memory area where the {size} vector can
    be stored.
    
    The procedure sets the size vector address {T->size} to {*pendP} and
    updates {*pendP} to point just past the end of that vector. It also
    sets the {T->d} field to {d}, but leaves all other fields undefined.
    It fails if there is not enough space from{*pendP} to the end of the 
    record.  */
    
char *kdtom_alloc_internal_vector(kdtom_t *T, size_t sztot, ppv_dim_t d, size_t elsz, char **pendP);
  /* Assumes that {T} is pointing to an area with {sztot} bytes and
    {*pendP} is the address within that area. Returns {*pendP} after
    making sure that there is enough space from there to the end of the
    record to hold a vector of {d} elements each with {elsz} bytes. */

size_t kdtom_bytesize(kdtom_t *T, bool_t total);
  /* If {total} is false, returns the total size in bytes used by the
    record {*T}, taking into account its kind. The result includes the
    storare used by the {T->size} vector, as well as any internal
    type-specific fields. In case of an array node, however, it does NOT
    include its voxels storage area.

    If {total} is true, returns instead the total size in bytes used by
    all the k-d-tree nodes reached from {T}, includin the node records
    themselves and the nominal size of the storage area of every array
    node (as returned by {kdtom_array_total_bytesize}). */

#endif
