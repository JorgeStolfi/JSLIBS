/* Multidimensional sample arrays stored as k-d-trees. */
/* Last edited on 2021-06-25 01:31:33 by jstolfi */

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
    in some operations. */
    
ppv_sample_t kdtom_get_sample(kdtom_t  *T, ppv_index_t ix[]);
  /* Obtains the sample {T.v[ix]} if {ix} is a valid index vector, otherwise
    bombs out. */
    
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

kdtom_t *kdtom_alloc(ppv_dim_t d, size_t rec_bytes, char **pendP);
  /* Allocates an area {T} on the heap for a {kdtom_t} record, 
    and sets the size vector address {T->size} to point to 
    a place within the same record.  The {T->d field is set to {d},
    but all other fields are left undefined.
    
    The parameter {rec_bytes} must be at least {hdr_bytes =
    sizof(kdtom_t)}. The record will actually have bytesize {rec_bytes +
    szv_bytes}, where {szv_bytes = d*sizeof(ppv_size_t)} is the
    bytesize used by by the {T.size} vector. The latter will be located
    right after the first {hdr_bytes} bytes. 
    
    The address of the allocated record is returned as the result. if
    {pendP} is not {NULL}, and {rec_bytes} is strictly greater than
    {hdr_bytes}, the address of the "free" area (if any) in the record,
    right after the {size} vector, is reurned in {*pendP}.
    
    Thus the command {free(T)} will reclaim the space used by the
    {kdtom_t} record AND by the vector {T->size}, plus any other fields
    of the record that are specific to the variant. One should NEVER
    call {free(T->size)}. */
    
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


#endif
