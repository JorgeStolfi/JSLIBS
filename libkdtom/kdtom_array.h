/* An internal k-d-tree node is basically a {ppv_array_t}. */
/* Last edited on 2024-12-05 10:32:39 by stolfi */

#ifndef kdtom_array_H
#define kdtom_array_H

#include <stdio.h>

#include <kdtom.h>
#include <ix.h>
#include <ppv_array.h>
#include <bool.h>

typedef struct kdtom_array_t 
  { kdtom_t h;               /* Header with general fields. */
    ppv_step_t *step;        /* Addressing increments of array. */
    ppv_pos_t base;          /* Base position (of element with zero index vector). */
    ppv_nbits_t bps;         /* Bits per sample of the array. */
    ppv_nbits_t bpw;         /* Bits per word of the array storage. */
    void *el;                /* Address of storage area. */
  } kdtom_array_t;
  /* A record {T} of this type describes an infinite grid of voxels
    {T.V} whose core part {T[k]} is the contents of a {ppv_array_t}
    object {A}. The voxel {T.V[T.ixlo]} will
    the voxel of {A} with zero index.
    
    The array descriptor consists of the {.d,.maxsmp,.size} fields of
    the head {h} plus the {.step,.base,.bps,.bpw,.el} fields specific of
    this node type.
    
    The {h} field must be the first field in the record, and {h.kind}
    must be {kdtom_kind_ARRAY}.
    
    The field {kdtom_kind_t} must be the first field in the record. */

kdtom_array_t *kdtom_array_make(ppv_array_t *A, ppv_sample_t fill);
  /* Creates a {kdtom_array_t} node {T} whose core is the given array.
  
    The fields {T.h.d,T.h.maxsmp,T.h.size,T.step,T.bps,T.bpw,T.el} will
    be copied from {A}. The number of axes {A.d} should not be zero. The
    vector {T.ixlo} will be set to all zeros, and {T.fill} will be the
    specified {fill} value.
    
    Note that, since the storage area {*(A.el)} is shared, neither {A} nor
    {T} should be completely reclaimed wile the other is still in
    use. */

kdtom_array_t *kdtom_array_do_make
  ( ppv_dim_t d,
    ppv_sample_t  maxsmp,
    ppv_sample_t fill,
    ppv_index_t ixlo[],
    ppv_size_t size[],
    ppv_pos_t base,
    ppv_step_t step[],
    ppv_nbits_t bps,
    ppv_nbits_t bpw,
    void *el
  );
  /* Allocates a new {kdtom_array_t} node and fills it with the given
    fields. The contents of {ixlo,size,step} are copied into
    {T.h.ixlo,T.h.size,T.step}, which are allocated in the record
    {*T} itself. */

kdtom_array_t *kdtom_array_clone(kdtom_array_t *T);
  /* Returns a copy of the node {T}. Copies the fixed fields
    and the internally allocated vectors {T.h.ixlo,T.h.size,T.step},
    but does NOT copy the voxel storage area {*(T.el)}. */

ppv_sample_t kdtom_array_get_core_sample(kdtom_array_t *T, ppv_index_t dx[]);
  /* Returns the sample {T.V[T.ixlo + dx]}. IMPORTANT: Assumes that this
    index is insde the core domain {T.DK}; that is, {dx} is in
    {0..T.size[k]-1}. The outcome is undefined otherwise. */

kdtom_t *kdtom_array_clip_core(kdtom_array_t *T, ppv_index_t ixlo[], ppv_size_t size[]);
  /* Returns a {kdtom_t} structure {S} that describes the same voxel
    grid as {T}, except that its core {S.K} is the box {B} with
    low corner {ixlo[0..T.d-1]} and size {size[0..T.d-1]}.
    
    IMPORTANT: The procedures assumes that {B} is not empty and is
    contained in {T.DK}.

    The parameters {S.maxsmp} and {S.fill} will be copied from {T}.
    
    The result {S} will be a newly allocated {kdtom_array_t} node, even
    it describes the same grid as {T}. The node {T} itself is not
    changed. Note that {T} and some of its descendants may not be
    reachable from {S}. */

size_t kdtom_array_node_bytesize(ppv_dim_t d);
  /* Size in bytes of a {kdtom_array_t} node {T} with dimension {d}, including the {T.h.size}
    and {T.step} vectors, but NOT including the sample storage area {*(T.el)}.
    The number of axes {d} should not be zero. */

size_t kdtom_array_bytesize(kdtom_array_t *T, bool_t total);
  /* If {total} is false, returns the size in bytes used by the array
    record {*T}. The result includes the storage used by the {T.h}
    header and the vectors {T.h.ixlo}, {T.h.size}, and {T.step};
    but NOT the sample storage area {*(T.el)}.
    
    If {total} is true, the result also includes the NOMINAL size in
    bytes of the sample storage area. This number is derived from the
    total count of distinct sample positions in the array, not counting
    replicated samples, as in {ppv_sample_count} with {reptoo=FALSE};
    assuming that samples are packed as tightly as allowed by the
    packing parameters {T.bps} and {T.bpw}.  
    
    Note that the actual storage area {*(T.el)} may be larger than this
    amount, e.g. if the elements of {A} are NOT packed as tightly as
    possible. */

ppv_array_t *kdtom_array_make_descr(kdtom_array_t *T);
  /* Allocates and returns a {ppv_array_t} descriptor {A}
    that describes the array of voxels that is the core {T.K} of {T}.
    
    The storage area {*(A.el)} will be shared with {*(T.el)}. The
    vectors {A.size} and {A.step} will be allocated in the same reacord
    {*A} an their contents will be copied from {T.h.size} and
    {T.step}. */

/* DEBUGGING */

void kdtom_array_print_fields(FILE *wr, kdtom_array_t *T);
  /* Prints to {wr} the fields that are specific of {T}, 
    as ".{field} = {value}". */

#endif
