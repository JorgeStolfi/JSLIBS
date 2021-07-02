/* An internal k-d-tree node is basically a {ppv_array_t}. */
/* Last edited on 2021-07-01 15:49:31 by jstolfi */

#ifndef kdtom_array_H
#define kdtom_array_H
#define _GNU_SOURCE
#include <stdio.h>

#include <kdtom.h>
#include <indexing.h>
#include <ppv_array.h>
#include <bool.h>

typedef struct kdtom_array_t 
  { kdtom_t h;            /* General fields. */
    ppv_step_t *step;        /* Addressing increments. */
    ppv_pos_t base;          /* Base position. */
    ppv_nbits_t bpw;         /* Bits per word. */
    void *el;                /* Address of storage area. */
  } kdtom_array_t;
  /* A record {T} of this type describes a block of voxels as a
    {ppv_array_t} object {A}.  The array descriptor consists of the 
    {d,bps,size} fields of the head {h} plus the {step,base,bpw,el} fields
    specific of this node type.
    
    The {h} field must be the first field in the
    record, and {h.kind} must be {kdtom_kind_ARRAY}.
    
    The field {kdtom_kind_t} must be the first field in the record. */

kdtom_array_t *kdtom_array_make(ppv_array_t *A);
  /* Creates a {kdtom_array_t} node {T} from the given array. The
    storage area {T->el} will be the same as {A-el}, so neither of {A}
    or {T} should be completely reclaimed wile the other is still in
    use. */

ppv_sample_t kdtom_array_get_sample(kdtom_array_t *T, ppv_index_t ix[]);
  /* Obtains the sample {T.v[ix]} if {ix} is a valid index vector for {T}, otherwise
    bombs out. */

size_t kdtom_array_node_bytesize(ppv_dim_t d);
  /* Size in bytes of a {kdtom_array_t} node {T} with dimension {d}, including the {T.h.size}
    and {T.step} vectors, but NOT including the sample storage area {*T->el}. */

size_t kdtom_array_bytesize(kdtom_array_t *T, bool_t total);
  /* If {total} is false, returns the size in bytes used by the array
    record {*T}. The result includes the storage used by the {T.h}
    header and its {T.h.size} vector, as well as the {T->step}
    vector; but NOT the sample storage area {T->el}. The {reptoo} parameter
    is ignored.
    
    If {total} is true, the result also includes the nominal size in
    bytes of the sample storage area. This number is derived from the
    total count of distinct sample positions in the array, not counting
    replicated samples, as in {ppv_sample_count} with {reptoo=FALSE};
    assuming that samples are packed as tightly as allowed by the
    packing parameters {T.h.bps} and {T.bpw}.  The actual storage 
    area {T->el} may be larger than this amount. */

#endif
