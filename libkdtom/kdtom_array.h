/* An internal k-d-tree node is basically a {ppv_array_t}. */
/* Last edited on 2021-06-25 06:37:35 by jstolfi */

#ifndef kdtom_array_H
#define kdtom_array_H
#define _GNU_SOURCE
#include <stdio.h>

#include <kdtom.h>
#include <indexing.h>
#include <ppv_array.h>
#include <bool.h>

typedef struct kdtom_array_t 
  { kdtom_t head;            /* General fields. */
    ppv_step_t *step;        /* Addressing increments. */
    ppv_pos_t base;          /* Base position. */
    ppv_nbits_t bpw;         /* Bits per word. */
    void *el;                /* Address of storage area. */
  } kdtom_array_t;
  /* A record {T} of this type describes a block of voxels as a
    {ppv_array_t} object {A}.  The array descriptor consists of the 
    {d,bps,size} fields of {head} plus the {step,base,bpw,el} fields
    of the record.
    
    The {head} must be the first field in the
    record, and {head.kind} must be {kdtom_kind_ARRAY}.
    
    The field {kdtom_kind_t} must be the first field in the record. */

kdtom_array_t *kdtom_array_make(ppv_array_t *A);
  /* Creates a {kdtom_array_t} node {T} from the given array. The
    storage area {T->el} will be the same as {A-el}, so neither of {A}
    or {T} should be completely reclaimed wile the other is still in
    use. */

ppv_sample_t kdtom_array_get_sample(kdtom_array_t *T, ppv_index_t ix[]);
  /* Obtains the sample {T.v[ix]} if {ix} is a valid index vector for {T}, otherwise
    bombs out. */

size_t kdtom_array_node_size(ppv_dim_t d);
  /* Size in bytes of a {kdtom_array_t} node {T}, including the {T.head.size}
    and {T.step} vectors, but NOT including the sample storage area {*T->el}. */

#endif
