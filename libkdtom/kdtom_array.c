/* See {kdtom_array.h}. */
/* Last edited on 2021-06-25 06:49:01 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <ppv_array.h>
#include <affirm.h>

#include <kdtom.h>
#include <kdtom_array.h>

kdtom_array_t *kdtom_array_alloc(ppv_dim_t d);
  /* Allocats a {kdtom_array_t} record {T}, including the internal
    {T.head.size} and {T.step} vectors. Initializes only the {T.head.d}
    and {T.head.kind} fields. */

kdtom_array_t *kdtom_array_make(ppv_array_t *A)
  { 
    ppv_dim_t d = A->d;
    kdtom_array_t *T = kdtom_array_alloc(A->d);
    for (ppv_axis_t ax = 0; ax < d; ax++)
      { T->head.size[ax] = A->size[ax];
        T->step[ax] = A->step[ax];
      }
    
    T->head.bps = A->bps;
    T->base = A->base;
    T->bpw = A->bpw;
    T->el = A->el;
    return T;
  }
  
size_t kdtom_array_node_size(ppv_dim_t d)
  {
    size_t fix_bytes = sizeof(kdtom_array_t);   /* Bytesize of all fixed fields incl. {head} */
    size_t size_bytes = d * sizeof(ppv_size_t); /* Bytesize of the {head.size} vector. */
    size_t step_bytes = d * sizeof(ppv_step_t); /* Bytesize of the {step} vector. */
    size_t tot_bytes = fix_bytes + size_bytes + step_bytes;
    return tot_bytes;
  }

kdtom_array_t *kdtom_array_alloc(ppv_dim_t d)
  {
    size_t fix_bytes = sizeof(kdtom_array_t);  /* Size of fixed felds incl. {head}. */
    size_t stp_bytes = d * sizeof(ppv_step_t); /* Size of {T.step} vector.*/
    size_t rec_bytes = fix_bytes + stp_bytes; /* Record size minus {size} vector. */
    char *pend;
    kdtom_array_t *T = (kdtom_array_t *)kdtom_alloc(d, rec_bytes, &pend);
    assert(T->head.d == d);
    T->head.kind = kdtom_kind_ARRAY;
    T->step = (ppv_step_t *)pend; pend += stp_bytes;
    return T;
  }

ppv_sample_t kdtom_array_get_sample(kdtom_array_t *T, ppv_index_t ix[])
  {
    assert(T->head.kind == kdtom_kind_ARRAY);
    ppv_pos_t pos = ix_position(T->head.d, ix, T->base, T->step);
    ppv_sample_t v = ppv_get_sample_at_pos(T->el, T->head.bps, T->bpw, pos);
    return v;
  }
