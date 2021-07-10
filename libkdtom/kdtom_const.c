/* See {kdtom_const.h}. */
/* Last edited on 2021-07-08 15:46:56 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <ppv_array.h>
#include <jsmath.h>
#include <affirm.h>

#include <kdtom.h>
#include <kdtom_const.h>

kdtom_const_t *kdtom_const_alloc(ppv_dim_t d);
  /* Allocats a {kdtom_const_t} record {T}, including the internal {T.h.size} vector.
    Initializes only the {T.h.d} and {T.h.kind} fields. */

kdtom_const_t *kdtom_const_make(ppv_dim_t d, ppv_nbits_t bps, ppv_size_t size[], ppv_sample_t smp)
  { 
    demand((0 <= d) && (d <= ppv_MAX_DIM), "invalid num of axes"); 
    demand((0 <= bps) && (bps <= ppv_MAX_BPS), "invalid bits per sample"); 
    ppv_sample_t maxsmp = ppv_max_sample(bps);
    demand((0 <= smp) && (smp <= maxsmp), "invalid sample value"); 
    
    kdtom_const_t *T = kdtom_const_alloc(d);

    T->h.bps = bps; 
    T->smp = smp; 
    for (ppv_axis_t k = 0; k < d; k++)
      { ppv_size_t szk = size[k];
        demand((0 <= szk) && (szk <= ppv_MAX_SIZE), "invalid size");
        T->h.size[k] = szk;
      }
    return T;
  }
   
size_t kdtom_const_node_bytesize(ppv_dim_t d)
  { 
    size_t fixf_bytes = sizeof(kdtom_const_t);   /* Fixed fields incl those of {T.h}. */
    size_t tot_bytes = iroundup(fixf_bytes, 8);  /* Account for address sync. */

    size_t sizv_bytes = d * sizeof(ppv_size_t);  /* Bytesize for {T.h.size} vector. */
    tot_bytes += iroundup(sizv_bytes, 8);        /* Paranoia, account for address sync. */

    return tot_bytes;
  }

kdtom_const_t *kdtom_const_alloc(ppv_dim_t d)
  { 
    /* Allocate the full record including all internal vectors: */
    size_t tot_bytes = kdtom_const_node_bytesize(d); /* Total node bytesize including all vectors. */
    kdtom_const_t *T = (kdtom_const_t *)notnull(malloc(tot_bytes), "no mem");
    
    /* Set {pend} to the free space inside the record, past all fixed fields: */
    size_t fix_bytes = sizeof(kdtom_const_t);  /* Size of fixed felds incl. those of {T.h}. */
    char *pend = addrsync(((char*)T) + fix_bytes, 8);
    
    /* Initialize the {T.h} fields, including the {T.h.size} vector: */
    kdtom_node_init((kdtom_t *)T, tot_bytes, d, &pend);
    T->h.kind = kdtom_kind_CONST;

    return T;
  }

ppv_sample_t kdtom_const_get_sample(kdtom_const_t *T, ppv_index_t ix[])
  {
    assert(T->h.kind == kdtom_kind_CONST);
    ppv_dim_t d = T->h.d;
    /* Check indices, just to be chato: */
    for (ppv_axis_t k = 0; k < d; k++)
      { ppv_size_t szk = T->h.size[k];
        demand((0 <=ix[k]) && (ix[k] < szk), "invalid index");
      }
    return T->smp;
  }

size_t kdtom_const_bytesize(kdtom_const_t *T)
  {
    size_t tot_bytes = kdtom_const_node_bytesize(T->h.d);
    return tot_bytes;
  }
