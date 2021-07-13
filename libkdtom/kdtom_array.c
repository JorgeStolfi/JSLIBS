/* See {kdtom_array.h}. */
/* Last edited on 2021-07-12 11:21:27 by jstolfi */

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
#include <kdtom_array.h>

kdtom_array_t *kdtom_array_make(ppv_array_t *A, ppv_sample_t fill)
  { 
    ppv_dim_t d = A->d;
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid array dimension {d}");
    ppv_sample_t maxsmp = A->maxsmp;
    demand(fill <= maxsmp, "invalid {fill} value");
    
    /* Allocate the full record including all internal vectors: */
    size_t tot_bytes = kdtom_array_node_bytesize(d); /* Total node bytesize including all vectors. */
    kdtom_array_t *T = (kdtom_array_t *)notnull(malloc(tot_bytes), "no mem");
    
    /* Set {pend} to the free space inside the record, past all fixed fields: */
    size_t fix_bytes = sizeof(kdtom_array_t);  /* Size of fixed felds incl. those of {h}. */
    char *pend = addrsync(((char*)T) + fix_bytes, 8);
    
    /* Initialize the {h} fields, including the {T.h.size} vector: */
    kdtom_node_init((kdtom_t *)T, kdtom_kind_ARRAY, d, maxsmp, fill, NULL, A->size, tot_bytes, &pend);
    
    /* Allocate and fill the {T.step} vector: */
    size_t elsz = sizeof(ppv_step_t);
    pend = addrsync(pend, 8);
    T->step = (ppv_step_t *)kdtom_alloc_internal_vector((kdtom_t *)T, tot_bytes, d, elsz, &pend);
    for (ppv_axis_t ax = 0; ax < d; ax++) { T->step[ax] = A->step[ax]; }
    
    /* Initialize all remaining fields: */
    T->base = A->base;
    T->bps = A->bps;
    T->bpw = A->bpw;
    T->el = A->el;
    return T;
  }
  
kdtom_t *kdtom_array_clip(kdtom_array_t *T, ppv_index_t ixlo[], ppv_size_t size[])
  {
    /* !!! Implement !!! */
    fprintf(stderr, "!! {kdtom_array_clip} not implemented\n");
    return NULL;
  }

size_t kdtom_array_node_bytesize(ppv_dim_t d)
  {
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid dimension {d}");

    size_t fixf_bytes = sizeof(kdtom_array_t); /* Fixed fields incl those of head part {h}. */
    size_t tot_bytes = iroundup(fixf_bytes, 8); /* Account for address sync. */

    size_t ixlo_bytes = d * sizeof(ppv_size_t);  /* Bytesize for {T.h.ixlo} vector. */
    tot_bytes += iroundup(ixlo_bytes, 8);        /* Paranoia, account for address sync. */

    size_t sizv_bytes = d * sizeof(ppv_size_t); /* Bytesize for {h.size} vector. */
    tot_bytes += iroundup(sizv_bytes, 8);       /* Paranoia, account for address sync. */

    size_t step_bytes = d * sizeof(ppv_step_t); /* Bytesize of the {step} vector. */
    tot_bytes += iroundup(step_bytes, 8);       /* Paranoia, account for address sync. */

    return tot_bytes;
  }

ppv_sample_t kdtom_array_get_core_sample(kdtom_array_t *T, ppv_index_t dx[])
  {
    assert(T->h.kind == kdtom_kind_ARRAY);
    ppv_pos_t pos = ix_position(T->h.d, dx, T->base, T->step);
    ppv_sample_t smp = ppv_get_sample_at_pos(T->el, T->bps, T->bpw, pos);
    return smp;
  }

size_t kdtom_array_bytesize(kdtom_array_t *T, bool_t total)
  {
    size_t node_bytes = kdtom_array_node_bytesize(T->h.d);
    size_t tot_bytes = node_bytes;
    if (total)
      { /* Create a temporary array descriptor: */
        ppv_array_t A = (ppv_array_t)
          { .d = T->h.d,
            .size = T->h.size,
            .step = T->step,
            .base = 0,
            .bps = T->bps,
            .bpw = T->bpw,
            .el = T->el
          };
        /* Count samples ignoring replicated ones: */
        ppv_sample_count_t npos = ppv_sample_count(&A, FALSE);
        /* Compute nominal bytes needed to store them: */
        size_t samp_bytes = ppv_tot_sample_bytes(npos, T->bps, T->bpw);
        tot_bytes += samp_bytes;
      }
    return tot_bytes;
  }
