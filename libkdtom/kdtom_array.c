/* See {kdtom_array.h}. */
/* Last edited on 2024-12-05 10:32:35 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <ppv_array.h>
#include <ix_io.h>
#include <jsmath.h>
#include <affirm.h>
#include <ixbox.h>

#include <kdtom.h>
#include <kdtom_const.h>

#include <kdtom_array.h>

#define ixFMT ppv_index_t_FMT
#define smpFMT ppv_sample_t_FMT
#define szFMT ppv_size_t_FMT
#define posFMT ppv_pos_t_FMT

kdtom_array_t *kdtom_array_make(ppv_array_t *A, ppv_sample_t fill)
  { 
    kdtom_array_t *T = kdtom_array_do_make
      ( A->d, A->maxsmp, fill, NULL, A->size, 
        A->base, A->step, A->bps, A->bpw, A->el
      );
    return T;
  }

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
  )
  { 
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid array dimension {d}");
    demand(fill <= maxsmp, "invalid {fill} value");
    
    /* Allocate the full record including all internal vectors: */
    size_t tot_bytes = kdtom_array_node_bytesize(d); /* Total node bytesize including all vectors. */
    kdtom_array_t *T = (kdtom_array_t *)notnull(malloc(tot_bytes), "no mem");
    
    /* Set {pend} to the free space inside the record, past all fixed fields: */
    size_t fix_bytes = sizeof(kdtom_array_t);  /* Size of fixed felds incl. those of {h}. */
    char *pend = addrsync(((char*)T) + fix_bytes, 8);
    
    /* Initialize the {h} fields, including the {T.h.size} vector: */
    kdtom_node_init((kdtom_t *)T, kdtom_kind_ARRAY, d, maxsmp, fill, ixlo, size, tot_bytes, &pend);
    
    /* Allocate and fill the {T.step} vector: */
    size_t elsz = sizeof(ppv_step_t);
    pend = addrsync(pend, 8);
    T->step = (ppv_step_t *)kdtom_alloc_internal_vector((kdtom_t *)T, tot_bytes, d, elsz, &pend);
    for (ppv_axis_t k = 0; k < d; k++) { T->step[k] = step[k]; }
    
    /* Initialize all remaining fields: */
    T->base = base;
    T->bps = bps;
    T->bpw = bpw;
    T->el = el;
    return T;
  }
  
kdtom_array_t *kdtom_array_clone(kdtom_array_t *T)
  {
    /* Create a temporary array descriptor: */
    kdtom_array_t *S = kdtom_array_do_make
       ( T->h.d, T->h.maxsmp, T->h.fill, T->h.ixlo, T->h.size, 
         T->base, T->step, T->bps, T->bpw, T->el
       );
    return S;
  }

kdtom_t *kdtom_array_clip_core(kdtom_array_t *T, ppv_index_t ixlo[], ppv_size_t size[])
  {
    ppv_dim_t d = T->h.d;
    
    demand(! ixbox_is_empty(d, size), "clipping box is empty");
    assert(size != NULL);

    ppv_array_t *A = kdtom_array_make_descr(T);
    for (ppv_axis_t k = 0; k < d; k++) 
      { /* Determine how many elements to skip along axis {k}: */
        ppv_index_t ixlok = (ppv_index_t)(ixlo == NULL ? 0 : ixlo[k]);
        assert(ixlok >= T->h.ixlo[k]);
        ppv_size_t skip = (ppv_size_t)(ixlok - T->h.ixlo[k]);
        
        /* Determine how many elements to keep along axis {k}: */
        ppv_size_t keep = size[k];
        assert(keep > 0);
        assert(skip + keep <= T->h.size[k]);
        
        /* Crop the array: */
        if ((skip > 0) || (keep < T->h.size[k])){ ppv_crop(A, k, skip, keep); }
        assert(A->size[k] == keep);
      }

    kdtom_t *S = (kdtom_t *)kdtom_array_make(A, T->h.fill);
    kdtom_translate(S, ixlo);
    free(A);
    return S;
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
      { /* Count samples ignoring replicated ones: */
        ppv_sample_count_t npos = ix_num_positions(T->h.d, T->h.size, T->step);
        /* Compute nominal bytes needed to store them: */
        size_t samp_bytes = ppv_tot_sample_bytes(npos, T->bps, T->bpw);
        tot_bytes += samp_bytes;
      }
    return tot_bytes;
  }
        
ppv_array_t *kdtom_array_make_descr(kdtom_array_t *T)
  { ppv_dim_t d = T->h.d;
    ppv_array_t *A = ppv_array_new_desc(d);
    A->maxsmp = T->h.maxsmp;
    for (ppv_axis_t k = 0; k < d; k++) 
      { A->size[k] = T->h.size[k];
        A->step[k] = T->step[k];
      }
    A->base = 0;
    A->bps = T->bps;
    A->bpw = T->bpw;
    A->el = T->el;
    return A;
  }

void kdtom_array_print_fields(FILE *wr, kdtom_array_t *T)
  {
    ix_print_steps(wr, ".step = [ ", T->h.d, T->step, 0, " ", " ]");
    fprintf(wr, " .base = " posFMT, T->base);
    fprintf(wr, " .bps = %d", T->bps);
    fprintf(wr, " .bpw = %d", T->bpw);
    fprintf(wr, " .el = %lu", ((uint64_t)T->el));

    return;
  }
