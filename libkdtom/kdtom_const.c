/* See {kdtom_const.h}. */
/* Last edited on 2021-07-19 05:05:44 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <ppv_array.h>
#include <jsmath.h>
#include <affirm.h>
#include <ixbox.h>

#include <kdtom.h>

#include <kdtom_const.h>

#define ixFMT ppv_index_t_FMT
#define smpFMT ppv_sample_t_FMT
#define szFMT ppv_size_t_FMT

ppv_sample_t kdtom_const_get_core_sample(kdtom_const_t *T, ppv_index_t dx[])
  {
    assert(T->h.kind == kdtom_kind_CONST);
    ppv_dim_t d = T->h.d;
    /* Check indices, just to be chato: */
    for (ppv_axis_t k = 0; k < d; k++) { assert((0 <= dx[k]) && (dx[k] < T->h.size[k])); }
    return T->smp;
  }

kdtom_const_t *kdtom_const_make
  ( ppv_dim_t d,
    ppv_sample_t maxsmp, 
    ppv_sample_t fill,
    ppv_index_t ixlo[],
    ppv_size_t size[], 
    ppv_sample_t smp
  )
  { 
    demand((0 <= d) && (d <= ppv_MAX_DIM), "invalid num of axes"); 
    demand((0 <= maxsmp) && (maxsmp <= ppv_MAX_SAMPLE_VAL), "invalid {maxsmp}"); 
    demand(fill <= maxsmp, "invalid {fill} value"); 
    demand(smp <= maxsmp, "invalid core sample value"); 
    
    /* Allocate the full record including all internal vectors: */
    size_t tot_bytes = kdtom_const_node_bytesize(d); /* Total node bytesize including all vectors. */
    kdtom_const_t *T = (kdtom_const_t *)notnull(malloc(tot_bytes), "no mem");
    
    /* Set {pend} to the free space inside the record, past all fixed fields: */
    size_t fix_bytes = sizeof(kdtom_const_t);  /* Size of fixed felds incl. those of {T.h}. */
    char *pend = addrsync(((char*)T) + fix_bytes, 8);
    
    /* Initialize the {T.h} fields: */
    kdtom_node_init((kdtom_t*)(&T->h), kdtom_kind_CONST, d, maxsmp, fill, ixlo, size, tot_bytes, &pend); 
    assert(T->h.kind == kdtom_kind_CONST);

    if (T->h.size[0] == 0) 
      { /* Core domain is empty, so {T->smp} is irrelevant: */
        T->smp = fill; 
      }
    else
      { T->smp = smp; }
    return T;
  }
  
kdtom_const_t *kdtom_const_clone(kdtom_const_t *T)
  {
    kdtom_const_t *S = kdtom_const_make(T->h.d, T->h.maxsmp, T->h.fill, T->h.ixlo, T->h.size, T->smp);
    return S;
  }

kdtom_const_t *kdtom_const_clip_core(kdtom_const_t *T, ppv_index_t ixlo[], ppv_size_t size[])
  {
    ppv_dim_t d = T->h.d;
    demand(! ixbox_is_empty(d, size), "new core domain is empty");
    demand(ixbox_is_contained(d, ixlo, size, T->h.ixlo, T->h.size), "new core domain not contained");
    kdtom_const_t *S = kdtom_const_make(d, T->h.maxsmp, T->h.fill, ixlo, size, T->smp);
    return S;
  }

size_t kdtom_const_node_bytesize(ppv_dim_t d)
  { 
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid dimension {d}");
    
    size_t fixf_bytes = sizeof(kdtom_const_t);   /* Fixed fields incl those of {T.h}. */
    size_t tot_bytes = iroundup(fixf_bytes, 8);  /* Account for address sync. */

    size_t ixlo_bytes = d * sizeof(ppv_size_t);  /* Bytesize for {T.h.ixlo} vector. */
    tot_bytes += iroundup(ixlo_bytes, 8);        /* Paranoia, account for address sync. */

    size_t sizv_bytes = d * sizeof(ppv_size_t);  /* Bytesize for {T.h.size} vector. */
    tot_bytes += iroundup(sizv_bytes, 8);        /* Paranoia, account for address sync. */

    return tot_bytes;
  }

size_t kdtom_const_bytesize(kdtom_const_t *T)
  {
    size_t tot_bytes = kdtom_const_node_bytesize(T->h.d);
    return tot_bytes;
  }

bool_t kdtom_const_is_all_fill(kdtom_const_t *T, ppv_sample_t fill)
  {
    if (T->h.fill != fill) { return FALSE; }
    if (T->h.size[0] == 0) { /* Empty core: */ return TRUE; }
    if (T->smp != fill) { /* Non-trivial core, not {fill}: */ return FALSE; }
    return TRUE;
  }

kdtom_const_t *kdtom_const_join_nodes
  ( ppv_size_t size[],
    ppv_axis_t ax,
    kdtom_const_t *T0, 
    ppv_size_t size0, 
    kdtom_const_t *T1, 
    ppv_size_t size1
  )
  {
    /* Basic compatibility: */
    if (T0->h.d != T1->h.d) { return NULL; }
    ppv_dim_t d = T0->h.d;

    demand(ax < d, "invalid {ax}");
    demand(size0 + size1 == size[ax], "inconsistent {size0,size1}");

    if (T0->h.maxsmp != T1->h.maxsmp) { return NULL; }
    ppv_sample_t maxsmp = T0->h.maxsmp;

    if (T0->h.fill != T1->h.fill) { return NULL; }
    ppv_sample_t fill = T0->h.fill;
    assert(fill <= maxsmp);
    
    /* Tests for trivially empty cores: */
    if ((T0->h.size[0] == 0) || (T0->smp == T0->h.fill)) 
      { /* {T0} is all fill: */ 
        /* Translate {T} by {size0} along axis {ax}: */
        kdtom_translate_one((kdtom_t *)T1, ax, size0);
        return T1;
      }
    if ((T1->h.size[0] == 0) || (T1->smp == T1->h.fill)) 
      { /* {T1} is all fill: */
        return T0;
      }
    if ((T0->h.ixlo[ax] >= size0) || (T0->h.ixlo[ax] + T0->h.size[ax] <= 0))
      { /* Core of {T0} is outside clip area: */ 
        return T1;
      }
    if ((T1->h.ixlo[ax] >= size1) || (T1->h.ixlo[ax] + T1->h.size[ax] <= 0))
      { /* Core of {T1} is outside clip area: */
        return T0;
      }
    
    /* There will be pieces of both cores, so: */
    if (T0->smp != T1->smp) { /* Distinct core values: */ return NULL; }
    ppv_sample_t smp = T0->smp;
    
    /* Check whether clipped cores can join into a single box: */
    if (T0->h.ixlo[ax] + T0->h.size[ax] < size0)
      { /* Core of {T0} does not reach to the joint: */ return NULL; }
    if (T1->h.ixlo[ax] > 0)
      { /* Core of {T1} starts after the joint: */ return NULL; }
    for (ppv_axis_t k = 0; k < d; k++)
      { if (k != ax)
          { /* Check if core has same extent along axis {k}: */
            if (T0->h.ixlo[k] != T1->h.ixlo[k]) { return NULL; }
            if (T0->h.size[k] != T1->h.size[k]) { return NULL; }
          }
      }
    
    /* No further objections, your honor: */
    kdtom_const_t *T = kdtom_const_make(d, maxsmp, fill, NULL, size, smp);
    return T;
  }

void kdtom_const_print_fields(FILE *wr, kdtom_const_t *T)
  {
    fprintf(wr, " .smp = " smpFMT, T->smp);

    return;
  }
