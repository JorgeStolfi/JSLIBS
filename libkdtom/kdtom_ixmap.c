/* See {kdtom_ixmap.h}. */
/* Last edited on 2021-06-28 03:31:25 by jstolfi */

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
#include <kdtom_ixmap.h>

kdtom_ixmap_t *kdtom_ixmap_alloc(ppv_dim_t d);
  /* Allocats a {kdtom_ixmap_t} record {T}, including the internal 
    vectors {T.head.size},  {T->pmix}, {T->rvix}, {T->ixlo}.
    Initializes only the {T.head.d} and {T.head.kind} fields. */

kdtom_ixmap_t *kdtom_ixmap_make
  ( ppv_dim_t d,
    ppv_size_t size[],
    ppv_axis_t prax[],
    bool_t rvix[],
    ppv_index_t ixlo[],
    ppv_sample_t fill, 
    kdtom_t *sub
  )
  { 
    ppv_dim_t d_sub = sub->d; 
    ppv_dim_t bps_sub = sub->bps; 
    
    kdtom_ixmap_t *T = kdtom_ixmap_alloc(d);
    T->head.bps = bps_sub;
    T->fill = fill;
    T->sub = sub;

    for (ppv_axis_t k = 0; k < d; k++)
      { ppv_size_t szk = size[k];
        demand((0 <= szk) && (szk <= ppv_MAX_SIZE), "invalid size");
        T->head.size[k] = szk;
        ppv_axis_t praxk = prax[k];
        demand((0 <= praxk) && (praxk < d_sub), "invalid index remap");
        T->prax[k] = praxk;
        ppv_index_t ixlok = ixlo[k];
        demand((-ppv_MAX_INDEX <= ixlok) && (ixlok <= ppv_MAX_INDEX), "invalid index shift");
        T->ixlo[k] = ixlok;
        T->rvix[k] = rvix[k];
      }
    return T;
  }
   
size_t kdtom_ixmap_node_bytesize(ppv_dim_t d)
  { 
    size_t fixf_bytes = sizeof(kdtom_ixmap_t);   /* Fixed fields incl those of {head}. */
    size_t tot_bytes = iroundup(fixf_bytes, 8);  /* Account for address sync. */

    size_t sizv_bytes = d * sizeof(ppv_size_t);  /* Bytesize for {head.size} vector. */
    tot_bytes += iroundup(sizv_bytes, 8);        /* Paranoia, account for address sync. */

    size_t ixlo_bytes = d * sizeof(ppv_index_t); /* Bytesize of the {ixlo} vector. */
    tot_bytes += iroundup(ixlo_bytes, 8);        /* Paranoia, account for address sync. */

    size_t prax_bytes = d * sizeof(ppv_axis_t);  /* Bytesize of the {prax} vector. */
    tot_bytes += iroundup(prax_bytes, 8);        /* Paranoia, account for address sync. */

    size_t rvix_bytes = d * sizeof(bool_t);      /* Bytesize of the {rvix} vector. */
    tot_bytes += iroundup(rvix_bytes, 8);        /* Paranoia, account for address sync. */

    return tot_bytes;
  }

kdtom_ixmap_t *kdtom_ixmap_alloc(ppv_dim_t d)
  { 
    /* Allocate the full record including all internal vectors: */
    size_t tot_bytes = kdtom_ixmap_node_bytesize(d); /* Total node bytesize including all vectors. */
    kdtom_ixmap_t *T = (kdtom_ixmap_t *)notnull(malloc(tot_bytes), "no mem");
    
    /* Set {pend} to the free space inside the record, past all fixed fields: */
    size_t fix_bytes = sizeof(kdtom_ixmap_t);  /* Size of fixed felds incl. those of {head}. */
    char *pend = addrsync(((char*)T) + fix_bytes, 8);

    /* Initialize the {head} fields, including the {T.head.size} vector: */
    kdtom_node_init((kdtom_t *)T, tot_bytes, d, &pend);
    assert(T->head.d == d);
    T->head.kind = kdtom_kind_IXMAP;

    /* Allocate the variant-specific vectors{T.prax,T.rvix,T.ixlo}: */
    
    /* Start with the {ixlo} vector to avoid multiple {addrsync}s: */
    size_t ixlo_elsz = sizeof(ppv_index_t);
    assert(ixlo_elsz == 8);
    assert(((uint64_t)pend) % 8 == 0); /* Because we just allocated the {T.head.size} vector. */
    T->ixlo = (ppv_index_t *)kdtom_alloc_internal_vector((kdtom_t *)T, tot_bytes, d, ixlo_elsz, &pend);
    
    size_t prax_elsz = sizeof(ppv_axis_t);
    assert(prax_elsz == 1);  /* So no {addrsync} is necessary. */
    T->prax = (ppv_axis_t *)kdtom_alloc_internal_vector((kdtom_t *)T, tot_bytes, d, prax_elsz, &pend);

    size_t rvix_elsz = sizeof(bool_t);
    assert(rvix_elsz == 1);  /* So no {addrsync} is necessary. */
    T->rvix = (bool_t *)kdtom_alloc_internal_vector((kdtom_t *)T, tot_bytes, d, rvix_elsz, &pend);

    /* Check for overflow. Note that some bytes may be left over due to rounding: */
    assert (pend - (char*)T <= tot_bytes);
    
    return T;
  }

ppv_sample_t kdtom_ixmap_get_sample(kdtom_ixmap_t *T, ppv_index_t ix[])
  {
    assert(T->head.kind == kdtom_kind_IXMAP);
    // ppv_dim_t d = T->head.d;
    // ppv_axis_t ax = T->ax; assert((0 <= ax) && (ax < d));
    // 
    // /* Save original index of axis {ax}: */
    // ppv_index_t ix_save = ix[ax];
    // 
    // /* Decide where to go, and adjust index: */
    // ppv_size_t sz0 = T->sub[0]->head.size[ax]; /* Size of low block on axis {ax}. */
    // int32_t ksub;     /* Which subnode (0 or 1) has the desired voxel. */
    // if (ix[ax] < sz0)
    //   { ksub = 0; }
    // else
    //   { ksub = 1; ix[ax] = ix[ax] - sz0; }
    // 
    // /* Get the sample: */
    ppv_sample_t v = kdtom_get_sample(T->sub, ix);
    //
    // /* Restore original index vector: */
    // ix[ax] = ix_save;
    return v;
  }

size_t kdtom_ixmap_bytesize(kdtom_ixmap_t *T, bool_t total)
  {
    size_t node_bytes = kdtom_ixmap_node_bytesize(T->head.d);
    size_t tot_bytes = node_bytes;
    if (total)
      { /* Add size of subtree: */
        size_t subt_bytes = kdtom_bytesize(T->sub, TRUE);
        tot_bytes += subt_bytes;
      }
    return tot_bytes;
  }
