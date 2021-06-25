/* See {kdtom_ixmap.h}. */
/* Last edited on 2021-06-24 01:33:50 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <ppv_array.h>
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
   
kdtom_ixmap_t *kdtom_ixmap_alloc(ppv_dim_t d)
  { 
    size_t fix_bytes = sizeof(kdtom_ixmap_t); /* Bytesize of all fixed fields incl. {head} */
    kdtom_ixmap_t *T = (kdtom_ixmap_t *)kdtom_alloc(d, fix_bytes, NULL);
    T->head.kind = kdtom_kind_IXMAP;
    assert(T->head.d == d); 
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
