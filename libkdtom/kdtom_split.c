/* See {kdtom_split.h}. */
/* Last edited on 2021-06-24 01:37:15 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <ppv_array.h>
#include <affirm.h>

#include <kdtom.h>
#include <kdtom_split.h>

kdtom_split_t *kdtom_split_alloc(ppv_dim_t d);
  /* Allocats a {kdtom_split_t} record {T}, including the internal {T.head.size} vector.
    Initializes only the {T.head.d} and {T.head.kind} fields. */

kdtom_split_t *kdtom_split_make(ppv_axis_t ax, kdtom_t *sub0, kdtom_t *sub1)
  { 
    ppv_dim_t d = sub0->d; 
    
    kdtom_split_t *T = kdtom_split_alloc(d);
    
    demand(sub1->d == T->head.d, "{d} mismatch");
    T->head.bps = sub0->bps; 
    demand(sub1->bps = T->head.bps, "{bps} mismatch");
    for (ppv_axis_t k = 0; k < d; k++)
      { if (k != ax)
          { ppv_size_t szk = sub0->size[k];
            T->head.size[k] = szk;
            demand(sub1->size[k] == szk, "size mismatch");
          }
        else
          { ppv_size_t szk = sub0->size[k] + sub1->size[k];
            T->head.size[k] = szk;
          }
      }
    T->ax = ax;
    T->sub[0] = sub0;
    T->sub[1] = sub1;
    return T;
  }
   
kdtom_split_t *kdtom_split_alloc(ppv_dim_t d)
  { 
    size_t fix_bytes = sizeof(kdtom_split_t); /* Bytesize of all fixed fields incl. {head} */
    kdtom_split_t *T = (kdtom_split_t *)kdtom_alloc(d, fix_bytes, NULL);
    T->head.kind = kdtom_kind_SPLIT;
    assert(T->head.d == d); 
    return T;
  }

ppv_sample_t kdtom_split_get_sample(kdtom_split_t *T, ppv_index_t ix[])
  {
    assert(T->head.kind == kdtom_kind_SPLIT);
    ppv_dim_t d = T->head.d;
    ppv_axis_t ax = T->ax; assert((0 <= ax) && (ax < d));
    
    /* Save original index of axis {ax}: */
    ppv_index_t ix_save = ix[ax];
    
    /* Decide where to go, and adjust index: */
    ppv_size_t sz0 = T->sub[0]->size[ax]; /* Size of low block on axis {ax}. */
    int32_t ksub;     /* Which subnode (0 or 1) has the desired voxel. */
    if (ix[ax] < sz0)
      { ksub = 0; }
    else
      { ksub = 1; ix[ax] = ix[ax] - sz0; }
    
    /* Get the sample: */
    ppv_sample_t v = kdtom_get_sample(T->sub[ksub], ix);
    
    /* Restore original index vector: */
    ix[ax] = ix_save;
    return v;
  }
