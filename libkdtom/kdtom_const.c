/* See {kdtom_const.h}. */
/* Last edited on 2021-06-25 06:50:16 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <ppv_array.h>
#include <affirm.h>

#include <kdtom.h>
#include <kdtom_const.h>

kdtom_const_t *kdtom_const_alloc(ppv_dim_t d);
  /* Allocats a {kdtom_const_t} record {T}, including the internal {T.head.size} vector.
    Initializes only the {T.head.d} and {T.head.kind} fields. */

kdtom_const_t *kdtom_const_make(ppv_dim_t d, ppv_nbits_t bps, ppv_size_t size[], ppv_sample_t val)
  { 
    demand((0 <= d) && (d <= ppv_MAX_DIM), "invalid num of axes"); 
    demand((0 <= bps) && (bps <= ppv_MAX_BPS), "invalid bits per sample"); 
    ppv_sample_t maxval = (ppv_sample_t)((1<<bps)- 1);
    demand((0 <= val) && (val <= maxval), "invalid sample value"); 
    
    kdtom_const_t *T = kdtom_const_alloc(d);

    T->head.bps = bps; 
    T->val = val; 
    for (ppv_axis_t k = 0; k < d; k++)
      { ppv_size_t szk = size[k];
        demand((0 <= szk) && (szk <= ppv_MAX_SIZE), "invalid size");
        T->head.size[k] = szk;
      }
    return T;
  }
   
size_t kdtom_const_node_size(ppv_dim_t d)
  { 
    size_t fix_bytes = sizeof(kdtom_const_t);   /* Bytesize of all fixed fields incl. {head} */
    size_t size_bytes = d * sizeof(ppv_size_t); /* Bytesize of the {head.size} vector. */
    size_t tot_bytes = fix_bytes + size_bytes;
    return tot_bytes;
  }

kdtom_const_t *kdtom_const_alloc(ppv_dim_t d)
  { 
    size_t fix_bytes = sizeof(kdtom_const_t); /* Bytesize of all fixed fields incl. {head} */
    kdtom_const_t *T = (kdtom_const_t *)kdtom_alloc(d, fix_bytes, NULL);
    T->head.kind = kdtom_kind_CONST;
    assert(T->head.d == d); 
    return T;
  }

ppv_sample_t kdtom_const_get_sample(kdtom_const_t *T, ppv_index_t ix[])
  {
    assert(T->head.kind == kdtom_kind_CONST);
    ppv_dim_t d = T->head.d;
    /* Check indices, just to be chato: */
    for (ppv_axis_t k = 0; k < d; k++)
      { ppv_size_t szk = T->head.size[k];
        demand((0 <=ix[k]) && (ix[k] < szk), "invalid index");
      }
    return T->val;
  }
