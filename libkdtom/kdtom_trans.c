/* See {kdtom_trans.h}. */
/* Last edited on 2021-07-08 15:51:13 by jstolfi */

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
#include <kdtom_trans.h>

kdtom_trans_t *kdtom_trans_alloc(ppv_dim_t d);
  /* Allocats a {kdtom_trans_t} record {T}, including the internal 
    vectors {T.h.size} and {T->ixlo}.
    Initializes only the {T.h.d} and {T.h.kind} fields. */

kdtom_trans_t *kdtom_trans_make
  ( ppv_dim_t d,
    ppv_size_t size[],
    ppv_index_t ixlo[],
    ppv_sample_t fill, 
    kdtom_t *sub
  )
  { 
    demand(sub->d == d, "axis counts mismatch");
    ppv_dim_t bps = sub->bps; 
    
    kdtom_trans_t *T = kdtom_trans_alloc(d);
    T->h.bps = bps;
    T->fill = fill;
    T->sub = sub;

    for (ppv_axis_t k = 0; k < d; k++)
      { ppv_size_t szk = size[k];
        demand((0 <= szk) && (szk <= ppv_MAX_SIZE), "invalid size");
        T->h.size[k] = szk;
        ppv_index_t ixlok = ixlo[k];
        demand((-ppv_MAX_INDEX <= ixlok) && (ixlok <= ppv_MAX_INDEX), "invalid index shift");
        T->ixlo[k] = ixlok;
      }
    return T;
  }
   
size_t kdtom_trans_node_bytesize(ppv_dim_t d)
  { 
    size_t fixf_bytes = sizeof(kdtom_trans_t);   /* Fixed fields incl those of {T.h}. */
    size_t tot_bytes = iroundup(fixf_bytes, 8);  /* Account for address sync. */

    size_t sizv_bytes = d * sizeof(ppv_size_t);  /* Bytesize for {T.h.size} vector. */
    tot_bytes += iroundup(sizv_bytes, 8);        /* Paranoia, account for address sync. */

    size_t ixlo_bytes = d * sizeof(ppv_index_t); /* Bytesize of the {ixlo} vector. */
    tot_bytes += iroundup(ixlo_bytes, 8);        /* Paranoia, account for address sync. */

    return tot_bytes;
  }

kdtom_trans_t *kdtom_trans_alloc(ppv_dim_t d)
  { 
    /* Allocate the full record including all internal vectors: */
    size_t tot_bytes = kdtom_trans_node_bytesize(d); /* Total node bytesize including all vectors. */
    kdtom_trans_t *T = (kdtom_trans_t *)notnull(malloc(tot_bytes), "no mem");
    
    /* Set {pend} to the free space inside the record, past all fixed fields: */
    size_t fix_bytes = sizeof(kdtom_trans_t);  /* Size of fixed felds incl. those of {T.h}. */
    char *pend = addrsync(((char*)T) + fix_bytes, 8);

    /* Initialize the {T.h} fields, including the {T.h.size} vector: */
    kdtom_node_init((kdtom_t *)T, tot_bytes, d, &pend);
    assert(T->h.d == d);
    T->h.kind = kdtom_kind_TRANS;

    /* Allocate the variant-specific vector {T.ixlo}: */
    size_t ixlo_elsz = sizeof(ppv_index_t);
    assert(ixlo_elsz == 8);
    assert(((uint64_t)pend) % 8 == 0); /* Because we just allocated the {T.h.size} vector. */
    T->ixlo = (ppv_index_t *)kdtom_alloc_internal_vector((kdtom_t *)T, tot_bytes, d, ixlo_elsz, &pend);

    /* Check for overflow. Note that some bytes may be left over due to rounding: */
    assert (pend - (char*)T <= tot_bytes);
    
    return T;
  }

ppv_sample_t kdtom_trans_get_sample(kdtom_trans_t *T, ppv_index_t ix[])
  {
    assert(T->h.kind == kdtom_kind_TRANS);
    ppv_dim_t d = T->h.d;
    
    /* Check index vector {ix}, compute vector {jx} for {T.sub}: */
    ppv_index_t jx[d];
    for (ppv_axis_t k = 0; k < d; k++)
      { ppv_index_t ixk = ix[k];
        demand((0 <= ixk) && (ixk < T->h.size[k]), "invalid index");
        ppv_index_t lok = T->ixlo[k];
        if ((ixk < lok) || (ixk - lok >= T->sub->size[k]))
          { /* Voxel is not in {T.sub}: */
            return T->fill;
          }
        else
          { jx[k] = ixk - lok; }
      }
    /* Fetch sample from {sub} at shifted indices: */
    ppv_sample_t smp = kdtom_get_sample(T->sub, jx);
    return smp;
  }

size_t kdtom_trans_bytesize(kdtom_trans_t *T, bool_t total)
  {
    size_t node_bytes = kdtom_trans_node_bytesize(T->h.d);
    size_t tot_bytes = node_bytes;
    if (total)
      { /* Add size of subtree: */
        size_t subt_bytes = kdtom_bytesize(T->sub, TRUE);
        tot_bytes += subt_bytes;
      }
    return tot_bytes;
  }
