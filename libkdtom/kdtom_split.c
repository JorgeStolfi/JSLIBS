/* See {kdtom_split.h}. */
/* Last edited on 2021-07-12 11:20:52 by jstolfi */

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
#include <kdtom_split.h>
#include <kdtom_const.h>

kdtom_t *kdtom_split_make
  ( ppv_sample_t fill,
    ppv_index_t ixlo[], 
    ppv_size_t size[], 
    ppv_axis_t ax, 
    ppv_index_t ixcut,
    kdtom_t *sub0,
    kdtom_t *sub1
  )
  { 
    demand(sub1->d == sub0->d, "{d} mismatch");
    ppv_dim_t d = sub0->d; 
    
    demand(sub1->maxsmp = sub0->maxsmp, "{maxsmp} mismatch");
    ppv_sample_t maxsmp = sub0->maxsmp;
    
    bool_t empty = FALSE;
    for (ppv_axis_t k = 0; k < d; k++) { if (size[k] == 0) { empty = TRUE; break; } }
    
    kdtom_t *T = NULL; 
    if (empty)
      { /* Return a constant node that is everywhere {fill}: */
        kdtom_const_t *Tc = kdtom_const_make(d, maxsmp, fill, NULL, NULL, fill);
        T = (kdtom_t *)Tc;
      }
    else
      { /* Set {T} to the proper node, but translated by {-ixlo}: */
      
        /* Save {ixlo[ax],size[ax]}: */
        ppv_index_t iax = ixlo[ax];
        ppv_index_t sax = size[ax];

        /* Compute the sizes {sz0,sz1} of each sub-domain along axis {ax}: */
        ppv_size_t sz0 = (ixcut <= iax ? 0 : (ixcut >= iax + sax ? sax : ixcut - iax));
        ppv_size_t sz1 = (ixcut <= iax ? sax : (ixcut >= iax + sax ? 0 : sax - (ixcut - iax)));
        
        /* Clip the two sub-nodes to the boxes with the sizes of the sub-domains, */
        /* then translate them so that the low corner of each sub-domain is at origin: */
        size[ax] = sz0;
        kdtom_t *T0 = kdtom_clip(sub0, ixlo, size);
        for (ppv_axis_t k = 0; k < d; k++) { T0->ixlo[k] -= ixlo[k]; }
        
        size[ax] = sz1;
        ixlo[ax] = iax;
        kdtom_t *T1 = kdtom_clip(sub1, ixlo, size);
        for (ppv_axis_t k = 0; k < d; k++) { T0->ixlo[k] -= ixlo[k]; }
        
        /* Restore {ixlo[ax],size[ax]}: */
        size[ax] = sax; ixlo[ax] = iax;

        /* See if the split can be simplified: */
        if ((ixcut <= iax) || (kdtom_is_all_fill(T0, fill) && (T1->fill == fill)))
          { /* We can return just {T1} shifted by {ixlo}: */
            T = T1;
          }
        else if ((ixcut >= iax + size[ax]) || (kdtom_is_all_fill(T1, fill) && (T0->fill == fill)))
          { /* We can return just {T0} shifted by {ixlo}: */
            T = T0;
          }
        else
          { /* Try to join them info a single node: */
            T = kdtom_join_nodes(size, ax, T0, sz0, T1, sz1);
            if (T == NULL)
              { /* Must create a real {kdtom_split_t} node: */
                T = (kdtom_t *)kdtom_split_do_make(fill, size, ax, T0, sz0, T1, sz1);
              }
          }
        /* Shift {T} to {ixlo}: */
        for (ppv_axis_t k = 0; k < d; k++) { T->ixlo[k] += ixlo[k]; }
      }
    assert(T != NULL);
    return T;
  }
  
kdtom_split_t *kdtom_split_do_make
  ( ppv_sample_t fill,
    ppv_size_t size[], 
    ppv_axis_t ax, 
    kdtom_t *T0,
    ppv_size_t sz0,
    kdtom_t *T1,
    ppv_size_t sz1
  )
  {
    demand(T0->d == T1->d, "incompatible axis count");
    ppv_dim_t d = T0->d;
    assert((d > 0) && (d <= ppv_MAX_DIM));
    demand(T0->maxsmp == T1->maxsmp, "incompatible {maxsmp}");
    ppv_sample_t maxsmp = T0->maxsmp;
    demand(fill <= maxsmp, "bad {fill}");
    demand(ax < d, "invalid axis {ax}");
    demand(sz0 + sz1 == size[ax], "inconsistent {sz0,sz1}");
    
    /* Allocate the full record including all internal vectors: */
    size_t tot_bytes = kdtom_split_node_bytesize(d); /* Total node bytesize including all vectors. */
    kdtom_split_t *T = (kdtom_split_t *)notnull(malloc(tot_bytes), "no mem");

    /* Set {pend} to the free space inside the record, past all fixed fields: */
    size_t fix_bytes = sizeof(kdtom_split_t);  /* Size of fixed felds incl. those of {T.h}. */
    char *pend = addrsync(((char*)T) + fix_bytes, 8);

    /* Initialize the {T.h} fields, including the {T.h.size} vector: */
    kdtom_node_init((kdtom_t*)T, kdtom_kind_SPLIT, d, maxsmp, fill, NULL, size, tot_bytes, &pend); 
    T->ax = ax;
    T->ixcut = sz0;
    
    T->sub[0] = T0;
    T->sub[1] = T1;
    
    return T;
  }
  
kdtom_t *kdtom_split_clip(kdtom_split_t *T, ppv_index_t ixlo[], ppv_size_t size[])
  {
    /* !!! Implement !!! */
    fprintf(stderr, "!! {kdtom_split_clip} not implemented\n");
    return NULL;
  }

size_t kdtom_split_node_bytesize(ppv_dim_t d)
  {
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid dimension {d}");
    
    size_t fixf_bytes = sizeof(kdtom_split_t);   /* Fixed fields incl those of {T.h}. */
    size_t tot_bytes = iroundup(fixf_bytes, 8);  /* Account for address sync. */

    size_t ixlo_bytes = d * sizeof(ppv_index_t); /* Bytesize for {T.h.ixlo} vector. */
    tot_bytes += iroundup(ixlo_bytes, 8);        /* Paranoia, account for address sync. */

    size_t size_bytes = d * sizeof(ppv_size_t);  /* Bytesize for {T.h.size} vector. */
    tot_bytes += iroundup(size_bytes, 8);        /* Paranoia, account for address sync. */
 
    return tot_bytes;
  }  

ppv_sample_t kdtom_split_get_core_sample(kdtom_split_t *T, ppv_index_t dx[])
  {
    assert(T->h.kind == kdtom_kind_SPLIT);
    ppv_dim_t d = T->h.d;
    /* Check indices, just to be chato: */
    for (ppv_axis_t k = 0; k < d; k++) { assert((0 <= dx[k]) && (dx[k] < T->h.size[k])); }

    ppv_axis_t ax = T->ax; assert((0 <= ax) && (ax < d));
    /* Save original index of axis {ax}: */
    ppv_index_t dx_save = dx[ax];
    
    /* Decide where to go, and adjust index: */
    ppv_size_t sz0 = T->ixcut - T->h.ixlo[ax]; /* Size of low sub-domain on axis {ax}. */
    int32_t ksub;     /* Which subnode (0 or 1) has the desired voxel. */
    if (dx[ax] < sz0)
      { ksub = 0; }
    else
      { ksub = 1; dx[ax] = dx[ax] - sz0; }
    
    /* Get the sample: */
    ppv_sample_t smp = kdtom_get_sample(T->sub[ksub], dx);
    
    /* Restore original index vector: */
    dx[ax] = dx_save;
    return smp;
  }

size_t kdtom_split_bytesize(kdtom_split_t *T, bool_t total)
  {
    size_t node_bytes = kdtom_split_node_bytesize(T->h.d);
    size_t tot_bytes = node_bytes;
    if (total)
      { /* Add size of subtrees: */
        for (int32_t i = 0; i < 2; i++)
          { size_t subt_bytes = kdtom_bytesize(T->sub[i], TRUE);
            tot_bytes += subt_bytes;
          }
      }
    return tot_bytes;
  }
