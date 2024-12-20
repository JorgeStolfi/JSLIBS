/* See {kdtom_split.h}. */
/* Last edited on 2024-12-05 10:32:56 by stolfi */

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
#include <kdtom_split.h>
#include <kdtom_const.h>
#include <kdtom_test.h>

#define szFMT ppv_size_t_FMT
#define ixFMT ppv_index_t_FMT
#define smpFMT ppv_sample_t_FMT

kdtom_split_t *kdtom_split_make
  ( ppv_dim_t d,
    ppv_sample_t maxsmp,
    ppv_sample_t fill,
    ppv_index_t ixlo[],
    ppv_size_t size[], 
    ppv_axis_t ax, 
    kdtom_t *sub0,
    ppv_size_t size0,
    kdtom_t *sub1,
    ppv_size_t size1
  )
  {
    assert((d > 0) && (d <= ppv_MAX_DIM));
    demand(fill <= maxsmp, "bad {fill}");
    demand(ax < d, "invalid axis {ax}");
    
    /* Check attributes of children: */
    if (sub0 != NULL)
      { demand(sub0->d == d, "incompatible {sub0.d}");
        demand(sub0->maxsmp == maxsmp, "incompatible {sub0.maxsmp}");
      }
    if (sub1 != NULL)
      { demand(sub1->d == d, "incompatible {sub1.d}");
        demand(sub1->maxsmp == maxsmp, "incompatible {sub1.maxsmp}");
      }
    
    /* Check domain geometry: */
    demand(! ixbox_is_empty(d, size), "domain is empty");
    demand(size0 >= 1, "invalid {size0}");
    demand(size0 + size1 == size[ax], "inconsistent {size0,size1}");
    demand((size0 == 0) == (sub0 == NULL), "inconsistent {size0,sub0}");
    demand((size1 == 0) == (sub1 == NULL), "inconsistent {size1,sub1}");
    
    bool_t empty_sub0_core = (sub0 == NULL) || (sub0->size[0] == 0);
    bool_t empty_sub1_core = (sub1 == NULL) || (sub1->size[0] == 0);
    
    /* The children nodes must be already clipped and relative: */
    if (! empty_sub0_core)
      { for (ppv_axis_t k = 0; k < d; k++) 
          { ppv_size_t size0k = (k == ax ? size0 : size[k]); /* Max extent of {sub0.DK}. */
            demand(sub0->ixlo[k] >= 0, "{sub0.ixlo} is negative");
            demand(sub0->ixlo[k] + sub0->size[k] <= size0k, "{sub0.size} too large");
          }
      }
    if (! empty_sub1_core)
      { for (ppv_axis_t k = 0; k < d; k++) 
          { ppv_size_t size1k = (k == ax ? size1 : size[k]); /* Max extent of {sub1.DK}. */
            demand(sub1->ixlo[k] >= 0, "{sub1.ixlo} is negative");
            demand(sub1->ixlo[k] + sub1->size[k] <= size1k, "{sub1.size} too large");
          }
      }

    /* Allocate the full record including all internal vectors: */
    size_t tot_bytes = kdtom_split_node_bytesize(d); /* Total node bytesize including all vectors. */
    kdtom_split_t *T = (kdtom_split_t *)notnull(malloc(tot_bytes), "no mem");

    /* Set {pend} to the free space inside the record, past all fixed fields: */
    size_t fix_bytes = sizeof(kdtom_split_t);  /* Size of fixed felds incl. those of {T.h}. */
    char *pend = addrsync(((char*)T) + fix_bytes, 8);

    /* Initialize the {T.h} fields, including the {T.h.size} vector: */
    kdtom_node_init((kdtom_t*)T, kdtom_kind_SPLIT, d, maxsmp, fill, ixlo, size, tot_bytes, &pend); 
    T->ax = ax;
    T->size0 = size0;
    T->sub[0] = sub0;
    T->sub[1] = sub1;
    
    return T;
  }
  
kdtom_t *kdtom_split_make_clipping
  ( ppv_dim_t d,
    ppv_sample_t maxsmp,
    ppv_sample_t fill,
    ppv_index_t ixlo[], 
    ppv_size_t size[], 
    ppv_axis_t ax, 
    ppv_index_t ixcut,
    kdtom_t *T0,
    kdtom_t *T1
  )
  { 
    bool_t debug = TRUE;
    #define dbtag "{kdtom_split_make_clipping}: "
    
    demand(T0 != NULL, "{T0} is null");
    demand(T0->d == d, "{T0.d} mismatch");
    demand(T0->maxsmp == maxsmp, "{T0.maxsmp} mismatch");

    demand(T1 != NULL, "{T1} is null");
    demand(T1->d == d, "{T1.d} mismatch");
    demand(T1->maxsmp == maxsmp, "{T1.maxsmp} mismatch");
    
    kdtom_t *T = NULL; /* Result node: */
    if (ixbox_is_empty(d, size))
      { /* Return a constant node that is everywhere {fill}: */
        if (debug) { fprintf(stderr, dbtag "core is empty\n"); }
        kdtom_const_t *Tc = kdtom_const_make(d, maxsmp, fill, NULL, NULL, fill);
        T = (kdtom_t *)Tc;
      }
    else 
      { /* Non-trivial core. */
      
        /* Clip the two sub-nodes to the sub-domains, */
        /* then translate them so that the low corner of each sub-domain is at origin: */
        
        /* Compute the low index {ixlo0} and size {size0} of {T.DK0} along axis {ax}: */
        ppv_index_t ixlo0;
        ppv_size_t size0;
        if (ixcut < ixlo[ax])
          { ixlo0 = 0; size0 = 0; }
        else if (ixcut >= ixlo[ax] + size[ax])
          { ixlo0 = ixlo[ax]; size0 = size[ax]; }
        else
          { ixlo0 = ixlo[ax]; size0 = (ppv_size_t)(ixcut - ixlo[ax]); }
        assert(size0 <= size[ax]);
        if (debug) { fprintf(stderr, dbtag "size0 = " szFMT "\n", size0); }
        
        /* Compute the low index {ixlo1} and size {size1} of {T.DK1} along axis {ax}: */
        ppv_index_t ixlo1 = ixlo0 + size0;
        ppv_size_t size1 = (ppv_size_t)(size[ax] - size0);
        
        /* Clip and shift the children nodes: */

        /* Save {ixlo[ax],size[ax]}: */
        ppv_index_t iax = ixlo[ax];
        ppv_index_t sax = size[ax];

        /* Clip core of {T0} to {T.DK0}, shift by {-ixlo0}: */
        kdtom_t *sub0 = NULL; /* {T0} clipped to {T.DK0} and shifted to be relative. */
        if (size0 != 0)
          { demand ((T0 != NULL), "{T0} cannot be null");
            ixlo[ax] = ixlo0;
            size[ax] = size0;
            sub0 = kdtom_clip_core(T0, ixlo, size);
            if (sub0 == T0) { sub0 = kdtom_clone(sub0); }  /* Can't modify {T0}. */
            for (ppv_axis_t k = 0; k < d; k++) { sub0->ixlo[k] -= ixlo[k]; }
            if (debug) { ixbox_print(stderr, d, dbtag "sub0.K =  [ ", sub0->ixlo, sub0->size, " ]\n"); }
          }
          
        /* Clip core of {T1} to {T.DK1}, shift by {-ixlo1}: */
        kdtom_t *sub1 = NULL; /* {T1} clipped to {T.DK1} and shifted to be relative. */
        if (size1 != 0)
          { demand ((T1 != NULL), "{T1} cannot be null");
            ixlo[ax] = ixlo1;
            size[ax] = size1;
            sub1 = kdtom_clip_core(T1, ixlo, size);
            if (sub1 == T1) { sub1 = kdtom_clone(sub1); }  /* Can't modify {T1}. */
            for (ppv_axis_t k = 0; k < d; k++) { sub1->ixlo[k] -= ixlo[k]; }
            if (debug) { ixbox_print(stderr, d, dbtag " sub1.K =  [ ", sub1->ixlo, sub1->size, " ]\n"); }
          }

        /* Restore {ixlo[ax],size[ax]}: */
        ixlo[ax] = iax; size[ax] = sax; 

        /* Create the record: */
        T = (kdtom_t *)kdtom_split_make(d, maxsmp, fill, ixlo, size, ax, sub0, size0, sub1, size1);
      }
    assert(T != NULL);
    return T;
    #undef dbtag
  }
  
kdtom_t *kdtom_split_simplify(kdtom_split_t *T)
  {
    /* !!! IMPLEMENT IT !!! */
    demand(FALSE, "{kdtom_split_simplify} not implemented");
  }

kdtom_split_t *kdtom_split_clone(kdtom_split_t *T)
  { 
    ppv_size_t size0 = T->size0;
    ppv_size_t size1 = T->h.size[T->ax] - size0;
    kdtom_split_t *S = kdtom_split_make
      ( T->h.d, T->h.maxsmp, T->h.fill, T->h.ixlo, T->h.size, 
        T->ax, T->sub[0], size0, T->sub[1], size1
      ); 
    return S;
  }

kdtom_t *kdtom_split_clip_core(kdtom_split_t *T, ppv_index_t ixlo[], ppv_size_t size[])
  {
    ppv_dim_t d = T->h.d;
    ppv_sample_t maxsmp = T->h.maxsmp;
    ppv_sample_t fill = T->h.fill;
    ppv_axis_t ax = T->ax;

    demand(! ixbox_is_empty(d, size), "clipping box is empty");
    assert(size != NULL);
    
    /* Obtain the two sub-trees {T0,T1}  translated to their places rel to {T}: */
    kdtom_t *T0 = kdtom_clone(T->sub[0]);
    kdtom_translate(T0, T->h.ixlo);

    kdtom_t *T1 = kdtom_clone(T->sub[1]);
    T->h.ixlo[ax] += T->size0;
    kdtom_translate(T1, T->h.ixlo);
    T->h.ixlo[ax] -= T->size0;

    /* Make the record: */
    ppv_index_t ixcut = T->h.ixlo[ax] + T->size0;
    kdtom_t *S = (kdtom_t *)kdtom_split_make_clipping(d, maxsmp, fill, ixlo, size, ax, ixcut, T0, T1);
    assert(S != NULL);
    return S;
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
    /* Check relative indices, just to be chato: */
    for (ppv_axis_t k = 0; k < d; k++) { assert((0 <= dx[k]) && (dx[k] < T->h.size[k])); }

    ppv_axis_t ax = T->ax; assert((0 <= ax) && (ax < d));
    /* Save original relative index of axis {ax}: */
    ppv_index_t dx_save = dx[ax];
    
    /* Decide where to go, and adjust index: */
    ppv_size_t size0 = T->size0; /* Size of low sub-domain on axis {ax}. */
    int32_t ksub;     /* Which subnode (0 or 1) has the desired voxel. */
    if (dx[ax] < size0)
      { ksub = 0; }
    else
      { ksub = 1; dx[ax] = dx[ax] - size0; }
    
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
        for (uint32_t i = 0;  i < 2; i++)
          { size_t subt_bytes = kdtom_bytesize(T->sub[i], TRUE);
            tot_bytes += subt_bytes;
          }
      }
    return tot_bytes;
  }

void kdtom_split_print_fields(FILE *wr, kdtom_split_t *T)
  {
    fprintf(wr, " .ax = %d", T->ax);
    fprintf(wr, " .size0 = " szFMT, T->size0);
    return;
  }

void kdtom_split_print_subtrees(FILE *wr, int32_t ind, kdtom_split_t *T)
  {
    kdtom_print_node(wr, ind, "s0", T->sub[0], TRUE);
    kdtom_print_node(wr, ind, "s1", T->sub[1], TRUE);
    return;
  }
