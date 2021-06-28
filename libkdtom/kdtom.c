/* See {kdtom.h}. */
/* Last edited on 2021-06-28 08:41:23 by jstolfi */

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
#include <kdtom_array.h>
#include <kdtom_const.h>
#include <kdtom_ixmap.h>

ppv_sample_t kdtom_get_sample(kdtom_t *T, ppv_index_t ix[])
  {
    /* Check index bounds.  Early warning is good. */
    for (ppv_axis_t ax = 0; ax < T->d; ax++)
      { demand((0 <= ix[ax]) && (ix[ax]< T->size[ax]), "invalid voxel index"); }
      
    /* Switch to proper handler: */
    switch (T->kind)
      { 
        case kdtom_kind_ARRAY:
          return kdtom_array_get_sample((kdtom_array_t *)T, ix);
        case kdtom_kind_CONST:
          return kdtom_const_get_sample((kdtom_const_t *)T, ix);
        case kdtom_kind_SPLIT:
          return kdtom_split_get_sample((kdtom_split_t *)T, ix);
        case kdtom_kind_IXMAP:
          return kdtom_ixmap_get_sample((kdtom_ixmap_t *)T, ix);
        default:
          assert(FALSE);
      }
  }

size_t kdtom_bytesize(kdtom_t *T, bool_t total)
  {
    switch (T->kind)
      { 
        case kdtom_kind_ARRAY:
          return kdtom_array_bytesize((kdtom_array_t *)T, total);
        case kdtom_kind_CONST:
          return kdtom_const_bytesize((kdtom_const_t *)T);
        case kdtom_kind_SPLIT:
          return kdtom_split_bytesize((kdtom_split_t *)T, total);
        case kdtom_kind_IXMAP:
          return kdtom_ixmap_bytesize((kdtom_ixmap_t *)T, total);
        default:
          assert(FALSE);
      }
  }

char *kdtom_alloc_internal_vector(kdtom_t *T, size_t sztot, ppv_dim_t d, size_t elsz, char **pendP)
  {
    char *pend = (*pendP);
    
    /* Check where the address {pend} is properly aligned for the element size: */
    demand(((uint64_t)pend) % elsz == 0, "{*pendP} sync error");
    
    /* Check for enough space, alocate the vector, and update {*pendP}: */
    size_t vec_bytes = d*elsz;
    assert (sztot - (pend - (char*)T) >= vec_bytes);
    char *vaddr = pend; pend += vec_bytes;
    (*pendP) = pend;
    
    return vaddr;
  }

void kdtom_node_init(kdtom_t *T, size_t sztot, ppv_dim_t d, char **pendP)
  {
    /* Check for clobbering of the fixed head fields: */
    size_t fix_bytes = sizeof(kdtom_t);        /* Fixed fields, incl. {size} pointer. */
    demand((*pendP) - (char *)T >= fix_bytes, "bad {*pendP}");
    
    /* Allocate the {T.size} vector: */
    size_t elsz = sizeof(ppv_size_t);
    T->size = (ppv_size_t *)kdtom_alloc_internal_vector(T, sztot, d, elsz, pendP);
    
    /* Set the {d} field: */
    T->d = d;
  }

kdtom_t *kdtom_grind_array(ppv_array_t *A)
  { 
    auto bool_t small_enough(ppv_array_t *B);
      /* Returns true if and only if {B} as a single array node
        would take less memory than a k-d-tree of it. */
        
    auto kdtom_array_t *shared_array_node(ppv_array_t *B);
      /* Returns a {kdtom_array_t} node that is just a copy
        of {B}'s descriptor, sharing same sample storage area
        and sample indexng parameters. */
        
    auto kdtom_array_t *copied_array_node(ppv_array_t *B);
      /*  Returns a {kdtom_array_t} node that has the same size as {B}
        but a a newly alocated sample storage area, after copying into
        it the sample values of the voxels of {B}. */
        
    auto ppv_axis_t choose_split_axis(ppv_array_t *B);
      /* Finds the axis along which {B} has the largest extent. */
        
    auto void split_array(ppv_array_t *B, ppv_axis_t ax, ppv_array_t **B0P, ppv_array_t **B1P);
      /* Splits the array {B} into two approximate halves by a plane
        perpendicular to {ax}. Returns the addresses of the descriptors
        of the two halves in {*B0P,*B1P} */
        
    auto kdtom_t *grind(ppv_array_t *B);
      /* Same as {kdtom_grind_array(B)} but slightly more efficient.
        Assumes that {B} is not empty. */
      
    /* General parameters: */
    ppv_dim_t d = A->d;
    ppv_nbits_t bps = A->bps;
    size_t split_sz = kdtom_split_node_bytesize(d);
    size_t const_sz = kdtom_const_node_bytesize(d);
    size_t array_sz = kdtom_array_node_bytesize(d);
    ppv_nbits_t bpw = ppv_best_bpw(bps);
    
    kdtom_t *T;
    if (ppv_is_empty(A))
      { T = (kdtom_t*)kdtom_const_make(d, bps, A->size, 0); }
    else
      { T = grind(A); }
    return T;
      
     /* Internal procedures */
     
    kdtom_t *grind(ppv_array_t *B)
      { 
        /* Is it constant? */
        ppv_sample_t vmin, vmax;
        ppv_sample_range(B, &vmin, &vmax);
        if (vmin == vmax)
          { return (kdtom_t*)kdtom_const_make(d, bps, B->size, vmin); }

        /* Is it worth splitting? */
        ppv_axis_t ax = choose_split_axis(B);
        if ((B->size[ax] <= 1) || small_enough(B))
          { return (kdtom_t*)shared_array_node(B); }

        /* Try splitting: */
        ppv_array_t *B0;
        ppv_array_t *B1;
        split_array(B, ax, &B0, &B1);
        kdtom_t *T0 = kdtom_grind_array(B0);
        kdtom_t *T1 = kdtom_grind_array(B1);

        /* Was it worth splitting? */
        if ((T0->kind == kdtom_kind_ARRAY) && (T1->kind == kdtom_kind_ARRAY))
          { /* Both parts are array nodes; just use a single array node: */
            free(T0); free(T1);
            return (kdtom_t*)shared_array_node(B);
          }
        else
          { /* If both sides are constant they must be distinct constants. */
            /* Replace shared-element arrays by deep copies:*/
            if (T0->kind == kdtom_kind_ARRAY) 
              { T0 = (kdtom_t *)copied_array_node(B0); }
            if (T1->kind == kdtom_kind_ARRAY) 
              { T1 = (kdtom_t *)copied_array_node(B1); }
            return (kdtom_t*)kdtom_split_make(ax, T0,T1);
          }
        }
        
    bool_t small_enough(ppv_array_t *B)
      {
        /* Compute size of storage area needed for the voxels: */
        ppv_sample_count_t npos = ppv_compute_npos_steps(d, B->size, NULL);
        size_t voxel_sz = ppv_tot_sample_bytes(npos, bps, bpw);
        /* Splitting will take at least one split node and two const nodes: */
        return (array_sz + voxel_sz < split_sz + 2*const_sz); 
      }
        
    kdtom_array_t *shared_array_node(ppv_array_t *B)
      {
        kdtom_array_t *Ta = kdtom_array_make(B);
        return Ta;
      }
        
    kdtom_array_t *copied_array_node(ppv_array_t *B)
      {
        ppv_array_t *C = ppv_array_new(d, B->size, bps, bpw);
        ppv_array_assign(C, B);
        kdtom_array_t *Ta = kdtom_array_make(C);
        return Ta;
      }
        
    ppv_axis_t choose_split_axis(ppv_array_t *B)
      {
        ppv_size_t szmax = 0;
        ppv_axis_t axmax = 0;
        for (ppv_axis_t k = 0; k < d; k++)
          { ppv_size_t szk = B->size[k];
            assert(szk > 0);   /* Empty array must have been caught before. */
            if (szk > szmax) { szmax = szk; axmax = k; }
          }
        return axmax;
      }
        
    void split_array(ppv_array_t *B, ppv_axis_t ax, ppv_array_t **B0P, ppv_array_t **B1P)
      {
        ppv_size_t sz0 = B->size[ax] / 2;
        assert(sz0 >= 1);
        ppv_array_t *B0 = ppv_array_clone(B);
        ppv_crop(B0, ax, 0, sz0);
        ppv_array_t *B1 = ppv_array_clone(B);
        ppv_crop(B1, ax, sz0, B->size[ax]-sz0);
        (*B0P) = B0;
        (*B1P) = B1;
        return;
      }

   }
            
