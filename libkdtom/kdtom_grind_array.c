/* See {kdtom_grind_array.h}. */
/* Last edited on 2021-07-19 05:19:54 by jstolfi */

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

#include <kdtom_grind_array.h>

kdtom_t *kdtom_grind_array(ppv_array_t *A, ppv_sample_t fill)
  { 
    ppv_dim_t d = A->d;
    ppv_sample_t maxsmp = A->maxsmp;
    
    /* General parameters: */
    size_t split_bytes = kdtom_split_node_bytesize(d);
    size_t const_bytes = kdtom_const_node_bytesize(d);
    size_t array_bytes = kdtom_array_node_bytesize(d);
    
    auto bool_t small_enough(ppv_array_t *B);
      /* Returns true if and only if {B} as a single array node
        would take less memory than a k-d-tree of it. */
        
    auto kdtom_array_t *shared_array_node(ppv_array_t *B);
      /* Returns a {kdtom_array_t} node that is just a copy
        of {B}'s descriptor, sharing same sample storage area
        and sample indexing parameters. */
        
    auto ppv_axis_t choose_split_axis(ppv_array_t *B);
      /* Finds the axis along which {B} has the largest extent. */
        
    auto void split_array(ppv_array_t *B, ppv_axis_t ax, ppv_array_t **B0P, ppv_array_t **B1P);
      /* Splits the array {B} into two approximate halves by a plane
        perpendicular to {ax}. Returns the addresses of the descriptors
        of the two halves in {*B0P,*B1P} */
        
    auto kdtom_t *grind(ppv_array_t *B);
      /* Same as {kdtom_grind_array(B)} but slightly more efficient.
        Assumes that {B} is not empty. */
      
    kdtom_t *T;
    if (ppv_is_empty(A))
      { T = (kdtom_t*)kdtom_const_make(d, maxsmp, fill, NULL, NULL, fill); }
    else
      { T = grind(A); }
    return T;
      
     /* Internal procedures */
     
    kdtom_t *grind(ppv_array_t *B)
      { 
        /* Is it constant? */
        ppv_sample_t smp_min, smp_max;
        ppv_sample_range(B, &smp_min, &smp_max);
        assert(smp_max <= maxsmp);
        if (smp_min == smp_max)
          { return (kdtom_t*)kdtom_const_make(d, maxsmp, fill, NULL, B->size, smp_max); }

        /* Is it worth splitting? */
        ppv_axis_t ax = choose_split_axis(B);
        if ((B->size[ax] <= 1) || small_enough(B))
          { return (kdtom_t*)shared_array_node(B); }

        /* Try splitting: */
        ppv_array_t *B0;
        ppv_array_t *B1;
        split_array(B, ax, &B0, &B1);
        ppv_size_t size0 = B0->size[ax];
        ppv_size_t size1 = B1->size[ax];

        kdtom_t *T0 = grind(B0);
        kdtom_t *T1 = grind(B1);

        kdtom_t *T = NULL;
        if ((T0->kind == kdtom_kind_ARRAY) && (T1->kind == kdtom_kind_ARRAY))
          { /* Both parts are array nodes; just use a single array node: */
            T = (kdtom_t*)shared_array_node(B);
          }
        else 
          { T = (kdtom_t*)kdtom_split_make(d, maxsmp, fill, NULL, B->size, ax, T0, size0, T1, size1); }
        return T;
      }
        
    bool_t small_enough(ppv_array_t *B)
      {
        /* Compute size of storage area needed for the voxels: */
        ppv_sample_count_t npos = ppv_compute_npos_steps(d, B->size, NULL);
        size_t tot_voxel_bytes = ppv_tot_sample_bytes(npos, B->bps, B->bpw);
        /* Splitting will take at least one split node and two const nodes: */
        return (array_bytes + tot_voxel_bytes < split_bytes + 2*const_bytes); 
      }
        
    kdtom_array_t *shared_array_node(ppv_array_t *B)
      {
        kdtom_array_t *Ta = kdtom_array_make(B, fill);
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
        ppv_size_t size0 = B->size[ax] / 2;
        assert(size0 >= 1);
        ppv_array_t *B0 = ppv_array_clone(B);
        ppv_crop(B0, ax, 0, size0);
        ppv_array_t *B1 = ppv_array_clone(B);
        ppv_crop(B1, ax, size0, B->size[ax]-size0);
        (*B0P) = B0;
        (*B1P) = B1;
        return;
      }

   }
            
