/* See voxb_erolate.h */
/* Last edited on 2021-07-09 01:04:05 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <jsmath.h>
#include <indexing.h>
#include <indexing_io.h>
#include <affirm.h>
#include <ppv_array.h>

#include <voxb_erolate.h>

void voxb_erolate_with_brush_big(ppv_array_t *A, ppv_brush_t *b, bool_t erode);
 /* Same as {voxb_erolate_with_brush}, but uses a sliding buffer
   technique to avoid allocating a whole copy of {A} as a work area.
   
   It uses only {O(N1*H)} extra space where {N1} is the number of voxels
   in any slice of {A} perpedicular to axis 0 (that is, voxels whose
   index {ix[0]} is fixed) and {H} is the number of similar slices of
   the brush. The time is roughly proportional to {N*M} where {N} is the
   total number of pixels of {A} and {M} is the number of voxels in the
   brush. */

void voxb_erolate_with_brush_small(ppv_array_t *A, ppv_brush_t *b, bool_t erode);
 /* Same as {voxb_erolate_with_brush},but always uses a
   whole copy of {A} as a work area.
   
   The extra space is thus {O(N)} and the running time is roughly
   proportional to {N*M} where {N} is the total number of pixels of {A}
   and {M} is the number of voxels in the brush. */

ppv_sample_t voxb_erolate_compute_result
  ( const ppv_index_t ixA[],
    ppv_array_t *T, 
    ppv_index_t buf_lo,
    ppv_index_t buf_hi,
    ppv_size_t sz0A,
    ppv_brush_t *b,
    ppv_sample_t vexp
  );
  /* Computes the value {vres} of the sample with indices {ixA[0..d-1]}
    of the eroded/dilated version of some sample array {A},
    using the brush {b} as structuring element.  
    
    Assumes that the slices of the original array {A} perpendicular to
    axis {0} are indexed {0..sz0A-1}. The buffer {T} contains the
    original voxel values of the slices {buf_lo..buf_hi} of {A}. Those
    should include all the slices of {A} needded to compute the result.
    The slice {k0} of the original {A}, if it is in that ange, is slice
    {k0 % T->size[0]} of {T}.
    
    The result will be {vexp} if any of the original {A} voxels in the
    set of brush indices {V(b)} translated by {ixA} is {vexp}, otherwise
    it will be {1-vexp}. */

#define voxb_erolate_BIG_SIZE (200*1000*1000)
  /* Size (in bits) when starts to use sliding buffer technique */

void voxb_erolate_with_brush(ppv_array_t *A, ppv_brush_t *b, bool_t erode)
  {
    ppv_dim_t d = A->d;
    demand(ppv_brush_dimension(b) == d, "incompatible dimensions");   /* Could be smaller? */
    demand(A->bps == 1, "sample values are not binary");
    
    ppv_sample_count_t nvA = ppv_sample_count(A, TRUE);
    if (nvA == 0)
      { /* Empty array, nothing to do: */ return; }
    else if ((d >= 2) && (nvA >= voxb_erolate_BIG_SIZE/A->bps))
      { voxb_erolate_with_brush_big(A, b, erode); }
    else
      { voxb_erolate_with_brush_small(A, b, erode); }

    return;
  }
      
void voxb_erolate_with_brush_big(ppv_array_t *A, ppv_brush_t *b, bool_t erode)
  { 
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "entering {voxb_erolate_with_brush_big}\n"); }
    ppv_dim_t d = A->d;
    assert(d >= 2);
    /* To avoid creating a copy of the full array {A}, we process it 
      one /slice/  at a time, where a slice is all the voxels with the same
      index {ix[0]}. That common index is the /slice index/.  
      
      We use a buffer {T} to hold the original contents of the lines of
      {A} that are needed to compute one slice of the result. The buffer
      is a {ppv_array_t} with same dims as {A} except that it has
      only {sz0T} slices in the direction of axis 0. Namely, when
      computing the slice of {A} with index {ix0}, the buffer has the
      original contents of all the slices of {A} from {ix0+dx0lo} to
      {ix0+dx0hi}, where {dx0lo} and {ix0hi} are the minimum and maximum
      relative indices along axis 0 in the brush {b}.
      
      The buffer slices are used in circular fashion. Namely, if the
      slice {k0} of {A} is stored in the buffer {T}, it is stored as
      slice {k0 % T.size[0]} of {T}. */
    
    auto ppv_array_t *alloc_buffer(ppv_size_t sz0);
      /* Allocates a buffer to hold {sz0} slices of {A}. */
    
    /* Allocate the work buffer: */
    ppv_index_t dxlo[d], dxhi[d];
    ppv_brush_index_ranges(b, dxlo, dxhi);
    ppv_size_t sz0b = (ppv_size_t)(dxhi[0] - dxlo[0] + 1); /* Number of slices in the brush. */
    ppv_size_t sz0T = sz0b; /* Number of slices in the buffer. */
    ppv_array_t *T = alloc_buffer(sz0T);
    assert(T->size[0] == sz0T);
    
    auto bool_t copy_voxel_A_to_T(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pT, ppv_pos_t pC);
      /* Copies a voxel from {A} to {T}, given its positions {posA} and {posT}
        in the two arrays. The index vector {ix} and the position {posC} are ignored.
        Always returns {FALSE}. */
    
    auto void copy_A_slice_to_T(ppv_index_t ix0);
      /* Copies slice {ix0} of {A} into the corresponding slice of {T},
        namely slice {ix0 % sz0T}.  If {ix0} is not in {0..A.size[0]-1}, does nothing. */
    
    /* The variables {ixA0buf_lo} and {ixA0buf_hi} keep the min and max of the
      indices of the slices of {A} that have been copied to {T}.
      If there are no slices, {ixA0buf_lo} must be 0 and and {ixA0buf_hi}
      must be {-1}. */
    ppv_index_t ixA0buf_lo = 0, ixA0buf_hi = -1;  /* Slices of {A}  currently in {T}. */

    ppv_sample_t vexp = (erode ? 0 : 1);    /* Sample value to be propagated. */
    
    /* When computing the voxel with index vector {ix} of the result,
      we need the original voxels of {A} with index vectors
      {ix+dx} (if they exist) where each {dx} is an index tuple of the 
      brush {b}. */

    for (ppv_index_t ix0 = 0; ix0 < A->size[0]; ix0++)
      { 
        /* Compute slice {ixA[0]} of the result and store it into {A}: */
        if (debug) { fprintf(stderr, "computing slice %ld of result\n", ix0); }
        
        /* We need lines {ix0+dx[0]} for {dx[0]} in {dxlo[0]..dxhi[0]}. */
        /* Ensure those slices are in {T}: */
        while (TRUE)
          { ppv_index_t nbuf = ixA0buf_hi - ixA0buf_lo + 1; /* Number of lines in buffer. */
            if ((nbuf > 0) && (ixA0buf_hi >= ix0+dxhi[0])) { break; }
            ppv_index_t ixA0cp = (nbuf <= 0 ? 0 : ixA0buf_hi + 1); /* Next line to copy to {T}. */
            if ((nbuf > 0) && (ixA0cp - ixA0buf_lo >= sz0T))
              { /* Buffer is full -- drop the last buffered slice: */
                if (debug) { fprintf(stderr, "  dropping slice %ld\n", ixA0buf_lo); }
                ixA0buf_lo++; 
                assert(ixA0cp - ixA0buf_lo < sz0T);
              }
            if (debug) { fprintf(stderr, "  loading slice %ld\n", ixA0cp); }
            copy_A_slice_to_T(ixA0cp);
            ixA0buf_hi = ixA0cp; 
            if (nbuf <= 0) { ixA0buf_lo = ixA0buf_hi; }
          }
        if (debug) { fprintf(stderr, "  slice range in buffer = {%ld..%ld}\n", ixA0buf_lo, ixA0buf_hi); }
        assert(imax(0,ix0+dxlo[0]) >= ixA0buf_lo);
        assert(imin(ix0+dxhi[0], A->size[0]-1) <= ixA0buf_hi);
            
        /* Compute slice {ix[0]} of the answer from voxels in {T}, store in {A} */
        ppv_index_t ixA[d];
        ixA[0] = ix0; for (ppv_axis_t ax = 1; ax < d; ax++) { ixA[ax] = 0; }
        ppv_pos_t posA = ppv_sample_pos(A, ixA);
        while (TRUE) 
          { ppv_sample_t vres = voxb_erolate_compute_result(ixA, T, ixA0buf_lo, ixA0buf_hi, A->size[0], b, vexp);
            ppv_set_sample_at_pos(A->el, A->bps, A->bpw, posA, vres);
            /* Advance {ixA} to the next index vector in lex order: */
            bool_t done = ix_next(d, ixA, A->size, ix_order_L, A->step, &posA, NULL, NULL, NULL, NULL);
            if (done || (ixA[0] != ix0)) { break; }
          } 
      }
    /* Reclaim the buffer's storage: */
    free(T->el);
    free(T);
    
    return;
    
    /* Local procedures: */
    ppv_array_t *alloc_buffer(ppv_size_t sz0)
      { ppv_size_t szT[d];
        for (ppv_axis_t i = 0; i < d; i++) { szT[i] = A->size[i]; }
        szT[0] = sz0;
        ppv_array_t *T = ppv_array_new(d, szT, A->maxsmp);
        return T;
      }
    
    bool_t copy_voxel_A_to_T(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pT, ppv_pos_t pC)
      { ppv_sample_t av = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pA);
        ppv_set_sample_at_pos(T->el, T->bps, T->bpw, pT, av);
        return FALSE;
      }

    void copy_A_slice_to_T(ppv_index_t ix0)
      { assert(sz0T == T->size[0]);
        if ((ix0 >= 0) && (ix0 < A->size[0]))
          { ppv_array_t *Az = ppv_slice(A, 0, ix0);
            ppv_array_t *Tz = ppv_slice(T, 0, ix0 % sz0T);
            (void) ppv_enum(copy_voxel_A_to_T, FALSE, Az, Tz, NULL); 
            free(Az);
            free(Tz);
          }
      }
  }

      
void voxb_erolate_with_brush_small(ppv_array_t *A, ppv_brush_t *b, bool_t erode)
  { 
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "entering {voxb_erolate_with_brush_small}\n"); }
    
    ppv_dim_t d = A->d;
    
    /* We just copy {A} to a buffer array {T} and compute the result from {T}. */
    
    ppv_array_t *T = ppv_array_new(d, A->size, A->maxsmp);
    ppv_array_assign(T, A);
    
    ppv_index_t ixA0buf_lo = 0;             /* Low slice of {A} in buffer. */
    ppv_index_t ixA0buf_hi = A->size[0]-1;  /* High slice of {A} in buffer. */

    ppv_sample_t vexp = (erode ? 0 : 1);    /* Sample value to be propagated. */
    
    ppv_index_t ix0_last = -1; /* Last value of {ix[0]}  seen. */
    if (debug) { fprintf(stderr, "slices: "); }
    auto bool_t erolate1(const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC);
      /* Computes element {ix} of the result from the original {A} voxels saved in {T},
        and stores that result in {A} at those indices. */
    bool_t stop = ppv_enum(erolate1, FALSE, A, T, NULL);
    assert(! stop);
    if (debug) { fprintf(stderr, "\n"); }
    
    return;
    
    bool_t erolate1(const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC)
     { ppv_sample_t vres =  voxb_erolate_compute_result(ix, T, ixA0buf_lo, ixA0buf_hi, A->size[0], b, vexp);
       ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pA, vres);
       if (debug && (ix[0] != ix0_last))
         { fprintf(stderr, " %ld", ix[0]);
           ix0_last = ix[0];
         }
       return FALSE;
     }

    /* Reclaim the buffer's storage: */
    free(T->el);
    free(T);
    
    return;
  }

ppv_sample_t voxb_erolate_compute_result
  ( const ppv_index_t ixA[],
    ppv_array_t *T, 
    ppv_index_t buf_lo,
    ppv_index_t buf_hi,
    ppv_size_t sz0A,
    ppv_brush_t *b,
    ppv_sample_t vexp
  )
  {
    bool_t debug = FALSE;

    ppv_dim_t d = T->d;

    if (debug) 
      { ix_print_indices (stderr, "computing A[", d, ixA, 0, ",", "]");
        fprintf(stderr, " vexp = %d\n", vexp);
      }

    ppv_size_t sz0T = T->size[0];
    
    /* Enumerate brush voxels, stop if finds a {vexp}: */
    auto bool_t checkvoxel(const ppv_index_t dx[]);
    bool_t found = ppv_brush_enum(checkvoxel, b);
    ppv_sample_t vres = (found ? vexp : 1 - vexp);
    return vres;
    
    bool_t checkvoxel(const ppv_index_t dx[])
      { ppv_index_t ixT[d]; /* Index of translated brush pixel {T}. */
        for (ppv_axis_t ax = 0; ax < d; ax++)
          { ixT[ax] = ixA[ax] + dx[ax];
            ppv_size_t szA = (ax == 0 ? sz0A : T->size[ax]);
            if ((ixT[ax] < 0) || (ixT[ax] >= szA)) { return FALSE; }
          }
        assert((buf_lo <= ixT[0]) && (ixT[0] <= buf_hi));
        ixT[0] = ixT[0] % sz0T;
        ppv_sample_t vget = ppv_get_sample(T, ixT);
        return (vget == vexp);
      }
  }
