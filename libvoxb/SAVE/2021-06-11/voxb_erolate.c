/* See voxb_erolate.h */
/* Last edited on 2021-06-11 17:06:19 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <i2.h>
#include <affirm.h>
#include <ppv_array.h>

#include <voxb_erolate.h>

#define NAX (ppv_array_NAXES)
  /* Number of indices in a {ppv_array_t}. */

typedef struct voxb_erolate_ball_t
  { int32_t irad;      /* Planes in ball are indexed from {-irad..+irad}. */
    int32_t *nv;       /* Number of ball voxels in each plane. */
    i2_t **ind;         /* Indices {(ix,iy)} of those voxes in each plane. */
    ppv_step_t **step; /* Relative positions of those voxels in each plane of {T}. */
  } voxb_erolate_ball_t;
  /* Description of a digital ball in a specific {ppv_array_t}.  
    
    The ball consists of every voxel with indices {dz,dy,dx} in the
    standard bi-infinite integer voxel grid whose center is at most a
    certain radius {rad} away from the center of voxel {0,0,0}. Those
    voxels have {dz}, {dy}, and {dx} in {-irad..+irad}.
    
    The ball is assumed to be symmetric in {Z}, so only planes with {dz}
    in {0..irad} are stored. The value {nv[dz]} is the number of voxels
    in plane {dz} of the digital ball. Element {ind[dz]} is a vector of
    size {nv[dz]} with all integer pairs {(dz,dy)} of those voxels (note
    the order). Element {step[dz]} is the bit address displacement in
    the associated array {T} from a generic voxel with indices
    {kz,ky,kx} othe voxel with indices {kz,ky+dy,kx+dx}. Note that the
    displacement is within the plane {dz} of {T}. */ 

voxb_erolate_ball_t voxb_erolate_get_ball(double rad, ppv_step_t step[]);
  /* Computes the indices of all voxels of a digital ball of radius
    {rad} (which must be positive), and the respective position
    displacements in a {ppv_array_y} with the given {step} vector. The
    bounds of the array are ignored in this computaton, so each index
    may be negative or beyond the last valid index. */

void voxb_erolate_free_ball(voxb_erolate_ball_t b);
  /* Reclaims the storage used by the internal tables of {b} (but not of {b} itself). */

ppv_sample_t voxb_erolate_compute_result
  ( ppv_array_t *T, 
    ppv_size_t nzA,
    ppv_index_t kz,
    ppv_index_t ky,
    ppv_index_t kx, 
    voxb_erolate_ball_t *b,
    ppv_sample_t vexp
  );
  /* Computes the value {vres} of the sample with indices {kz,ky,kx} (in that order)
    of the eroded/dilated version of some {ppv_array_t} {A}.  
    
    Assumes that the planes of the original array {A} are indexed
    {0..nzA-1}. Assumes also that the plane {kz+dz} of the original {A},
    if it exists, is stored in plane {(kz+dz) \mod nzT} of the buffer
    array {T}, where {nzT} is the height of {T}. Assumes that the
    relative indices and displacements of the voxels in the digital ball
    are stored in {b}.
    
    The result will be {vexp} if any of the original {A} voxels in the
    ball centered at {kz,ky,kx} is {vexp}, otherwise it will be
    {1-vexp}. */

void voxb_erolate_with_ball(ppv_array_t *A, double rad)
  {
    demand(A->bps == 1, "sample values are not binary");
    
    if (rad == 0) { return; }
    
    /* Dimensions of {A}: */
    ppv_size_t nzA = A->size[0];
    ppv_size_t nyA = A->size[1];
    ppv_size_t nxA = A->size[2];
    
    auto ppv_array_t *alloc_buffer(ppv_size_t nzB);
      /* Allocates a buffer to hold {nzB} planes of {A}. */
    
    /* Allocate the work buffer: */
    int32_t irad  = (int32_t)ceil(fabs(rad)); /* Radius of digital ball. */
    ppv_sample_t vexp = (rad > 0 ? 0 : 1);    /* Sample value to be propagated. */

    int32_t nzT = 2*irad + 1; /* Number of planes in the buffer. */
    ppv_array_t *T = alloc_buffer(nzT);
    
    /* The buffer is a {ppv_array_t} with same {X,Y} dims as {A} and only
      {nzT} planes in the direction {Z}. 
      
      When computing plane {kz} of the result, the original planes of {A} from
      {kz-irad} to {kz+irad} (inclusive both) are stored in the buffer.
      Specifically, the original plane {kz+dz} of {A} is stored in plane
      {(kz+dz) % nzT} of {T}, for {dz} in {-irad..+irad}. */
    
    auto void copy_voxel_A_to_T(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pT, ppv_pos_t pC);
      /* Copies a voxel from {A} to {T}, given its positions {posA} and {posT}
        in the two arrays. The index vector {ix} and the position {posC} are ignored. */
    
    auto void copy_A_to_T(int32_t rz);
      /* Copies plane {rz} of {A} into the corresponding plane of {T},
        namely plane {rz % nzT}.  If {rz} is not in {0..nzA-1}, does nothing. */

    voxb_erolate_ball_t b = voxb_erolate_get_ball(fabs(rad), T->step);
    assert(b.irad == irad);
    
    /* The arrays {stepb} and {nvb} will have {irad+1} elements,
      and the array {stepb[dz]} will have {nvb[dz]} elements, for {dz} in {0..irad}.
      
      When computing the voxel with indices {(kz,ky,kx)} of the result,
      we need the original voxels of {A} with indices
      {(kz+dz,ky+dy,kx+dx)} (if they exist) where {(dz,dy,dx)} is a
      relative index displacement from the central voxel of the digital
      ball to each voxel of the ball. The displacement {dz} ranges in
      {-irad..+irad}, and the set of pairs {dy,dx} depends on {dz}.
      
      That will be the voxel with indices {((kz+dz)%nzT,ky+dy,kx+dx)} of
      {T}. Its position in {T} will be

        {T.base + ((kz+dz)%nzT)*T.step[0] + (ky+dy)*T.step[1] + (kx+dx)*T.step[2]}
        
      which is the position in {T} of the voxel with indices {((kz+dz)%nzT,ky,kx)}
      plus the displacement {dy*T.step[1] + dx*T.step[2]}. 
      
      This displacement is obtained as {b->step[iabs(dz)][iv]} where {iv} ranges
      in {0..b->nv[iabs(dz)]-1}.  The corresponding voxel index increment 
      pairs {(dx,dy)} for each of those pixels are stored in {indb[iabs(dz)][iv]}. */

    /* Preload the buffer {T} with planes {0..irad-1} of {A}: */
    for (int32_t dz = 0; dz < irad; dz++) { copy_A_to_T(dz); }
    
    /* Index vector to access voxels of {A}: */
    ppv_index_t indA[ppv_array_NAXES]; ppv_index_clear(indA);
                
    for (int32_t kz = 0; kz < nzA; kz++)
      { /* Now the original planes {kz-irad .. kz+irad-1} of {A} (if exist) are in {T}. */
        
        /* Ensure plane {kz+irad} of {A} is in {T}  too: */
        copy_A_to_T(kz + irad);
        
        /* Compute plane {kz} of the answer from voxels in {T}, store in {A} */
        for (int32_t ky = 0; ky < nyA; ky++)
          { for (int32_t kx = 0; kx < nxA; kx++)
              { ppv_sample_t vres = voxb_erolate_compute_result(T, nzA, kz,ky,kx, &b, vexp);
                indA[0] = kz; indA[1] = ky; indA[2] = kx;
                ppv_set_sample(A, indA, vres);
              }
          }
      }
    /* Reclaim the buffer's storage: */
    free(T->el);
    free(T);
    voxb_erolate_free_ball(b);
    
    return;
    
    /* Local procedures: */
    ppv_array_t *alloc_buffer(ppv_size_t nzB)
      { ppv_size_t szT[ppv_array_NAXES];
        for (int32_t i = 0; i < ppv_array_NAXES; i++) { szT[i] = A->size[i]; }
        szT[0] = nzB;
        ppv_array_t *T = notnull(malloc(sizeof(ppv_array_t)), "no mem");
        (*T) = ppv_new_array(szT, A->bps, A->bpw);
        return T;
      }
    
    void copy_voxel_A_to_T(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pT, ppv_pos_t pC)
      { ppv_sample_t av = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pA);
        ppv_set_sample_at_pos(T->el, T->bps, T->bpw, pT, av);
      }

    void copy_A_to_T(int32_t rz)
      { assert(nzT == T->size[0]);
        if ((rz >= 0) && (rz < nzA))
          { ppv_array_t Az = (*A); ppv_slice(&Az, 0, rz);
            ppv_array_t Tz = (*T); ppv_slice(&Tz, 0, rz % nzT);
            ppv_enum(copy_voxel_A_to_T, FALSE, &Az, &Tz, NULL); 
          }
      }
  }

voxb_erolate_ball_t voxb_erolate_get_ball(double rad, ppv_step_t step[])
  { 
    bool_t debug = TRUE;
    
    int32_t irad  = (int32_t)ceil(fabs(rad)); /* Radius of digital ball. */
    int32_t nzb = irad+1; /* Number of planes stored. */
    
    ppv_step_t stepy = step[1];
    ppv_step_t stepx = step[2];

    voxb_erolate_ball_t b;
    b.irad = irad;
    b.nv = notnull(malloc(nzb*sizeof(int32_t)), "no mem");
    b.ind = notnull(malloc(nzb*sizeof(i2_t *)), "no mem");
    b.step = notnull(malloc(nzb*sizeof(ppv_step_t *)), "no mem");
    if (debug) { fprintf(stderr, "voxels per plane = "); }
    int32_t nvtot = 0; /* Total number of voxels in ball. */
    for (int32_t dz = 0;  dz < nzb; dz++)
      { double radz = sqrt(rad*rad - dz*dz) + 1.0e-8; /* Radius of ball on this plane. */
        int32_t iradz = (int32_t)ceil(fabs(radz)); /* Max radius of digital ball on this plane. */
        int32_t maxnv = (2*iradz+1)*(2*iradz+1); /* Max voxels of ball on this plane. */
        b.ind[dz] = notnull(malloc(maxnv*sizeof(i2_t)), "no mem");
        b.step[dz] = notnull(malloc(maxnv*sizeof(ppv_step_t)), "no mem");
        int32_t nv = 0; /* Number of voxels in plane. */
        for (int32_t dy = -iradz; dy <= +iradz; dy++)
          { for (int32_t dx = -iradz; dx <= +iradz; dx++)
             { /* The following computaton should be exact for reasonable {rad}: */
               double d2 = dx*(double)dx + dy*(double)dy + dz*(double)dz; 
               if (d2 <= rad*rad) 
                 { assert(nv < maxnv);
                   b.ind[dz][nv] = (i2_t){{ dx,dy }};
                   b.step[dz][nv] = dy*stepy + dx*stepx;
                   nv++;
                 }
             }
           }
         if (debug) { fprintf(stderr, " %d", nv); }
         b.nv[dz] = nv;
         b.ind[dz] = realloc(b.ind[dz], nv*sizeof(i2_t));
         b.step[dz] = realloc(b.step[dz], nv*sizeof(ppv_step_t));
         nvtot += nv;
       }
     if (debug) { fprintf(stderr, "\n"); }
     if (debug) { fprintf(stderr, "total voxels in full ball = %d\n", 2*nvtot - b.nv[0]); }
     return b;
   }

void voxb_erolate_free_ball(voxb_erolate_ball_t b)
  { for(int32_t dz = 0; dz <= b.irad; dz++)
      { free(b.ind[dz]);
        free(b.step[dz]);
      }
    free(b.nv);
    free(b.ind);
    free(b.step);
  }

ppv_sample_t voxb_erolate_compute_result
  ( ppv_array_t *T, 
    ppv_size_t nzA,
    ppv_index_t kz,
    ppv_index_t ky,
    ppv_index_t kx, 
    voxb_erolate_ball_t *b,
    ppv_sample_t vexp
  )
  {
    bool_t debug = FALSE;
    
    int32_t irad = b->irad;
    
    /* Dimensions of buffer: */
    int32_t nzT = (int32_t)T->size[0];
    int32_t nyT = (int32_t)T->size[1];
    int32_t nxT = (int32_t)T->size[2];
    
    /* Assumes that {T} has the same dimensions as {A} except in {Z}: */
    int32_t nyA = nyT;
    int32_t nxA = nxT;
    
    assert(nzT >= 2*irad+1);
    ppv_index_t indT[ppv_array_NAXES];
    ppv_index_clear(indT);
    if (debug) { fprintf(stderr, "computing A[%ld,%ld,%ld] vexp = %d\n", kz,ky,kx,vexp); }
    for (int32_t dz = -irad; dz <= irad; dz++)
      { ppv_index_t rz = kz + dz; 
        if (debug) { fprintf(stderr, "  checking plane %ld of original A\n", rz); }
        if ((rz >= 0) && (rz < nzA))
          { /* Scan ball neighborhood elements of original {a} in plane {rz}. */
            int32_t dza = (dz < 0 ? -dz : dz);
            int32_t nvz = b->nv[dza];
            i2_t *indz = b->ind[dza];
            ppv_step_t *stepz = b->step[dza];
            indT[0] = rz % nzT;
            indT[1] = ky;
            indT[2] = kx;
            if (debug) { fprintf(stderr, "  assuming it is in plane %ld of T\n", indT[0]); }
            ppv_pos_t basez = ppv_sample_pos(T, indT); /* Address of voxel {[ky,kx]} on this plane of {T}. */ 
            for (int32_t i = 0; i < nvz; i++)
              { i2_t indk = indz[i];
                /* Compute {X} and {Y} indices {rx,ry} of this ball voxel in the array: */
                ppv_index_t rx = kx + indk.c[0]; 
                ppv_index_t ry = ky + indk.c[1]; 
                if ((rx >= 0) && (rx < nxA) && (ry >= 0) && (ry < nyA))
                  { /* Get the ball voxel from {T}: */
                    ppv_pos_t pos = basez + stepz[i];
                    ppv_sample_t vsmp = ppv_get_sample_at_pos(T->el,T->bps,T->bpw, pos);
                    if (debug) { fprintf(stderr, "    T[%ld,%ld,%ld] = %d\n", rz,ry,rx, vsmp); }
                    if (vsmp == vexp) { return vexp; }
                  }
              }
          }
      }
    return 1-vexp;
  }
