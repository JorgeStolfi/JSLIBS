/* See voxm_splat.h */
/* Last edited on 2021-06-12 12:15:50 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r3.h>
#include <r3_motion.h>
#include <r3_extra.h>
#include <r3x3.h>
#include <i3.h>
#include <affirm.h>
#include <ppv_array.h>
#include <vec.h>

#include <voxm_splat.h>

void voxm_splat_object_multi
  ( ppv_array_desc_t *A,
    r3_double_func_t *obj,
    int32_t ns,
    r3_motion_state_t S[],
    double maxR,
    bool_t sub
  )
  {
    bool_t debug = TRUE;
    int32_t k;
    for (k = 0; k < ns; k++) 
      { if (debug) { fprintf(stderr, "."); }
        voxm_splat_object(A, obj, &(S[k]), maxR, sub);
      }
    if (debug) { fprintf(stderr, "\n"); }
  }

void voxm_splat_object
  ( ppv_array_desc_t *A, 
    r3_double_func_t *obj,
    r3_motion_state_t *S,
    double maxR,
    bool_t sub
  )
  {
    
    /* Get the cell counts along each axis. */
    i3_t N = (i3_t){{ (int32_t)A->size[2], (int32_t)A->size[1], (int32_t)A->size[0] }}; 
    
    /* Compute the bounding box of the brush: */
    i3_t kmin;
    i3_t kmax;
    int32_t j;
    for (j = 0; j < 3; j++)
      { kmin.c[j] = (int32_t)floor(S->p.c[j] - maxR - 1.0);
        if (kmin.c[j] < 0) { kmin.c[j] = 0; }
        kmax.c[j] = (int32_t)ceil (S->p.c[j] + maxR + 1.0);
        if (kmax.c[j] >= N.c[j]) { kmax.c[j] = N.c[j] - 1; }
      }
      
    /* Get inverse of pose matrix: */
    r3x3_t Minv;
    r3x3_inv(&(S->M), &Minv);
    
    /* Enumerate voxels in bounding box: */
    int32_t kx, ky, kz;
    for (kz = kmin.c[2]; kz <= kmax.c[2]; kz++)
      { for (ky = kmin.c[1]; ky <= kmax.c[1]; ky++)
          { for (kx = kmin.c[0]; kx <= kmax.c[0]; kx++)
              { /* Get coordinates of pixel center: */
                r3_t qvox = (r3_t){{ kx + 0.5, ky + 0.5, kz + 0.5 }};
                /* Map to reference object coordinates: */
                r3_t qobj;
                r3_sub(&qvox, &(S->p), &qobj);
                r3x3_map_row(&qobj, &Minv, &qobj);
                /* Evaluate object there: */
                double val = obj(&qobj);
                if (val > 0.0)
                  { /* Modify voxel value according to coverage: */
                    voxm_splat_voxel(A, kx, ky, kz, val, sub);
                  }
              }
          }
      }
  }
  
void voxm_splat_voxel(ppv_array_desc_t *A, int32_t kx, int32_t ky, int32_t kz, double val, bool_t sub)
  {
    /* Quantize the value: */
    ppv_sample_t maxsmp = (ppv_sample_t)((1u << A->bps) - 1); /* Max sample value. */
    ppv_sample_t smp; /* Value {val} quantized to {0..maxsmp}. */
    if (val <= 0.0)
      { smp = 0; }
    else if (val >= 1.0)
      { smp = maxsmp; }
    else
      { smp = (ppv_sample_t)floor(val * maxsmp + 0.4999999); }
      
    /* Fetch the current sample {osmp}: */
    ppv_dim_t d = A->d;
    assert(d == 3);
    ppv_index_t ix[d];
    ix[0] = kz; 
    ix[1] = ky; 
    ix[2] = kx; 
    ppv_pos_t pos = ppv_sample_pos(A, ix);
    ppv_sample_t osmp = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pos);
            
    if (!sub)
      { /* Set voxel value to maximum of the two: */
        if (smp > osmp)
          { /* Update value in array: */
            ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, smp);
          }
      }
    else
      { /* Set voxel value to min of old value and complement of new value: */
        smp = (ppv_sample_t)(maxsmp - smp);
        if (smp < osmp)
          { /* Update value in array: */
            ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, smp);
          }
      }
  }
