/* See voxb_splat.h */
/* Last edited on 2024-11-11 07:38:26 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r3.h>
#include <r3_motion.h>
#include <r3x3.h>
#include <i3.h>
#include <affirm.h>
#include <ppv_array.h>
#include <vec.h>

#include <voxb_splat.h>

void voxb_splat_object_multi
  ( ppv_array_t *A,
    r3_pred_t *obj,
    int32_t ns,
    r3_motion_state_t S[],
    double maxR,
    voxb_op_t op
  )
  {
    bool_t debug = TRUE;
    int32_t k;
    for (k = 0; k < ns; k++) 
      { if (debug) { fprintf(stderr, "."); }
        voxb_splat_object(A, obj, &(S[k]), maxR, op, FALSE);
      }
    if (debug) { fprintf(stderr, "\n"); }
  }

void voxb_splat_object
  ( ppv_array_t *A, 
    r3_pred_t *obj,
    r3_motion_state_t *S,
    double maxR,
    voxb_op_t op,
    bool_t debug
  )
  {
    demand(A->bps == 1, "tomogram is not binary");
    
    /* Get the voxel counts along each axis. */
    i3_t N = (i3_t){{ (int32_t)A->size[2], (int32_t)A->size[1], (int32_t)A->size[0] }}; 
    
    /* Compute the bounding box {kmin..kmax} of the affected voxels: */
    i3_t kmin, kmax; /* Elements are {X,Y,Z}. */
    for (uint32_t j = 0;  j < 3; j++)
      { if (op == voxb_op_AND)
          { /* Can affect anywhere: */
            fprintf(stderr, "<--> op = %d\n", op);
            kmin.c[j] = 0; kmax.c[j] = N.c[j] - 1;
          }
        else
          { /* Will affect only where there is object: */
            kmin.c[j] = (int32_t)floor(S->p.c[j] - maxR - 1.0);
            if (kmin.c[j] < 0) { kmin.c[j] = 0; }
            kmax.c[j] = (int32_t)ceil (S->p.c[j] + maxR + 1.0);
            if (kmax.c[j] >= N.c[j]) { kmax.c[j] = N.c[j] - 1; }
          }
      }
      
    if (debug)
       { fprintf(stderr, "splatting object with op = %d maxR = %.3f coord ranges = ", op, maxR);
         for (uint32_t j = 0;  j < 3; j++)
           { fprintf(stderr, " %d..%d (%d)", kmin.c[j], kmax.c[j], kmax.c[j] - kmin.c[j] + 1); }
         fprintf(stderr, "\n");
       }
    
    /* Get inverse of pose matrix: */
    r3x3_t Minv;
    r3x3_inv(&(S->M), &Minv);
    
    bool_t null_val; /* Null value for {op}. */
    if ((op == voxb_op_OR) || (op == voxb_op_SUB) || (op == voxb_op_XOR))
      { null_val = FALSE; }
    else if (op == voxb_op_AND)
      { null_val = TRUE; }
    else
      { demand(FALSE, "invalid {op}"); }
    
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
                bool_t val = obj(&qobj);
                if (val != null_val)
                  { /* Modify voxel value: */
                    voxb_splat_voxel(A, kx, ky, kz, val, op, debug);
                  }
              }
          }
      }
  }

void voxb_splat_voxel
  ( ppv_array_t *A, 
    int32_t kx, 
    int32_t ky, 
    int32_t kz, 
    bool_t val, 
    voxb_op_t op, 
    bool_t debug
  )
  {
    if (debug) 
      { fprintf(stderr, "  {voxb_splat_voxel}: ix = [ %d %d %d ]\n", kz,ky,kx); }
    
    /* Fetch the current sample {osmp}: */
    ppv_dim_t d = A->d;
    ppv_index_t ix[d];
    for (ppv_axis_t j = 0; j < d; j++) { ix[j] = 0; }
    ix[0] = kz; 
    ix[1] = ky; 
    ix[2] = kx; 
    ppv_pos_t pos = ppv_sample_pos(A, ix);
    ppv_sample_t osmp = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pos);
    ppv_sample_t vsmp = (ppv_sample_t)(val & 1);
    ppv_sample_t smp; /* New sample value. */
    switch (op)
      { 
        case voxb_op_OR:   smp = osmp | vsmp; break;
        case voxb_op_AND:  smp = osmp & vsmp; break;
        case voxb_op_SUB:  smp = osmp & (~vsmp); break;
        case voxb_op_XOR:  smp = osmp + vsmp; break;
        default: assert(FALSE);
      }
    smp = smp & 1; /* Get rid of any spurious bits. */   
    if (smp != osmp) { ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, smp); }
  }

