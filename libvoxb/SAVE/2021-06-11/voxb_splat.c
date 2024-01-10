/* See voxb_splat.h */
/* Last edited on 2021-06-10 19:21:32 by jstolfi */

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

#include <voxb_splat.h>

void voxb_splat_object_multi
  ( ppv_array_t *a,
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
        voxb_splat_object(a, obj, &(S[k]), maxR, op);
      }
    if (debug) { fprintf(stderr, "\n"); }
  }

void voxb_splat_object
  ( ppv_array_t *a, 
    r3_pred_t *obj,
    r3_motion_state_t *S,
    double maxR,
    voxb_op_t op
  )
  {
    demand(a->bps == 1, "tomogram is not binary");
    
    /* Get the voxel counts along each axis. */
    i3_t N = (i3_t){{ (int32_t)a->size[2], (int32_t)a->size[1], (int32_t)a->size[0] }}; 
    
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
                    voxb_splat_voxel(a, kx, ky, kz, val, op);
                  }
              }
          }
      }
  }
  
void voxb_splat_voxel(ppv_array_t *a, int32_t kx, int32_t ky, int32_t kz, bool_t val, voxb_op_t op)
  {
    /* Fetch the current sample {osmp}: */
    int32_t NA = ppv_array_NAXES;
    ppv_index_t ix[NA];
    int32_t j;
    for (j = 0; j < NA; j++) { ix[j] = 0; }
    ix[0] = kz; 
    ix[1] = ky; 
    ix[2] = kx; 
    ppv_pos_t pos = ppv_sample_pos(a, ix);
    ppv_sample_t osmp = ppv_get_sample_at_pos(a->el, a->bps, a->bpw, pos);
    ppv_sample_t vsmp = ((ppv_sample_t)val) & 1;
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
    if (smp != osmp) { ppv_set_sample_at_pos(a->el, a->bps, a->bpw, pos, smp); }
  }
