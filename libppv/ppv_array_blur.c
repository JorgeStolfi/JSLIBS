/* See ppv_array_blur.h */
/* Last edited on 2021-07-09 01:41:20 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <jsmath.h>
#include <ix.h>
#include <ix_io.h>
#include <sample_conv.h>
#include <affirm.h>

#include <ppv_types.h>
#include <ppv_array.h>

#include <ppv_array_blur.h>

ppv_array_t *ppv_array_blur
  ( ppv_array_t *A, 
    ppv_sample_floatize_proc_t *floatize,
    ppv_sample_t maxsmpG, 
    int32_t radius,
    double wt[],
    int32_t stride,
    ppv_sample_quantize_proc_t *quantize
  )
  {
    demand(radius >= 1, "blur radius must be positive");
    demand(maxsmpG >= 1, "invalid output max sample");
    
    demand((stride >= 1) && (stride <= radius+1), "invalid probing stride");
    demand((radius + 1) % stride == 0, "probing stride must divide {radius+1}");
    
    ppv_dim_t d = A->d;
    
    /* Compute the size of {G} and allocate it: */
    ppv_size_t szG[d];
    for (ppv_axis_t k = 0; k < d; k++) { szG[k] = (A->size[k]+stride-1)/stride; }
    ppv_array_t *G = ppv_array_new(d, szG, maxsmpG);
    
    /* Now compute the samples of {G}: */
    ppv_index_t ixA[d]; /* Index vector for {A}. */
    auto bool_t compav(const ppv_index_t *ix, ppv_pos_t pG, ppv_pos_t pB, ppv_pos_t pC);
    ppv_enum(compav, FALSE, G, NULL, NULL);
    
    return G;
    
    bool_t compav(const ppv_index_t *ix, ppv_pos_t pG, ppv_pos_t pB, ppv_pos_t pC)
      { for (ppv_axis_t k = 0; k < d; k++) { ixA[k] = (ppv_index_t)(stride*ix[k]); }
        double avg = ppv_array_blur_single(A, ixA, radius, floatize, wt);
        ppv_sample_t smpG;
        assert(! isnan(avg));
        if (quantize != NULL)
          { smpG = quantize(avg, maxsmpG); }
        else
          { smpG = (ppv_sample_t)sample_conv_quantize
              ( (float)avg, maxsmpG, FALSE, 0.0, 1.0, NULL,NULL, NULL,NULL, NULL,NULL );
          }
        ppv_set_sample_at_pos(G->el, G->bps, G->bpw, pG, smpG);
        return FALSE;
      }
  }

double ppv_array_blur_single
  ( ppv_array_t *A, 
    const ppv_index_t ix[],
    int32_t radius,
    ppv_sample_floatize_proc_t *floatize,
    double wt[]
  )
  {
    ppv_dim_t d = A->d;
    
    bool_t debug = FALSE;

    int32_t szw = 2*radius + 1; /* Neighborhood size along each axis. */
    ppv_sample_t maxsmpA = A->maxsmp; /* Max sample of {A}. */

    if (debug) { ix_print_indices (stderr, "avg around A[", d, ix, 0, ",", "]\n"); }

    double avg; /* Neighborhood average. */
    if (A->maxsmp == 0)
      { /* All samples are zero anyway: */ avg = 0.0; }
    else
      { /* Enumerate neighborhood voxels: */
        ppv_index_t jx[d];  /* Index of voxel in neighborhood */
        for (ppv_axis_t k = 0; k < d; k++) { jx[k] = ix[k] - radius; }
        double sum_wx = 0.0;  /* Weighted sum of samples in neighborhood. */
        double sum_w = 0.0;   /* Sum of weights used. */
        while (TRUE)
          { /* Compute weight of sample {A[jx]} (0 if it does not exist): */
            if (debug) { ix_print_indices (stderr, "  A[", d, jx, 0, ",", "]\n"); }
            double wj = 1.0; /* Sample weight in neighborhood. */
            for (ppv_axis_t k = 0; k < d; k++) 
              { ppv_index_t jxk = jx[k];
                if ((jxk >= 0) && (jxk < A->size[k]))
                  { /* Voxel {A[jx]} exists: */
                    ppv_index_t rix = jx[k] - ix[k] + radius; /* Index withing neighborhood. */
                    assert((rix >= 0) && (rix < szw));
                    wj = wj * wt[rix];
                  }
                else
                  { /* Voxel {A[jx]} does not exist: */
                    wj = 0.0; break;
                  }
              }
            if (wj > 0.0)
              { /* Fetch and accumulate sample {A[jx]}: */
                ppv_sample_t vj = ppv_get_sample(A, jx);
                double xj; /* Floatized sample value. */
                if (floatize != NULL)
                  { xj = floatize(vj, maxsmpA); }
                else
                  { xj = sample_conv_floatize(vj, maxsmpA, TRUE, 0.0,1.0, NULL,NULL, NULL,NULL); }
                if (debug) { fprintf(stderr, "  wj = %16.12f xj = %8.6f\n", wj, xj); }
                sum_wx += wj*xj;
                sum_w += wj;
              }
            /* Advance {jx} to the next index in neigborhood: */
            ppv_axis_t axf = 0;  /* Axis to try to increment next. */
            while((axf < d) && (jx[axf] >= ix[axf] + radius))
              { jx[axf] = ix[axf] - radius;
                axf++;
              }
            if (axf >= d) { /* Scanned the whole neighborhood */ break; }
            jx[axf]++;
          }
        /* Compute and return the average: */
        if (debug) { fprintf(stderr, " = %24.16e/%24.16e", sum_wx, sum_w); } 
        avg = (sum_w == 0.0 ? NAN : sum_wx/sum_w);
      }
    if (debug) { fprintf(stderr, " = %18.16f\n", avg); } 
    return avg;
  }
