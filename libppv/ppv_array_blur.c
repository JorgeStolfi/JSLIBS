/* See ppv_array_blur.h */
/* Last edited on 2021-07-03 08:09:32 by jstolfi */

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
#include <sample_conv.h>
#include <affirm.h>

#include <ppv_types.h>
#include <ppv_array.h>

#include <ppv_array_blur.h>

ppv_array_t *ppv_array_blur
  ( ppv_array_t *A, 
    ppv_sample_floatize_proc_t *floatize,
    ppv_nbits_t bps, 
    ppv_size_t radius,
    double wt[],
    ppv_size_t hstep,
    ppv_sample_quantize_proc_t *quantize
  )
  {
    demand(radius >= 1, "blur radius must be positive");
    ppv_size_t t = 2*radius + 1;
    
    demand((hstep >= 0) && (hstep <= radius), "invalid hstep");
    ppv_size_t s = 2*hstep + 1;
    demand(t % s == 0, "probing step must divide window size");
    
    ppv_dim_t d = A->d;
    
    demand(A->bps > 0, "invalid input bits per sample");
    demand(bps > 0, "invalid output bits per sample");
    
    ppv_sample_t mG = (ppv_sample_t)((1LU << bps) - 1); /* Max sample of {G}. */
    
    /* Compute the size of {G} and allocate it: */
    ppv_size_t szG[d];
    for (ppv_axis_t k = 0; k < d; k++) { szG[k] = A->size[k]/s; }
    ppv_nbits_t bpw = ppv_best_bpw(bps);
    ppv_array_t *G = ppv_array_new(d, szG, bps, bpw);
    
    /* Now compute the samples of {G}: */
    ppv_index_t ixA[d]; /* Index vector for {A}. */
    auto bool_t compav(const ppv_index_t *ix, ppv_pos_t pG, ppv_pos_t pB, ppv_pos_t pC);
    ppv_enum(compav, FALSE, G, NULL, NULL);
    
    return G;
    
    bool_t compav(const ppv_index_t *ix, ppv_pos_t pG, ppv_pos_t pB, ppv_pos_t pC)
      { for (ppv_axis_t k = 0; k < d; k++) { ixA[k] = s*ix[k] + hstep; }
        double avg = ppv_array_blur_single(A, ixA, radius, floatize, wt);
        ppv_sample_t vG;
        assert(! isnan(avg));
        if (quantize != NULL)
          { vG = quantize(avg, mG); }
        else
          { vG = (ppv_sample_t)sample_conv_quantize
              ( (float)avg, mG, FALSE, 0.0, 1.0, NULL,NULL, NULL,NULL, NULL,NULL );
          }
        ppv_set_sample_at_pos(G->el, bps, bpw, pG, vG);
        return FALSE;
      }
  }

double ppv_array_blur_single
  ( ppv_array_t *A, 
    const ppv_index_t ix[],
    ppv_size_t radius,
    ppv_sample_floatize_proc_t *floatize,
    double wt[]
  )
  {
    bool_t debug = FALSE;

    ppv_dim_t d = A->d;
    ppv_size_t t = 2*radius + 1; /* Neighborhood size along each axis. */
    ppv_sample_t mA = (ppv_sample_t)((1LU << A->bps) - 1); /* Max sample of {A}. */

    if (debug) { ix_print_indices (stderr, "avg arount A[", d, ix, 0, ",", "]"); }

    double avg; /* Neighborhood average. */
    if (A->bps == 0)
      { /* All samples are zero anyway: */ avg = 0.0; }
    else
      { /* Enumerate neighborhood voxels: */
        ppv_index_t jx[d];  /* Index of voxel in neighborhood */
        for (ppv_axis_t k = 0; k < d; k++) { jx[k] = ix[k] - radius; }
        double sum_wx = 0.0;  /* Weighted sum of samples in neighborhood. */
        double sum_w = 0.0;   /* Sum of weights used. */
        while (TRUE)
          { /* Compute weight of sample {A[jx]} (0 if it does not exist): */
            double wj = 1.0; /* Sample weight in neighborhood. */
            for (ppv_axis_t k = 0; k < d; k++) 
              { ppv_index_t jxk = jx[k];
                if ((jxk >= 0) && (jxk < A->size[k]))
                  { /* Voxel {A[jx]} exists: */
                    ppv_index_t rix = jx[k] - ix[k] + radius; /* Index withing neighborhood. */
                    assert((rix >= 0) && (rix < t));
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
                  { xj = floatize(vj, mA); }
                else
                  { xj = sample_conv_floatize(vj, mA, TRUE, 0.0,1.0, NULL,NULL, NULL,NULL); }
                sum_wx += wj*xj;
                sum_w += wj;
              }
            /* Advance {jx} to the next index in neigborhood: */
            ppv_axis_t axf = 0;  /* Axis to try to increment next. */
            while((axf < d) && (jx[axf] >= ix[axf] + radius))
              { jx[axf] = ix[axf] + radius;
                axf++;
              }
            if (axf >= d) { /* Scanned the whole neighborhood */ break; }
            jx[axf]++;
          }
        /* Compute and return the average: */
        if (debug) { fprintf(stderr, " = %24.16e/%24.1e", sum_wx, sum_w); } 
        avg = (sum_w == 0.0 ? NAN : sum_wx/sum_w);
      }
    if (debug) { fprintf(stderr, " = %18.16f\n", avg); } 
    return avg;
  }
