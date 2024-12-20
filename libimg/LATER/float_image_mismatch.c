/* See {float_image_mismatch.h}. */
/* Last edited on 2024-12-04 23:21:30 by stolfi */

#include <float_image_mismatch.h>

#include <assert.h>
#include <limits.h>
#include <math.h>
 
#include <affirm.h>

/* INTERNAL PROTOTYPES */

void fimm_est_avg_var(int32_t ni, double v[], double m[], int32_t k, int32_t dk, double *Mp, double *Vp, double *Wp);
  /* Estimates the mean {*Mp} and variance {*Vp} of the samples
    {v[k+i*dk]}, for {i} in {0..ni-1}, taking into account the mask
    values {m[k+i*dk]}. Assumes that the noise variance in each {
    
    Also computes the weight {*Wp} as the
    harmonic mean of all the weights. Returns {*Mp=0} and {*Vp=0} if
    all the {w[k+i*dk]} are zero. */

/* IMPLEMENTATIONS */

double float_image_mismatch_var
  ( int32_t ni,                          /* Number of images being compared. */
    float_image_eval_sample_t *eval, /* Single-channel image evaluator. */
    int32_t hwx,                         /* Half-width of comparison window. */
    double wx[],                     /* Horizontal weight table. */
    int32_t hwy,                         /* Half-height of comparison window. */
    double wy[],                     /* Vertical weight table. */
    r2_t p[]                         /* Sampling grid center for each image. */
  )  
  { 
    int32_t nwx = 2*hwx+1;
    int32_t nwy = 2*hwy+1;
    int32_t nw = nwx*nwy; /* Number of samples in window grid. */
    
    /* Get the image samples {y[i,k] = y[i*nw + k]}, and their mask weights: */
    /* One plane with {nw} samples from each image. */
    double y[ni*nw];
    double m[ni*nw];
    int32_t i;
    for (i = 0; i < ni; i++)
      { /* Get samples of image {i}: */
        double *vi = &(v[i*nw]);
        double *wi = &(w[i*nw]);
        eval(i, &(p[i]), hwx, hwy, yi, mi);
      }
    
    /* Sums over all {k}: */
    double sum_WV2 = 0; /* Sum of mass-weighted variances. */
    double sum_W = 0;  /* Sum of masses. */

    /* Enumerate sampling points: */
    int32_t jx, jy;
    for (jx = 0; jx < nwx; jx++)
      { for (jy = 0; jy < nwy; jy++)
          { int32_t k = jx + nwx*jy;
            double M, V2, W;
            fimm_est_avg_var(ni, v, m, k, nw, &M, &V2, &W);
            W *= wx[jx]*wy[jy];
            sum_WV2 += W*V2;
            sum_W += W;
          }
      }
    return (sum_W == 0 ? 0.0 : sum_WV / sum_W);
  }

void fimm_est_avg_var(int32_t ni, double v[], double m[], int32_t k, int32_t dk, double *Mp, double *V2p, double *Wp)
  {
    /* Compute the weight : */
    double sum_wy = 0; /* Sum of {w[ki]*y[ki]} */
    double sum_w = 0;   /* Sum of {w[ki]} */
    double sum_m = 0;  /* Sum of {[ki]} */
    ki = k;
    for (i = 0; i < ni; i++)
      { double vi = v[ki]; 
        double wi = w[ki]; 
        sum_wv += wi*vi;
        sum_w += wi;
        sum_m += 1.0/wi;
        ki += dk;
      }
    double M, V2, W;
    int32_t i, ki;
    /* Compute the avg {M} and the harmonic mean weight {W}: */
    double sum_wy = 0; /* Sum of {w[ki]*y[ki]} */
    double sum_w = 0;   /* Sum of {w[ki]} */
    double sum_m = 0;  /* Sum of {[ki]} */
    ki = k;
    for (i = 0; i < ni; i++)
      { double vi = v[ki]; 
        double wi = w[ki]; 
        sum_wv += wi*vi;
        sum_w += wi;
        sum_m += 1.0/wi;
        ki += dk;
      }
    if (sum_w == 0)
      { M = V = W = 0.0; }
    else
      { M = sum_wv/sum_w; 
        W = 1.0/sum_m; 
        /* Now compute the var {V}: */
        double sum_we2 = 0;
        ki = k;
        for (i = 0; i < ni; i++)
          { double vi = v[ki]; 
            double wi = w[ki]; 
            double ei = vi - M;
            sum_we2 += wi*ei*ei;
            ki += dk;
          }
        V = sum_we2/sum_w; 
      }

    (*Mp) = M;
    (*Vp) = V;
    (*Wp) = W;
  }
    
            
    
