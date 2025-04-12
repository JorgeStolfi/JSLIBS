/* See {r2_align_multiscale.h}. */
/* Last edited on 2025-03-19 12:46:16 by stolfi */

#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <math.h>
 
#include <bool.h>
#include <sign.h>
#include <r2.h>
#include <i2.h>
#include <jsmath.h>
#include <affirm.h>
#include <sve_minn.h>
#include <sve_minn_iterate.h>

#include <r2_align.h>
#include <r2_align_enum.h>
#include <r2_align_quadopt.h>

#include <r2_align_multiscale.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

void r2_align_multiscale
  ( int32_t ni,                          /* Number of images to align. */
    r2_align_multiscale_mismatch_t *F2,  /* Function that evaluates the mismatch between the images. */
    bool_t quadopt,                      /* Use quadratic optimization? */
    r2_t arad[],                         /* Max delta vector coordinates along each axis. */
    bool_t bal,                          /* True if alignment vector adjustments should be balanced. */
    double tol,                          /* Desired precision. */
    r2_t p[],                            /* (IN/OUT) Corresponding points in each image. */
    double *F2val_P                      /* (OUT) Mismatch for the computed alignment vector. */
  )
  {
    /* Find the maximum search radius: */
    double rmax = 0;
    for (uint32_t i = 0;  i < ni; i++)
      { for (uint32_t j = 0;  j < 2; j++)
          { double rij = arad[i].c[j];
            demand(rij >= 0, "invalid search radius");
            if (rij > rmax) { rmax = rij; }
          }
      }
      
    if (rmax == 0)
      { /* Nothing to do: */ 
        i2_t iscale0 = (i2_t){{ 0, 0 }};
        (*F2val_P) = F2(ni, p, iscale0);
      }
    else
      { 
        /* Compute the initial scale {smax} and its reduction factor {fscale}: */
        int32_t smax = 0;
        double fscale = 1.0;
        while (rmax > 0.5) { smax = smax+1; rmax = rmax/2; fscale = fscale/2; }
        /* Reduce the problem to scale {smax}: */
        r2_t srad[ni];  /* Search radius at current scale. */
        for (uint32_t i = 0;  i < ni; i++)
          { for (uint32_t j = 0;  j < 2; j++) 
              { srad[i].c[j] = arad[i].c[j]*fscale;
                p[i].c[j] = p[i].c[j]*fscale;
              }
          }
        /* Now solve the problem at increasing scales: */
        int32_t scale = smax;
        while(TRUE)
          { /* Solve the problem at scale {scale}: */
            i2_t iscale = (i2_t){{ scale, scale }};
            if (quadopt)
              { r2_align_multiscale_single_scale_quadopt(ni, iscale, F2, srad, bal, tol, p, F2val_P); }
            else
              { r2_align_multiscale_single_scale_enum(ni, iscale, F2, srad, bal, tol, p, F2val_P); }
            /* Are we done? */
            if (scale == 0) { break; }
            /* Expand to the next finer scale: */
            for (uint32_t i = 0;  i < ni; i++)
              { for (uint32_t j = 0;  j < 2; j++) 
                  { srad[i].c[j] = fmin(0.5, 2 * srad[i].c[j]);
                    p[i].c[j] = 2 * p[i].c[j];
                  }
              }
          
          }
      }
  }

void r2_align_multiscale_single_scale_enum
  ( int32_t ni,                         /* Number of objects to align. */
    i2_t iscale,                        /* Object scaling exponent along each axis. */  
    r2_align_multiscale_mismatch_t *F2, /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],                        /* Max delta vector coordinates for each object. */
    bool_t bal,                         /* True if alignment vector adjustments should be balanced. */
    double tol,                         /* Desired precision. */
    r2_t p[],                           /* (IN/OUT) Corresponding points in each object. */
    double *F2val_P                     /* (OUT) Mismatch for the computed alignment vector. */
  )
  {
    auto double single_F2(int32_t ni, r2_t p[]);
    
    r2_align_enum(ni, single_F2, arad, bal, tol, p, F2val_P);
    return;
    
    double single_F2(int32_t ni, r2_t p[])
      { return F2(ni, p, iscale); }
  }
  
void r2_align_multiscale_single_scale_quadopt
  ( int32_t ni,                         /* Number of objects to align. */
    i2_t iscale,                        /* Object scaling exponent along each axis. */  
    r2_align_multiscale_mismatch_t *F2, /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],                        /* Max delta vector coordinates for each object. */
    bool_t bal,                         /* True if alignment vector adjustments should be balanced. */
    double tol,                         /* Desired precision. */
    r2_t p[],                           /* (IN/OUT) Corresponding points in each object. */
    double *F2val_P                     /* (OUT) Mismatch for the computed alignment vector. */
  )
  {
    auto double single_F2(int32_t ni, r2_t p[]);
    
    r2_align_quadopt(ni, single_F2, arad, bal, tol, p, F2val_P);
    return;
    
    double single_F2(int32_t ni, r2_t p[])
      { return F2(ni, p, iscale); }
  }
