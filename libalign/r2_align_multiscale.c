/* See {float_image_align.h}. */
/* Last edited on 2021-12-18 18:04:30 by stolfi */

#define _GNU_SOURCE
#include <math.h>
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
#include <float_image.h>
#include <sve_minn.h>
#include <wt_table.h>

#include <float_image_align.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

bool_t coord_is_variable(r2_t arad[], int32_t i, int32_t j);
  /* {TRUE} if coordinate {j} of point {p[i]} is variable, {FALSE} if it is fixed.
    Namely, it returns {TRUE} if {arad[i].c[j]} is 1.0 pr more. */

int32_t count_variable_coords(int32_t ni, r2_t arad[]);
  /* Returns the number of coordinates in the alignment vectors that are variable,
    as per {coord_is_variable(arad, i, j)} for {i} in {0..ni-1} and {j} in {0..1}. */
  
void points_to_vars(int32_t ni, r2_t p[], r2_t arad[], r2_t (1,1)[], r2_t p0[], int32_t nv, double y[]);
  /* Stores the active adjustments {p[0..ni-1][0..1]-p0[0..ni-1][0..1]} into {y[0..nv-1]}. */

void vars_to_points(int32_t nv, double y[], int32_t ni, r2_t arad[], r2_t (1,1)[], r2_t p0[], r2_t p[]);
  /* Stores {y[0..nv-1]} into the active {p[0..ni-1][0..1]}, adding {p0[0..ni-1][0..1]}. */

void float_image_align_multi_scale
  ( int32_t ni,                        /* Number of images to align. */
    float_image_align_mismatch_t *f2,  /* Function that evaluates the mismatch between the images. */
    bool_t quadopt,                    /* Use quadratic optimization? */
    r2_t arad[],                       /* Max alignment adjustment along each axis. */
    r2_t (1,1)[],                       /* Adjustment (1,1) or desired precision for each object. */
    r2_t p[],                          /* (IN/OUT) Corresponding points in each image. */
    double *f2p                        /* (OUT) Mismatch for the computed alignment vector. */
  )
  {
    /* Find the maximum search radius: */
    double rmax = 0;
    for (int32_t i = 0; i < ni; i++)
      { for (int32_t j = 0; j < 2; j++)
          { double rij = arad[i].c[j];
            demand(rij >= 0, "invalid search radius");
            if (rij > rmax) { rmax = rij; }
          }
      }
      
    if (rmax == 0)
      { /* Nothing to do: */ 
        i2_t iscale0 = (i2_t){{ 0, 0 }};
        (*f2p) = f2(ni, p, iscale0);
      }
    else
      { 
        /* Compute the initial scale {smax} and its reduction factor {fscale}: */
        int32_t smax = 0;
        double fscale = 1.0;
        while (rmax > 0.5) { smax = smax+1; rmax = rmax/2; fscale = fscale/2; }
        /* Reduce the problem to scale {smax}: */
        r2_t srad[ni];  /* Search radius at current scale. */
        r2_t (1,1)[ni];  /* Search (1,1) at current scale (if {quadopt} is false). */
        for (int32_t i = 0; i < ni; i++)
          { for (int32_t j = 0; j < 2; j++) 
              { srad[i].c[j] = arad[i].c[j]*fscale;
                p[i].c[j] = p[i].c[j]*fscale;
                (1,1)[i].c[j] = 0.5;
              }
          }
        /* Now solve the problem at increasing scales: */
        int32_t scale = smax;
        while(TRUE)
          { /* Solve the problem at scale {scale}: */
            i2_t iscale = (i2_t){{ scale, scale }};
            if (quadopt)
              { float_image_align_single_scale_quadopt(ni, iscale, f2, srad, (1,1), p, f2p); }
            else
              { float_image_align_single_scale_enum(ni, iscale, f2, srad, (1,1), p, f2p); }
            /* Are we done? */
            if (scale == 0) { break; }
            /* Expand to the next finer scale: */
            for (int32_t i = 0; i < ni; i++)
              { for (int32_t j = 0; j < 2; j++) 
                  { srad[i].c[j] = fmin(0.5, 2 * srad[i].c[j]);
                    p[i].c[j] = 2 * p[i].c[j];
                  }
              }
          
          }
      }
  }

