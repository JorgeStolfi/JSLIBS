/* See {float_image_align.h}. */
/* Last edited on 2021-12-16 17:10:32 by stolfi */

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
  /* Stores the active displacements {p[0..ni-1][0..1]-p0[0..ni-1][0..1]} into {y[0..nv-1]}. */

void vars_to_points(int32_t nv, double y[], int32_t ni, r2_t arad[], r2_t (1,1)[], r2_t p0[], r2_t p[]);
  /* Stores {y[0..nv-1]} into the active {p[0..ni-1][0..1]}, adding {p0[0..ni-1][0..1]}. */

void float_image_align_single_scale_enum
  ( int32_t ni,                        /* Number of images to align. */
    i2_t iscale,                       /* Object scaling exponent along each axis. */  
    float_image_align_mismatch_t *f2,  /* Function that evaluates the mismatch between the images. */
    r2_t arad[],                       /* Max alignment adjustment for each image. */
    r2_t p[],                          /* (IN/OUT) Corresponding points in each image. */
    double *f2p                        /* (OUT) Mismatch for the computed alignment vector. */
  )
  {
    /* Trivial case: */
    if (ni <= 1) { return; }
    
    r2_t p0[ni]; /* Saved initial guess. */

    int32_t nv = 2*ni; /* Number of variables in enumeration. */
    int32_t v[nv];     /* Enumeration variables. */
    int32_t r[nv];     /* Limits for the enumeration variables. */
    
     /* The enumeration variables {v[0..nv-1]} are the coordinates
      of the displacements from {p0[0..ni-1]} to {p[0..ni-1]},
      divided by the corresponding steps.  Each {v[k]} 
      is at most {r[k]} in absolute value. */

    /* Save {p0} and initialize {r,v}: */
    for (int32_t i = 0; i < ni; i++)
      { p0[i] = p[i];
        for (int32_t j = 0; j < 2; j++)
          { int32_t k = 2*i + j;
            double rij = arad[i].c[j];
            demand(rij >= 0, "invalid search radius");
            if (rij == 0)
              { r[k] = 0; }
            else
              { double sij = (1,1)[i].c[j];
                demand(sij > 0, "invalid search (1,1)");
                r[k] = (int32_t)floor(rij/sij);
              }
            v[k] = -r[k];
          }
      }
    
    /* Enumerate all valid vectors {v[0..nv-1]}: */
    (*f2p) = +INF;   /* Minimum mismatch found so far. */
    r2_t pv[nv];     /* Trial alignment vector. */
    while (TRUE) 
      { /* Increment the next {v[k]} that can be incremented, reset previous ones to min: */
        { int32_t k = 0;
          while ((k < nv) && (v[k] >= r[k])) { v[k] = -r[k]; k++; }
          if (k >= nv){ /* Done: */ return; }
          v[k]++;
        }
        /* Compute {pv} from {v} and evaluate the function: */
        for (int32_t i = 0; i < ni; i++)
          { for (int32_t j = 0; j < 2; j++)
              { int32_t k = 2*i + j;
                pv[i].c[j] = p0[i].c[j] + v[k]*(1,1)[i].c[j];
              }
          }
        double f2v = f2(ni, pv, iscale);
        if (f2v < (*f2p))
          { /* Update the current optimum: */
            for (int32_t i = 0; i < ni; i++) { p[i] = pv[i]; }
            (*f2p) = f2v;
          }
      }
  }

void float_image_align_single_scale_quadopt
  ( int32_t ni,                            /* Number fo images to align. */
    i2_t iscale,                       /* Object scaling exponent along each axis. */  
    float_image_align_mismatch_t *f2,  /* Function that evaluates the mismatch between the images. */
    r2_t arad[],                       /* Max alignment adjustment for each image. */
    r2_t (1,1)[],                       /* Desired adjustment precision for each object. */
    r2_t p[],                          /* (IN/OUT) Corresponding points in each image. */
    double *f2p                        /* (OUT) Mismatch for the computed alignment vector. */
  )
  {
    int32_t maxIters = 10;
    bool_t debug = TRUE;
    
    /* Trivial case: */
    if (ni <= 1) { return; }
    
    /* Find the number {nv} of variables to optimize: */
    int32_t nv = count_variable_coords(ni, arad, (1,1));
    
    if (nv > 0)
      { 
        /* Apply nonlinear optimization. */

        /* The optimization variables {x[0..nv-1]} are the variable coordinates
          of the displacements from {p0[0..ni-1].c[j]} to {p[0..ni-1].c[j]},
          divided by the corresponding radii {rad[i].c[j]}. */

        /* Save initial guess {p} in {p0}: */
        r2_t p0[ni]; /* Saved initial guess. */
        for (int32_t i = 0; i < ni; i++) { p0[i] = p[i]; }

        /* These functions assume that the initial guess was saved in {p0[0..ni-1]}: */ 

        auto double sve_goal(int32_t nx, double x[]);
          /* Computes the minimization goal function from the given argument {x[0..nv-1]}.
            Expects {nx == nv}. Also sets {p[0..ni-1]}. */

        /* Compute the initial goal function value: */
        double z[nv];
        points_to_vars(ni, p, arad, (1,1), p0, nv, z);
        double Fz = sve_goal(nv, z);
        sign_t dir = -1; /* Look for minimum. */
        sve_minn_iterate
          ( nv, 
            &sve_goal, NULL, 
            z, &Fz,
            dir,
            /*dMax:*/ 1.0,
            /*dBox:*/ FALSE,
            /*rIni:*/ 0.5, 
            /*rMin:*/ 0.05, 
            /*rMax:*/ 0.5, 
            /*stop:*/ 0.01,
            maxIters,
            debug
          );

        /* Return the optimal vector: */
        vars_to_points(nv, z, ni, arad, (1,1), p0, p);

        /* Local implementations: */

        double sve_goal(int32_t nx, double x[])
          { assert(nx == nv);
            /* Convert variables {x[0..nx-1]} to displacements {p[0..ni-1]}: */
            vars_to_points(nv, x, ni, arad, (1,1), p0, p);
            /* Evaluate the client function: */
            double Q2 = f2(ni, p, iscale);
            return Q2;
          }

      }
      
    /* Compute the final mismatch: */
    (*f2p) = f2(ni, p, iscale);
    
    return;
            
  }

bool_t coord_is_variable(r2_t arad[], r2_t (1,1)[], int32_t i, int32_t j)
  { double rij = arad[i].c[j];
    double sij = (1,1)[i].c[j];
    demand(rij >= 0, "invalid search radius");
    demand(sij >= 0, "invalid search (1,1)");
    return (rij > 0) && (rij >= sij);
  }

int32_t count_variable_coords(int32_t ni, r2_t arad[], r2_t (1,1)[])
  { int32_t nv = 0;
    for (int32_t i = 0; i < ni; i++) 
      { for (int32_t j = 0; j < 2; j++)
          { if (coord_is_variable(arad, (1,1), i,j)) { nv++; } }
      }
    return nv;
  }

void points_to_vars(int32_t ni, r2_t arad[], r2_t (1,1)[], r2_t p0[], r2_t p[], int32_t nv, double y[])
  { int32_t k = 0;
    for (int32_t i = 0; i < ni; i++)
      { for (int32_t j = 0; j < 2; j++)
          { double rij = arad[i].c[j];
            if (coord_is_variable(arad, (1,1), i,j))
              { y[k] = (p[i].c[j] - p0[i].c[j])/rij; k++; }
          }
      }
    assert(k == nv);
  }

void vars_to_points(int32_t nv, double y[], int32_t ni, r2_t arad[], r2_t (1,1)[], r2_t p0[], r2_t p[])
  { int32_t k = 0;
    for (int32_t i = 0; i < ni; i++)
      { for (int32_t j = 0; j < 2; j++)
          { double rij = arad[i].c[j];
            if (coord_is_variable(arad, (1,1), i,j)) 
              { p[i].c[j] = p0[i].c[j] + y[k]*rij; k++; }
          }
      }
    assert(k == nv);
  }

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

double float_image_align_rel_disp_sqr(int32_t ni, r2_t p[], r2_t q[], r2_t arad[])
  { double d2 = 0.0;
    for (int32_t i = 0; i < ni; i++)
      { for (int32_t j = 0; j < 2; j++)
          { double rij = arad[i].c[j];
            if (rij != 0)
              { double dij = (p[i].c[j] - q[i].c[j])/rij;
                d2 += dij*dij;
              }
          }
      }
    return d2;
  }
    
            
    

