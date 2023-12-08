/* See {float_image_align.h}. */
/* Last edited on 2023-11-25 18:20:25 by stolfi */

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

#include <float_image_align.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

??void float_image_align_single_scale_enum
  ( int ni,                   /* Number of objects to align. */
    i2_t iscale,              /* Object scaling exponent along each axis. */  
    float_image_align_mismatch_t *f2,  /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],              /* Max alignment adjustment for each object. */
    r2_t step[],              /* Adjustment step for each object. */
    r2_t p[],                 /* (IN/OUT) Corresponding points in each object. */
    double *f2p               /* (OUT) Mismatch for the computed alignment vector. */
  )

void float_image_align_single_scale_enum
  ( int ni,                            /* Number fo images to align. */
    int scale,                         /* Image scale. */  
    float_image_align_mismatch_t *f2,  /* Function that evaluates the mismatch between the images. */
    r2_t rad[],                        /* Max alignment adjustment for each image. */
    r2_t step[],                       /* Granularity of alignment in each coordinate. */
    r2_t p[],                           /* (IN/OUT) Corresponding points in each image. */
    double *f2p                        /* (OUT) Mismatch for the computed alignment vector. */
  )
  {
    /* Trivial case: */
    if (ni <= 1) { return; }
    
    r2_t p0[ni]; /* Saved initial guess. */

    int nv = 2*ni; /* Number of variables in enumeration. */
    int v[nv];     /* Enumeration variables. */
    int r[nv];     /* Limits for the enumeration variables. */
    
     /* The enumeration variables {v[0..nv-1]} are the coordinates
      of the displacements from {p0[0..ni-1]} to {p[0..ni-1]},
      divided by the corresponding steps.  Each {v[k]} 
      is at most {r[k]} in absolute value. */

    /* Save {p0} and initialize {r,v}: */
    int i, j;
    for (i = 0; i < ni; i++)
      { p0[i] = p[i];
        for (j = 0; j < 2; j++)
          { int k = 2*i + j;
            double rij = rad[i].c[j];
            demand(rij >= 0, "invalid search radius");
            if (rij == 0)
              { r[k] = 0; }
            else
              { double sij = step[i].c[j];
                demand(sij > 0, "invalid search step");
                r[k] = (int)floor(rij/sij);
              }
            v[k] = -r[k];
          }
      }
    
    /* Enumerate all valid vectors {v[0..nv-1]}: */
    (*f2p) = +INF;   /* Minimum mismatch found so far. */
    r2_t pv[nv];     /* Trial alignment vector. */
    while (TRUE) 
      { /* Increment the next {v[k]} that can be incremented, reset previous ones to min: */
        { int k = 0;
          while ((k < nv) && (v[k] >= r[k])) { v[k] = -r[k]; k++; }
          if (k >= nv){ /* Done: */ return; }
          v[k]++;
        }
        /* Compute {pv} from {v} and evaluate the function: */
        for (i = 0; i < ni; i++)
          { for (j = 0; j < 2; j++)
              { int k = 2*i + j;
                pv[i].c[j] = p0[i].c[j] + v[k]*step[i].c[j];
              }
          }
        double f2v = f2(ni,scale,pv);
        if (f2v < (*f2p))
          { /* Update the current optimum: */
            for (i = 0; i < ni; i++) { p[i] = pv[i]; }
            (*f2p) = f2v;
          }
      }
  }

void float_image_align_single_scale_quadopt
  ( int ni,                            /* Number fo images to align. */
    int scale,                         /* Image scale. */  
    float_image_align_mismatch_t *f2,  /* Function that evaluates the mismatch between the images. */
    r2_t rad[],                        /* Max alignment adjustment for each image. */
    r2_t p[],                          /* (IN/OUT) Corresponding points in each image. */
    double *f2p                        /* (OUT) Mismatch for the computed alignment vector. */
  )
  {
    int maxIters = 10;
    bool_t debug = TRUE;
    
    /* Trivial case: */
    if (ni <= 1) { return; }
    
    /* Find the number {nv} of variables to optimize: */
    int nv = 0;
    { int i, j; 
      for (i = 0; i < ni; i++) 
        { for (j = 0; j < 2; j++)
            { double rij = rad[i].c[j];
              demand(rij >= 0, "invalid search radius");
              if (rij > 0) { nv++; }
            }
        }
    }
    
    if (nv > 0)
      { 
        /* Apply nonlinear optimization. */

        /* The optimization variables {x[0..nv-1]} are the coordinates
          of the displacements from {p0[0..ni-1].c[j]} to {p[0..ni-1].c[j]},
          divided by the corresponding radii {rad[i].c[j]}. */

        /* Save initial guess {p} in {p0}: */
        r2_t p0[ni]; /* Saved initial guess. */
        { int i; for (i = 0; i < ni; i++) { p0[i] = p[i]; } }

        /* These functions assume that the initial guess was saved in {p0[0..ni-1]}: */ 

        auto void points_to_vars(r2_t q[], double y[]);
          /* Stores the active displacements {q[0..ni-1]-p0[0..ni-1]} into {y[0..nv-1]}. */

        auto void vars_to_points(double y[], r2_t q[]);
          /* Stores {y[0..nv-1]} into the active {q[0..ni-1]}, adding {p0[0..ni-1]}. */

        auto double sve_goal(int nx, double x[]);
          /* Computes the minimization goal function from the given argument {x[0..nv-1]}.
            Expects {nx == nv}. Also sets {p[0..ni-1]}. */

        /* Compute the initial goal function value: */
        double z[nv];
        points_to_vars(p, z);
        double Fz = sve_goal(nv, z);
        sign_t dir = -1; /* Look for minimum. */
        sve_minn_iterate
          ( nv, 
            &sve_goal, NULL, 
            z, &Fz,
            dir,
            /*dMax:*/ 1.0,
            /*rIni:*/ 0.5, 
            /*rMin:*/ 0.05, 
            /*rMax:*/ 0.5, 
            /*stop:*/ 0.01,
            maxIters,
            debug
          );

        /* Return the optimal vector: */
        vars_to_points(z, p);

        /* Local implementations: */

        void points_to_vars(r2_t q[], double y[])
          { int i, j;
            int k = 0;
            for (i = 0; i < ni; i++)
              { for (j = 0; j < 2; j++)
                  { double rij = rad[i].c[j];
                    if (rij > 0) { y[k] = (q[i].c[j] - p0[i].c[j])/rij; k++; }
                  }
              }
            assert(k == nv);
          }

        void vars_to_points(double y[], r2_t q[])
          { int i, j;
            int k = 0;
            for (i = 0; i < ni; i++)
              { for (j = 0; j < 2; j++)
                  { double rij = rad[i].c[j];
                    if (rij > 0) { q[i].c[j] = p0[i].c[j] + y[k]*rij; k++; }
                  }
              }
            assert(k == nv);
          }

        double sve_goal(int nx, double x[])
          { assert(nx == nv);
            /* Convert variables {x[0..nx-1]} to displacements {p[0..ni-1]}: */
            vars_to_points(x, p);
            /* Evaluate the client function: */
            double Q2 = f2(ni, scale, p);
            return Q2;
          }

      }
      
    /* Compue the final mismatch: */
    (*f2p) = f2(ni,scale,p);
    
    return;
  }

void float_image_align_multi_scale
  ( int ni,                            /* Number of images to align. */
    float_image_align_mismatch_t *f2,  /* Function that evaluates the mismatch between the images. */
    bool_t quadopt,                    /* Use quadratic optimization? */
    r2_t rad[],                        /* Max alignment adjustment along each axis. */
    r2_t p[],                          /* (IN/OUT) Corresponding points in each image. */
    double *f2p                        /* (OUT) Mismatch for the computed alignment vector. */
  )
  {
    /* Find the maximum search radius: */
    int i, j;
    double rmax = 0;
    for (i = 0; i < ni; i++)
      { for (j = 0; j < 2; j++)
          { double rij = rad[i].c[j];
            demand(rij >= 0, "invalid search radius");
            if (rij > rmax) { rmax = rij; }
          }
      }
      
    if (rmax == 0)
      { /* Nothing to do: */ (*f2p) = f2(ni,0,p); }
    else
      { 
        /* Compute the initial scale {m} and its reduction factor {fscale}: */
        int m = 0;
        double fscale = 1.0;
        while (rmax > 0.5) { m = m+1; rmax = rmax/2; fscale = fscale/2; }
        /* Reduce the problem to scale {m}: */
        r2_t srad[ni];  /* Search radius at current scale. */
        r2_t step[ni];  /* Search step at current scale (if {quadopt} is false). */
        for (i = 0; i < ni; i++)
          { for (j = 0; j < 2; j++) 
              { srad[i].c[j] = rad[i].c[j]*fscale;
                p[i].c[j] = p[i].c[j]*fscale;
                step[i].c[j] = 0.5;
              }
          }
        /* Now solve the problem at increasing scales: */
        int s = m;
        while(TRUE)
          { /* Solve the problem at scale {s}: */
            if (quadopt)
              { float_image_align_single_scale_quadopt(ni, s, f2, srad, p, f2p); }
            else
              { float_image_align_single_scale_enum(ni, s, f2, srad, step, p, f2p); }
            /* Are we done? */
            if (s == 0) { break; }
            /* Expand to the next finer scale: */
            for (i = 0; i < ni; i++)
              { for (j = 0; j < 2; j++) 
                  { srad[i].c[j] = fmin(0.5, 2 * srad[i].c[j]);
                    p[i].c[j] = 2 * p[i].c[j];
                  }
              }
          
          }
      }
      
  }

double float_image_align_rel_disp_sqr(int ni, r2_t p[], int scale, r2_t q[], r2_t rad[])
  {
    double d2 = 0.0;
    double fscale = pow(2.0, scale);
    int i, j;
    for (i = 0; i < ni; i++)
      { for (j = 0; j < 2; j++)
          { double rij = rad[i].c[j];
            if (rij != 0)
              { double dij = (fscale*p[i].c[j] - q[i].c[j])/rij;
                d2 += dij*dij;
              }
          }
      }
    return d2;
  }
    
            
    

