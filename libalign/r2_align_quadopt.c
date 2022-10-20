/* See {r2_align_quadopt.h}. */
/* Last edited on 2022-10-20 06:16:29 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
 
#include <bool.h>
#include <sign.h>
#include <r2.h>
#include <i2.h>
#include <jsmath.h>
#include <affirm.h>
#include <sve_minn.h>
#include <r2_align.h>

#include <r2_align_quadopt.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

void r2_align_quadopt
  ( int32_t ni,                   /* Number of objects to align. */
    r2_align_mismatch_t *F2,  /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],              /* Max alignment adjustment for each object. */
    double tol,               /* Desired precision. */
    r2_t p[],                 /* (IN/OUT) Corresponding points in each object. */
    double *F2valP             /* (OUT) Mismatch for the computed alignment vector. */
  )
  {
    int32_t maxIters = 10;
    bool_t debug = TRUE;
    
    demand(tol > 0, "invalid {tol}");
    int32_t nd = r2_align_count_degrees_of_freedom(ni, arad); 

    /* Save the initial guess: */
    r2_t ctr[ni]; /* Center of ellipsoid (aved initial guess). */
    for (int32_t i = 0; i < ni; i++) { ctr[i] = p[i]; }

    if (nd > 0)
      { 
        /* Compute the axes and radii of the search ellipsoid {\RF}: */
        r2_t U[nd*ni]; /* Alignment adjustment vectrors. */
        double urad[nd];
        r2_align_compute_search_ellipsoid (ni, arad, nd, U, urad);

        /* The optimization variables {x[0..nd-1]} are the coordinates in the {U} basis
          of the adjustments from {ctr} to the sampling point {psmp},
          divided by the corresponding radii {urad[0..nd-1]}. */

        /* Compute the tolerance {xtol} in the {x} variables space: */
        double xtol = +INF;
        for (int32_t k = 0; k < nd; k++) 
          { double xtolk = tol/urad[k]; if (xtolk < xtol) { xtolk = xtol; } }

        auto void compute_sample_point(int32_t nz, double z[]);
          /* Converts the minimization variables {z[0..nz-1]} to the corresponding
            sampling alignment and stores it in {p[0..ni-1]}. Expects {nz == nd}. */
        
        auto double sve_goal(int32_t nz, double z[]);
          /* Computes the minimization goal function from the given argument {z[0..nv-1]}.
            Expects {nz == nd}. Also sets {p[0..ni-1]}. Assumes  that the initial guess was
            saved in {ctr[0..ni-1]}. */

        /* Compute the initial goal function value: */
        double x[nd]; for (int32_t k = 0; k < nd; k++) { x[k] = 0.0; }
        double Fx = sve_goal(nd, x);
        sign_t dir = -1; /* Look for minimum. */
        double rMax = 0.5;                /* Max simplex radius. */
        double rMin = fmin(0.05, 3*xtol); /* Min simplex radius. */
        sve_minn_iterate
          ( nd, 
            &sve_goal, NULL, 
            x, &Fx,
            dir,
            /*dMax:*/ 1.0,
            /*dBox:*/ FALSE,
            /*rIni:*/ 0.5, 
            /*rMin:*/ rMin, 
            /*rMax:*/ rMax, 
            /*stop:*/ xtol,
            maxIters,
            debug
          );

        /* Return the optimal vector and its mismatch value: */
        compute_sample_point(nd, x);
        (*F2valP) = Fx;

        /* Local implementations: */

        void compute_sample_point(int32_t nz, double z[])
          { assert(nz == nd);
            for (int32_t i = 0; i < ni; i++)
              { for (int32_t j = 0; j < 2; j++)
                  { p[i].c[j] = ctr[i].c[j];
                    double rij = arad[i].c[j];
                    if (rij != 0)
                      { double dij = 0;
                        for (int32_t k = 0; k < nd; k++) 
                          { r2_t *uk = &(U[k*ni]);
                            dij += z[k]*urad[k]*uk[i].c[j];
                          }
                        p[i].c[j] += dij;
                      }
                  }
              }
          }
        
        double sve_goal(int32_t nz, double z[])
          { assert(nz == nd);
            /* Convert variables {z[0..nx-1]} to adjustments {p[0..ni-1]}: */
            compute_sample_point(nz, z);
            /* Evaluate the client function: */
            double F2val = F2(ni, p);
            return F2val;
          }
      }
    else
      { (*F2valP) = F2(ni, ctr); }
    
    return;
            
  }
            
    

