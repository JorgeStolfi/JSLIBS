/* See {r2_align_quadopt.h}. */
/* Last edited on 2024-12-05 10:20:11 by stolfi */

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
  ( int32_t ni,               /* Number of objects to align. */
    r2_align_mismatch_t *F2,  /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],              /* Max delta vector coordinates for each object. */
    bool_t bal,               /* True if alignment vector adjustments should be balanced. */
    double tol,               /* Desired precision. */
    r2_t p[],                 /* (IN/OUT) Corresponding points in each object. */
    double *F2val_P           /* (OUT) Mismatch for the computed alignment vector. */
  )
  {
    int32_t maxIters = 10;
    bool_t debug = FALSE;
    
    demand(tol > 0, "invalid {tol}");
    i2_t nv = r2_align_count_variable_coords (ni, arad);
    int32_t nd = r2_align_count_degrees_of_freedom(nv, bal);

    /* Save the initial guess: */
    r2_t pctr[ni]; /* Center of {\RD} ellipsoid (saved initial guess). */
    for (uint32_t i = 0;  i < ni; i++) { pctr[i] = p[i]; }

    if (nd > 0)
      { 
        /* Compute the axes and radii of the search ellipsoid {\RF}. It
          is a section of the ellipsoid {\RD} by the hyperplanes that
          represent the X and Y balancing and conformation constraints
          for the deltas {p[i]-pctr[i]}. */
        r2_t U[nd*ni]; /* Delta vector basis. */
        double urad[nd];
        r2_align_compute_search_ellipsoid (ni, arad, bal, nd, U, urad);

        /* The optimization variables {x[0..nd-1]} are the coordinates in the {U} basis
          of balanced conformal delta vectors, divided by the corresponding radii
          {urad[0..nd-1]}. */

        /* Compute the tolerance {xtol} in the {x} variables space: */
        double xtol = +INF;
        for (uint32_t k = 0;  k < nd; k++) 
          { double xtolk = tol/urad[k]; if (xtolk < xtol) { xtolk = xtol; } }

        auto void compute_sample_alignment(int32_t nxt, double xt[]);
          /* Converts the minimization variables {xt[0..nxt-1]} to the corresponding
            sampling alignment and stores it in {p[0..ni-1]}. Expects {nxt == nd}. */
        
        auto double sve_goal(int32_t nxt, double xt[]);
          /* Computes the minimization goal function from the given
            minimization variables {xt[0..nv-1]}. Expects {nxt == nd}.
            Also sets {p[0..ni-1]}. Assumes that the initial alignment
            guess was saved in {pctr[0..ni-1]}. */

        /* Compute the initial goal function value: */
        double x[nd]; for (uint32_t k = 0;  k < nd; k++) { x[k] = 0.0; }
        double Fx = sve_goal(nd, x);
        
        sign_t dir = -1; /* Look for minimum. */
        double rMax = 0.5;                /* Max simplex radius. */
        double rMin = fmin(0.05, 3*xtol); /* Min simplex radius. */
        bool_t sve_debug = debug;
        bool_t sve_debug_probes = FALSE; 
        
        sve_minn_iterate
          ( /*n:*/       nd, 
            /*F:*/       &sve_goal,
            /*OK:*/      NULL, 
            /*Proj:*/    NULL,
            /*x:*/       x,
            /*FxP:*/     &Fx,
            /*dir:*/     dir,
            /*ctr:*/     NULL,
            /*dMax:*/    1.0,
            /*box:*/     FALSE,
            /*rIni:*/    0.5, 
            /*rMin:*/    rMin, 
            /*rMax:*/    rMax, 
            /*minStep:*/ xtol,
            maxIters,
            sve_debug, sve_debug_probes
          );

        /* Return the optimal vector and its mismatch value: */
        compute_sample_alignment(nd, x);
        (*F2val_P) = Fx;

        /* Local implementations: */

        void compute_sample_alignment(int32_t nxt, double xt[])
          { assert(nxt == nd);
            for (uint32_t i = 0;  i < ni; i++)
              { for (uint32_t j = 0;  j < 2; j++)
                  { p[i].c[j] = pctr[i].c[j];
                    double rij = arad[i].c[j];
                    if (rij != 0)
                      { double dij = 0;
                        for (uint32_t k = 0;  k < nd; k++) 
                          { r2_t *uk = &(U[k*ni]);
                            dij += xt[k]*urad[k]*uk[i].c[j];
                          }
                        p[i].c[j] += dij;
                      }
                  }
              }
          }
        
        double sve_goal(int32_t nxt, double xt[])
          { assert(nxt == nd);
            /* Convert variables {xt[0..nxt-1]} to the delta vector {p[0..ni-1]}: */
            compute_sample_alignment(nxt, xt);
            /* Evaluate the client function: */
            double F2val = F2(ni, p);
            return F2val;
          }
      }
    else
      { (*F2val_P) = F2(ni, pctr); }
    
    return;
            
  }
            
    

