/* See {sve_minn_iterate.h} */
/* Last edited on 2025-04-01 09:09:37 by stolfi */

#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <affirm.h>
#include <vec.h>
#include <jsmath.h>
#include <rmxn.h>
#include <rmxn_spin.h>
#include <rmxn_shift.h>
#include <rmxn_regular_simplex.h>
#include <gausol_print.h>
#include <gausol_solve.h>

#include <sve_minn.h>

#include <sve_minn_iterate.h>

#define Pr fprintf
#define Er stderr

/* INTERNAL PROTOTYPES */

void sve_minn_iterate
  ( uint32_t n, 
    sve_goal_t *F, 
    sve_pred_t *OK,
    sve_proj_t *Proj,
    double x[],
    double *FxP,
    sign_t dir, 
    double dCtr[],
    double dMax,
    bool_t dBox,
    double rIni,
    double rMin, 
    double rMax,
    double minStep,
    uint32_t maxIters,
    bool_t debug,
    bool_t debug_probes
  )
  { 
    if (debug_probes) { /* Implies {debug}: */ debug = TRUE; }
    if (debug) { Pr(Er, ">> enter %s >>\n", __FUNCTION__); }
    
    demand(rMin >= sve_minn_MIN_RADIUS, "invalid {rMin}, too small");
    demand(rMin <= rMax, "invalid {rMin > rMax}");
    demand(rMin <= dMax, "invalid {rMin > dMax}");
    demand((rMin <= rIni) && (rIni <= rMax), "invalid {rIni}");
    
    /* Allocate storage for the simplex: */
    uint32_t nv = n+1; /* Number of vertices in simplex. */
    double *v = talloc(nv*n, double); /* Cartesian coords of simplex vertices. */
    
    /* Allocate storage for the sample values: */
    uint32_t nf = (n+1)*(n+2)/2; /* Number of probe points. */
    double *Fv = talloc(nf, double); /* Sampled function values. */

    /* Iteration variables: */
    double radius = rIni; /* Current probe simplex radius: */
    double dPrev = dMax; /* Distance moved in previous iteration ({dMax} if none). */
    uint32_t nIters = 0; /* Counts quadratic step iterations. */
    uint32_t nEvals = 0; /* Counts function evaluations. */
    
    /* Get initial function value: */
    double Fx = (*FxP);
    
    while(TRUE) 
      { if (debug) { Pr(Er, "  iteration %4d  function = %22.16e\n", nIters, Fx); }
        double Fx_save = Fx; /* Before {Proj}. */
        if (Proj != NULL)
          { /* Apply client's projection: */
            Fx = Proj(n, x, Fx);
            if (debug) 
              { Pr(Er, " applyed {Proj}");
                if (Fx != Fx_save) { Pr(Er, " function changed %22.16e --> %22.16e", Fx_save, Fx); }
                Pr(Er, "\n");    
              }
          }
        /* Check budget: */
        if (nIters >= maxIters) 
          { if (debug) { Pr(Er, "  iteration limit exhausted\n"); }
            break;
          }
        
        /* Compute distance to domain center: */
        double dist;
        if (dCtr == NULL)
          { dist = (dBox ? rn_L_inf_norm(n, x) : rn_norm(n, x)); }
        else
          { dist = (dBox ? rn_L_inf_dist(n, dCtr, x) : rn_dist(n, dCtr, x)); }

        /* Check for termination: */
        if ((OK != NULL) && (OK(nIters, n, x, Fx, dist, dPrev, radius))) 
          { if (debug) { Pr(Er, "  client is satisfied\n"); }
            break;
          }
        
        double dStep; /* Distance moved. */
        uint32_t stepKind = sve_minn_single_step
          ( n, x, &Fx, F, dir, dCtr, dMax, dBox, &radius, rMin, 
            v, Fv, &nEvals, &dStep, debug, debug_probes
          );
        nIters++;
        if (debug) 
          { char *choice = ((char*[3]){"quadratic min", "sample point", "old center"})[stepKind];
            Pr(Er, "  chosing %s  step len = %22.16e function = %22.16e\n", choice, dStep, Fx);
            if (debug_probes)
              { Pr(Er, "    new guess = \n");
                rn_gen_print(Er, n, x, "%20.16f", "    [ ", "\n      ", " ]\n");
              }
          }
          
        /* Estimate the distance {dEst} to the minimum and choose new radius {rNew}: */
        double dEst, rNew;  
        if (stepKind == 0)
          { /* Quadratic minimization succeeded, assume geometric convergence: */
            dEst = (dPrev <= dStep ? +INF : dStep*dStep/(dPrev - dStep));
            rNew = dStep/2;
          }
        else if (stepKind == 1)
          { /* Moved to a simplex probe point: */
            dEst = +INF;
            rNew = radius/2;
          }
        else if (stepKind == 2)
          { /* Current guess is still the optimum: */
            dEst = radius;
            rNew = radius/2;
          }
        else
          { assert(FALSE); }
        if (rNew < 0.1*radius) { rNew = 0.1*radius; }
        if (rNew > 2.0*radius) { rNew = 2.0*radius; }
        if (debug)  { Pr(Er, "  next simplex radius = %12.8f", rNew); }
        if (rNew < rMin)
          { rNew = rMin;
            if (debug)  { Pr(Er, " augmented to {rMin} = %12.8f", rNew); }
          }
        if (rNew > rMax) 
          { rNew = rMax; 
            if (debug)  { Pr(Er, " reduced to {rMax} = %12.8f", rNew); }
          }
        if (debug) { Pr(Er, "\n"); }

        /* Check for convergence: */
        bool_t dEst_is_small = (dEst < minStep);
        bool_t no_progress = ((dStep < minStep) && (radius <= rMin));
        if (dEst_is_small || no_progress) 
          { if (debug) 
              { Pr(Er, "  seems to have converged:\n");
                Pr(Er, "    dPrev = %22.16e\n", dPrev);
                Pr(Er, "    dStep = %22.16e\n", dStep);
                Pr(Er, "    dEst =  %22.16e\n", dEst);
              }
            break;
          }
        /* The motion is at least the center-to-edge radius of the simplex: */
        /* Remember how much we moved in this iteration: */
        dPrev = dStep;
        radius = rNew;
      }
    if (debug) { Pr(Er, "  did %d iterations and %d function evaluations.\n", nIters, nEvals); }
    free(v);
    free(Fv);
    
    /* Return final function value: */
    (*FxP) = Fx;
    if (debug) { Pr(Er, "<< leave %s <<\n", __FUNCTION__); }
  }
