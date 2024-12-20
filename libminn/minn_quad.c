/* See {minn_quad.h}. */
/* Last edited on 2024-12-05 13:14:39 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
 
#include <bool.h>
#include <rn.h>
#include <sign.h>
#include <jsmath.h>
#include <affirm.h>
#include <minn.h>

#include <sve_minn.h>

#include <minn_quad.h>

/* !!! Add {OK} and {Proj} parameters. !!! */

void minn_quad
  ( uint32_t n,        /* Dimension of search space. */
    minn_goal_t *F,   /* Function to be minimized. */
    bool_t box,       /* True to search in the unit cube, false in the unit ball. */
    double tol,       /* Desired precision. */
    double v[],       /* (OUT) Minimum vector found. */
    double *Fval_P    /* (OUT) Goal function value at the minimum. */
  )
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, ">> enter %s >>\n", __FUNCTION__); }

    demand(tol > 0, "invalid {tol}");

    /* The initial guess is {(0,0,...0)}: */
    rn_zero(n, v);
    double Fv = F(n, v);
    
    if (n > 0)
      { 
        /* Optimize: */
        sign_t dir = -1; /* Look for minimum. */
        uint32_t maxIters = 10;
        double *ctr = NULL;        /* Search domain center is the origin. */
        double dMax = 1.0;         /* Search domain radius. */
        double rIni = 0.5;         /* Initial probe simplex radius. */
        double rMin = tol;         /* Minimum probe simplex radius. */
        double rMax = 0.70*dMax;   /* Maximum probe simplex radius. */
        double minStep = 0.25*tol; /* Stop when {x} moves less than this. */
        bool_t sve_debug = debug;
        bool_t sve_debug_probes = FALSE; 
        
        sve_minn_iterate
          ( n, F, NULL, NULL,
            v, &Fv, dir,
            ctr, dMax, box, rIni, rMin, rMax, minStep,
            maxIters,
            sve_debug, sve_debug_probes
          );

        (*Fval_P) = Fv;

      }
    else
      { /* Just compute the goal for the given vector: */
        (*Fval_P) = F(n, v);
      }
    
    if (debug) { fprintf(stderr, "<< leave %s <<\n", __FUNCTION__); }
  }
