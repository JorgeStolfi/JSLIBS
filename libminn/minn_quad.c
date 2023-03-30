/* See {minn_quad.h}. */
/* Last edited on 2023-03-27 15:09:50 by stolfi */

#define _GNU_SOURCE
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

void minn_quad
  ( int32_t n,        /* Dimension of search space. */
    minn_goal_t *F,   /* Function to be minimized. */
    double dMax,      /* Radius of search domain, or {+INF}. */
    bool_t box,       /* True to search in the unit cube, false in the unit ball. */
    double tol,       /* Desired precision. */
    double v[],       /* (OUT) Minimum vector found. */
    double *Fval_P    /* (OUT) Goal function value at the minimum. */
  )
  {
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, ">> enter %s >>\n", __FUNCTION__); }

    demand(tol > 0, "invalid {tol}");

    /* The initial guess is {(0,0,...0)}: */
    rn_zero(n, v);
    double Fv = F(n, v);
    
    if (n > 0)
      { 
        /* Optimize: */
        sign_t dir = -1; /* Look for minimum. */
        int32_t maxIters = 10;
        double rIni = 0.5*dMax;   /* Initial probe simplex radius. */
        double rMin = tol;        /* Minimum probe simplex radius. */
        double rMax = 0.70*dMax;  /* Maximum probe simplex radius. */
        double stop = 0.25*tol;   /* Stop when {x} moves less than this. */
        sve_minn_iterate
          ( n, 
            F, NULL, 
            v, &Fv,
            dir, dMax, box, rIni, rMin, rMax, stop,
            maxIters,
            debug
          );

        (*Fval_P) = Fv;

      }
    else
      { /* Just compute the goal for the given vector: */
        (*Fval_P) = F(n, v);
      }
    
    if (debug) { fprintf(stderr, "<< leave %s <<\n", __FUNCTION__); }
  }
