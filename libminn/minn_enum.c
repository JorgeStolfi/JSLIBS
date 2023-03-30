/* See {minn_enum.h}. */
/* Last edited on 2023-03-27 15:09:46 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
 
#include <bool.h>
#include <jsmath.h>
#include <affirm.h>
#include <rn.h>
#include <rmxn.h>

#include <minn.h>

#include <minn_enum.h>

void minn_enum
  ( int32_t n,        /* Dimension of search space. */
    minn_goal_t *F,   /* Function to be minimized. */
    double dMax,      /* Radius of domain, or {+INF}. */ 
    bool_t box,       /* True to search in the unit cube, false in the unit ball. */
    double tol[],     /* Desired precision. */
    double v[],       /* (OUT) Minimum vector found. */
    double *Fval_P    /* (OUT) Goal function value at the minimum. */
  )
  {
    demand(tol > 0, "invalid {tol}");
    demand(isfinite(dMax) && (dMax >= 0), "invalid {dMax}");
    bool_t debug = FALSE;
    
    /* Compute the number of samples {ns[i]} along each half-axis {i}: */
    double fudge = 1.0e-10; /* Fudge domain expansion to ensure grid edge is in. */
    int32_t ns[n];
    for (int32_t i = 0; i < n; i++) { ns[i] = (int32_t)floor(dMax/tol[i] + fudge); }
    
    /* Enumerate all integer tuples {t[0..n-1]} where {t[k]} ranges in {-ns..+ns}: */
    (*Fval_P) = +INF;  /* Minimum goal found so far. */
    int32_t t[n];      /* Enumeration variables. */
    for (int32_t i = 0; i < n; i++) { t[i] = -ns[i]; }
    double u[n];       /* Sampling vector. */
    int32_t knext = -1; /* Next tuple elem to be incremented is {t[knext]}. */
    while (knext < n)
      { if (debug)
          { fprintf(stderr, "  t = ( ");
            for (int32_t k = 0; k < n; k++) { fprintf(stderr, " %d", t[k]); }
            fprintf(stderr, " )\n");
          }
        /* !!! Should avoid generating tuples outside {\RF} !!! */

        /* Build the sample vector {u}: */
        for (int32_t i = 0; i < n; i++) 
          { u[i] = ((double)t[i])*tol[i];
            assert(fabs(u[i]) <= dMax + 2*fudge*tol[i]);
          }

        /* Determine whether the point is inside the domain: */
        bool_t inside = (box ? TRUE : rn_norm_sqr(n, u) <= 1 + fudge);
        if (inside)
          { /* Evaluate the goal function at {v}: */
            double Fval = F(n, u);
            if (Fval < (*Fval_P))
              { /* Update the current optimum: */
                rn_copy(n, u, v);
                (*Fval_P) = Fval;
              }
          }

        /* Get the next tuple {t[0..nd-1]}: */
        if ((knext < 0) || (t[knext] >= ns[knext]))
          { /* Can't increment this, try next one: */
            knext++;
          }
        else
          { /* Increment {t[knext]} and clear all previous {t[k]}: */
            t[knext]++;
            while (knext > 0) { knext--; t[knext] = -ns[knext]; }
          }
      }
  }
