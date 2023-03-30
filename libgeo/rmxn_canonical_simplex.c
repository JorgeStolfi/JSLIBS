/* See rmxn_canonical_simplex.h. */
/* Last edited on 2023-03-26 16:41:37 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <rn.h>
#include <rmxn.h>
#include <jsrandom.h>
#include <affirm.h>
#include <cmp.h>

#include <rmxn_extra.h>
#include <rmxn_canonical_simplex.h>

void rmxn_canonical_simplex(int32_t d, int32_t n, double V[])
  { 
    demand(0 <= d, "bad dimension {d}");
    demand(d < n, "space dimension {n} is too small for {d}");
    int32_t i, j;
    for (i = 0; i <= d; i++) 
      { for (j = 0; j < n; j++)
          { V[i*n + j] = (i == j ? 1 : 0); }
      }
  }

double rmxn_canonical_simplex_radius(int32_t d)
  { double D = (double)d;
    return sqrt(D/(D+1));
  }

double rmxn_canonical_simplex_subradius(int32_t d, int32_t k)
  { double D = (double)d;
    double K = (double)k;
    return sqrt((D-K)/((D+1)*(K+1)));
  }
  
double rmxn_canonical_simplex_edge(int32_t d)  
  { return M_SQRT2; }
  
double rmxn_canonical_simplex_height(int32_t d)
  { double D = (double)d;
    return sqrt((D+1)/D);
  }

double rmxn_canonical_simplex_measure(int32_t d)
  { double D = (double)d;
    return sqrt(D+1)*exp(-lgamma(D+1));
  }


void rmxn_canonical_simplex_throw(int32_t d, double x[])
  { /* Generate a random point in the unit cube {[0_1]^d}: */
    int32_t i;
    for (i = 0; i < d; i++) { x[i] = drandom(); }
    /* Sort {x} by increasing value: */
    auto int32_t cmp(const void *a, const void *b);
    int32_t cmp(const void *a, const void *b) { return cmp_double((double *)a, (double *)b); }
    qsort(x, d, sizeof(double), &cmp);
    /* Now map the ordered {d}-simplex to the canonical {d}-simplex: */
    x[d] = 1;
    for (i = d; i > 0; i--) { x[i] = x[i] - x[i-1]; }
  }
  
void rmxn_canonical_simplex_ball_throw(int32_t d, double x[]) 
  { /* Generate a random point in the unit {(d+1)}-dimensional ball: */
    rn_throw_ball(d+1, x);
    /* Compute projection {s} of {x} on {u = (1,1,..1)/sqrt(d+1)}: */
    double sum = 0.0;
    int32_t i;
    for (i = 0; i <= d; i++) { sum += x[i]; }
    double s = sum/sqrt(d+1);
    /* Ensure projection is in {[-1_+1]}: */
    if (s > +1.0) { s = +1.0; }
    if (s < -1.0) { s = -1.0; }
    /* Compute radius of ball slice at coordinate {s}: */
    double rx = sqrt(1.0 - s*s);
    /* Expand unit ball to cylinder with axis {u} and project normal to {u}: */
    if (rx == 0) { rx = 1; }
    double q = sum/(d+1);
    for (i = 0; i <= d; i++) { x[i] = (x[i] - q)/rx; }
    /* Scale to the radius of the canonical {d}-simplex: */
    double rc = rmxn_canonical_simplex_radius(d);
    rn_scale(d+1, rc, x, x);
    /* Shift to center in {(1,1,..1)/(d+1)}: */
    rn_shift(d+1, 1.0/(d+1), x, x);
  }
