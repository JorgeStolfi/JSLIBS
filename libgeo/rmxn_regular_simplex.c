/* See rmxn_regular_simplex.h. */
/* Last edited on 2023-03-26 16:43:24 by stolfi */

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

#include <rmxn_regular_simplex.h>

void rmxn_regular_simplex(int32_t n, double V[])
  { double N = (double)n;
    double SN1 = sqrt(N+1);
    double c = (SN1 - 1)/N;
    double d = 1 + (N-1)*c;
    int32_t i, j;
    /* Set the matrix {p}: */
    for (i = 0; i <= n; i++) 
      { int32_t ni = i*n;
        if (i == 0)
          { /* Set the first row to {(-1,-1,..-1)}: */
            for (j = 0; j < n; j++) { V[ni + j] = -1; }
          }
        else
          { /* Set row {i} to {(1+d+c)*u_{i-1} - (c,c,..c)}: */
            for (j = 0; j < n; j++) { V[ni + j] = (i == j+1 ? d : -c); }
          }
      }
    }

double rmxn_regular_simplex_radius(int32_t n)
  { double N = (double)n;
    return sqrt(N);
  }

double rmxn_regular_simplex_subradius(int32_t n, int32_t k)
  { double N = (double)n;
    double K = (double)k;
    return sqrt((N-K)/(K+1));
  }
  
double rmxn_regular_simplex_edge(int32_t n)  
  { double N = (double)n;
    return sqrt(2*(N+1));
  }
  
double rmxn_regular_simplex_height(int32_t n)
  { double N = (double)n;
    return (N+1)/sqrt(N);
  }
  
double rmxn_regular_simplex_measure(int32_t n)
  { double N = (double)n;
    return exp((N+1)*log(N+1)/2 - lgamma(N+1));
  }
