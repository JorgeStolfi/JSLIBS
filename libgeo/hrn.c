/* See {hrn.h}. */
/* Last edited on 2023-02-25 15:29:35 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <hrn.h>

#include <rn.h>
#include <rmxn.h>
#include <affirm.h>
#include <sign.h>
#include <sign_get.h>
#include <bool.h>

#include <math.h>
#include <stdlib.h>
#include <assert.h>

/* Based on HR2.m3 and HR3.m3, created 1994-05-04 by J. Stolfi. */

/* !!! We should use Kahan's summation for all scalar products. */

void rn_to_hrn(int32_t n, double P[], double w, double p[])
  { int32_t i;
    p[0] = w;
    for (i = 0; i < n; i++) { p[i+1] = w*P[i]; }
  }

void hrn_to_rn(int32_t n, double p[], double P[])
  { int32_t i;
    double w = p[0];
    demand(w != 0, "point is at infinity");
    for (i = 0; i < n; i++) { P[i] = p[i+1]/w; }
  }

sign_t hrn_side(int32_t n, double p[], double h[])
  { double dd = rn_dot(n, p, h);
    return sign_double(dd);
  }

hrn_pmap_t hrn_pmap_alloc(int32_t m, int32_t n)
  { hrn_pmap_t M;
    M.m = m; M.n = n;
    M.dir = rmxn_alloc(m+1,n+1);
    M.inv = rmxn_alloc(m+1,n+1);
    return M;
  }

void hrn_pmap_free(hrn_pmap_t M)
  { free(M.dir);
    free(M.inv);
  }

void hrn_map_point(double p[], hrn_pmap_t *M, double q[])
  { rmxn_map_row(M->m, M->n, p, M->dir, q); }

void hrn_pmap_comp(hrn_pmap_t *M, hrn_pmap_t *N, hrn_pmap_t *R)
  { demand(M->n == N->m, "incompatible domain/co-domain dimensions");
    int32_t m = M->m, p = M->n, n = N->n;
    demand((R->m == m) && (R->n == n), "result map has wrong dimensions");
    rmxn_mul(m, p, n, M->dir, N->dir, R->dir);
    rmxn_mul(n, p, m, N->inv, M->inv, R->inv);
  }

hrn_pmap_t hrn_pmap_inv(hrn_pmap_t *M)
  { hrn_pmap_t R;
    R.m = M->n;
    R.n = M->m;
    R.dir = M->inv;
    R.inv = M->dir;
    return R;
  }

void hrn_canonical_simplex(int32_t d, int32_t n, double p[])
  { int32_t n1 = n+1;
    int32_t i, j;
    for (i = 0; i <= d; i++) 
      { int32_t n1i = n1*i;
        for (j = 0; j <= n; j++)
          { p[n1i + j] = (i == j ? 1 : 0); }
      }
  }

void hrn_regular_simplex(int32_t n, double p[])
  { double N = (double)n;
    double SN1 = sqrt(N+1);
    double c = (SN1 - 1)/N;
    double d = 1 + (N-1)*c;
    int32_t i, j;
    int32_t n1 = n+1;
    /* Set the matrix {p}: */
    for (i = 0; i <= n; i++) 
      { int32_t n1i = i*n1;
        p[n1i] = 1; /* Weight. */
        if (i == 0)
          { /* Set the first row to {(-1,-1,..-1)}: */
            for (j = 1; j <= n; j++) { p[n1i + j] = -1; }
          }
        else
          { /* Set row {i} to {(1+d+c)*u_{i-1} - (c,c,..c)}: */
            for (j = 1; j <= n; j++) { p[n1i + j] = (i == j ? d : -c); }
          }
      }
    }
