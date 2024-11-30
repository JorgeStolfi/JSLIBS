/* See {hrn.h}. */
/* Last edited on 2024-11-23 07:58:40 by stolfi */

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

void rn_to_hrn(uint32_t n, double P[], double w, double p[])
  { p[0] = w;
    for (uint32_t i = 0;  i < n; i++) { p[i+1] = w*P[i]; }
  }

void hrn_to_rn(uint32_t n, double p[], double P[])
  { double w = p[0];
    demand(w != 0, "point is at infinity");
    for (uint32_t i = 0;  i < n; i++) { P[i] = p[i+1]/w; }
  }

sign_t hrn_side(uint32_t n, double p[], double h[])
  { double dd = rn_dot(n, p, h);
    return sign_double(dd);
  }

hrn_pmap_t hrn_pmap_alloc(uint32_t m, uint32_t n)
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

void hrn_pmap_compose(hrn_pmap_t *M, hrn_pmap_t *N, hrn_pmap_t *R)
  { demand(M->n == N->m, "incompatible domain/co-domain dimensions");
    uint32_t m = M->m, p = M->n, n = N->n;
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

void hrn_canonical_simplex(uint32_t d, uint32_t n, double p[])
  { uint32_t n1 = n+1;
    for (uint32_t i = 0;  i <= d; i++) 
      { uint32_t n1i = n1*i;
        for (uint32_t j = 0;  j <= n; j++)
          { p[n1i + j] = (i == j ? 1 : 0); }
      }
  }

void hrn_regular_simplex(uint32_t n, double p[])
  { double N = (double)n;
    double SN1 = sqrt(N+1);
    double c = (SN1 - 1)/N;
    double d = 1 + (N-1)*c;
    uint32_t n1 = n+1;
    /* Set the matrix {p}: */
    for (uint32_t i = 0;  i <= n; i++) 
      { uint32_t n1i = i*n1;
        p[n1i] = 1; /* Weight. */
        if (i == 0)
          { /* Set the first row to {(-1,-1,..-1)}: */
            for (uint32_t j = 1;  j <= n; j++) { p[n1i + j] = -1; }
          }
        else
          { /* Set row {i} to {(1+d+c)*u_{i-1} - (c,c,..c)}: */
            for (uint32_t j = 1;  j <= n; j++) { p[n1i + j] = (i == j ? d : -c); }
          }
      }
    }
