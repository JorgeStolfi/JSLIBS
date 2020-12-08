/* See rmxn_extra.h. */
/* Last edited on 2019-12-18 21:50:21 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <rn.h>
#include <rmxn.h>
#include <jsrandom.h>
#include <affirm.h>
#include <cmp.h>

#include <rmxn_extra.h>

void rmxn_perturb_unif(int m, int n, double mag, double M[])
  {
    int i, j;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++) 
          { double d = dabrandom(-mag, +mag);
            M[i*n + j] +=  d; 
          }
      }
  }

void rmxn_rot2(double c, double s, double *x, double *y);
  /* Rotates the vector {(x,y)} by an angle {theta}, given
     {x}, {y}, {c = cos(theta)}, and {s = sin(theta)}. */

void rmxn_throw_ortho(int n, double M[])
  { /* Start with the identity matrix: */
    rmxn_ident(n, n, M);
    /* Now apply mirroring ops to it: */
    int k, i, j;
    double s[n];
    int flip = 0;
    for (k = 0; k < n; k++)
      { /* Now the first {k} rows and cols of {M} are a random orthonormal {k×k} matrix. */
        double *Mk = &(M[n*k]); /* Row {k} of {M} */
        /* Pick a random direction {v} in {R^{k+1}}, store it in row {k} of {M}: */
        rn_throw_dir(k+1, Mk);
        /* We now apply to rows {0..k-1} a reflection that takes {u_k} to {v}. */
        /* Compute the unit vector {s[0..k]} normal to the bisector of {u_k} and {v}: */
        double s2 = 0;
        for (j = 0; j <= k; j++)
          { double sj = Mk[j] - (j == k ? 1 : 0); 
            s[j] = sj; 
            s2 += sj*sj;
          }
        if (s2 != 0)
          { /* Reflection is not trivial. */
            double sm = sqrt(s2);
            for (j = 0; j <= k; j++) { s[j] /= sm; }
            /* Apply reflection along {s} to each row: */
            for (i = 0; i < k; i++)
              { double *Mi = &(M[i*n]);
                (void)rn_mirror(k+1, Mi, s, Mi);
                flip = 1 - flip;
              }
          }
      }
    /* Now flip the first row if needed to keep the determinant positive: */
    if (flip != 0) { for (j = 0; j < n; j++) { M[j] = -M[j]; } }
  }

void rmxn_spin_rows(int m, int n, double A[], double M[])
  { /* Generate a random orthonormal {n×n} matrix {N}: */
    double N[n*n];
    rmxn_throw_ortho(n, N);
    /* Map each row of {A} by {N} (beware of aliasing between {A} and {M}): */
    double v[n];
    int i;
    for (i = 0; i < m; i++)
      { rmxn_map_row(n, n, &(A[i*n]), N, v);
        rn_copy(n, v, &(M[i*n]));
      }
  }

void rmxn_spin_cols(int m, int n, double A[], double M[])
  { /* Generate a random orthonormal {m×m} matrix {N}: */
    double N[m*m];
    rmxn_throw_ortho(m, N);
    /* Map each col of {A} by {N} (beware of aliasing between {A} and {M}): */
    double a[m], v[m];
    int j;
    for (j = 0; j < n; j++)
      { rmxn_get_col(m, n, A, j, a);
        rmxn_map_col(m, m, a, N, v);
        rmxn_set_col(m, n, M, j, v);
      }
  }
  
void rmxn_shift_rows(int m, int n, double A[], double v[], double M[])
  { int i, j;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++) 
          { int ij = i*n + j; M[ij] = A[ij] + v[j]; }
      }
  }
  
void rmxn_shift_cols(int m, int n, double v[], double A[], double M[])
  { int i, j;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++) 
          { int ij = i*n + j; M[ij] = A[ij] + v[i]; }
      }
  }

void rmxn_canonical_simplex(int d, int n, double V[])
  { 
    demand(0 <= d, "bad dimension {d}");
    demand(d < n, "space dimension {n} is too small for {d}");
    int i, j;
    for (i = 0; i <= d; i++) 
      { for (j = 0; j < n; j++)
          { V[i*n + j] = (i == j ? 1 : 0); }
      }
  }

double rmxn_canonical_simplex_radius(int d)
  { double D = (double)d;
    return sqrt(D/(D+1));
  }

double rmxn_canonical_simplex_subradius(int d, int k)
  { double D = (double)d;
    double K = (double)k;
    return sqrt((D-K)/((D+1)*(K+1)));
  }
  
double rmxn_canonical_simplex_edge(int d)  
  { return M_SQRT2; }
  
double rmxn_canonical_simplex_height(int d)
  { double D = (double)d;
    return sqrt((D+1)/D);
  }

double rmxn_canonical_simplex_measure(int d)
  { double D = (double)d;
    return sqrt(D+1)*exp(-lgamma(D+1));
  }

void rmxn_regular_simplex(int n, double V[])
  { double N = (double)n;
    double SN1 = sqrt(N+1);
    double c = (SN1 - 1)/N;
    double d = 1 + (N-1)*c;
    int i, j;
    /* Set the matrix {p}: */
    for (i = 0; i <= n; i++) 
      { int ni = i*n;
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

double rmxn_regular_simplex_radius(int n)
  { double N = (double)n;
    return sqrt(N);
  }

double rmxn_regular_simplex_subradius(int n, int k)
  { double N = (double)n;
    double K = (double)k;
    return sqrt((N-K)/(K+1));
  }
  
double rmxn_regular_simplex_edge(int n)  
  { double N = (double)n;
    return sqrt(2*(N+1));
  }
  
double rmxn_regular_simplex_height(int n)
  { double N = (double)n;
    return (N+1)/sqrt(N);
  }
  
double rmxn_regular_simplex_measure(int n)
  { double N = (double)n;
    return exp((N+1)*log(N+1)/2 - lgamma(N+1));
  }

void rmxn_throw_canonical_simplex(int d, double x[])
  { /* Generate a random point in the unit cube {[0_1]^d}: */
    int i;
    for (i = 0; i < d; i++) { x[i] = drandom(); }
    /* Sort {x} by increasing value: */
    auto int cmp(const void *a, const void *b);
    int cmp(const void *a, const void *b) { return cmp_double((double *)a, (double *)b); }
    qsort(x, d, sizeof(double), &cmp);
    /* Now map the ordered {d}-simplex to the canonical {d}-simplex: */
    x[d] = 1;
    for (i = d; i > 0; i--) { x[i] = x[i] - x[i-1]; }
  }
  
void rmxn_throw_canonical_simplex_ball(int d, double x[]) 
  { /* Generate a random point in the unit {(d+1)}-dimensional ball: */
    rn_throw_ball(d+1, x);
    /* Compute projection {s} of {x} on {u = (1,1,..1)/sqrt(d+1)}: */
    double sum = 0.0;
    int i;
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
