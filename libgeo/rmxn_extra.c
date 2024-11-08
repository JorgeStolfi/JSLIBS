/* See rmxn_extra.h. */
/* Last edited on 2024-11-07 22:37:56 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <rn.h>
#include <rmxn.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>
#include <cmp.h>

#include <rmxn_extra.h>

void rmxn_throw(int32_t m, int32_t n, double M[])
  { for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++) 
          { M[i*n + j] += dabrandom(-1.0, +1.0); }
      }
  }

void rmxn_perturb_unif(int32_t m, int32_t n, double pabs, double prel, double M[])
  {
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++) 
          { double *Mij = &(M[i*n + j]);
            double mag = pabs + prel*fabs(*Mij);
            double d = dabrandom(-mag, +mag);
            (*Mij) += d; 
          }
      }
  }

void rmxn_rot2(double c, double s, double *x, double *y);
  /* Rotates the vector {(x,y)} by an angle {theta}, given
     {x}, {y}, {c = cos(theta)}, and {s = sin(theta)}. */

void rmxn_throw_ortho(int32_t n, double M[])
  { /* Start with the identity matrix: */
    rmxn_ident(n, n, M);
    /* Now apply mirroring ops to it: */
    double s[n];
    int32_t flip = 0;
    for (int32_t k = 0; k < n; k++)
      { /* Now the first {k} rows and cols of {M} are a random orthonormal {k×k} matrix. */
        double *Mk = &(M[n*k]); /* Row {k} of {M} */
        /* Pick a random direction {v} in {R^{k+1}}, store it in row {k} of {M}: */
        rn_throw_dir(k+1, Mk);
        /* We now apply to rows {0..k-1} a reflection that takes {u_k} to {v}. */
        /* Compute the unit vector {s[0..k]} normal to the bisector of {u_k} and {v}: */
        double s2 = 0;
        for (int32_t j = 0; j <= k; j++)
          { double sj = Mk[j] - (j == k ? 1 : 0); 
            s[j] = sj; 
            s2 += sj*sj;
          }
        if (s2 != 0)
          { /* Reflection is not trivial. */
            double sm = sqrt(s2);
            for (int32_t j = 0; j <= k; j++) { s[j] /= sm; }
            /* Apply reflection along {s} to each row: */
            for (int32_t i = 0; i < k; i++)
              { double *Mi = &(M[i*n]);
                (void)rn_mirror(k+1, Mi, s, Mi);
                flip = 1 - flip;
              }
          }
      }
    /* Now flip the first row if needed to keep the determinant positive: */
    if (flip != 0) { for (int32_t j = 0; j < n; j++) { M[j] = -M[j]; } }
  }

void rmxn_throw_ortho_complement(int32_t n, int32_t p, double A[], int32_t q, double M[])
  {
    demand(p + q <= n, "invalid {p+q}");
    if (q == 0) { return; }
    for (int32_t k = 0; k < q; k++)
      { double *Mk = &(M[k*n]); /* Row {k} o f{M}. */
        while (TRUE)
          { /* Pick a random unit vector in {\RR^n}: */
            rn_throw_dir (n, Mk);
            /* Project it onto the space orthogonal to {A} and {M}: */
            for (int32_t r = 0; r < p; r++)
              { double *Ar = &(A[r*n]);
                double s = rn_dot(n, Ar, Mk);
                rn_mix_in(n, -s, Ar, Mk);
              }
            for (int32_t r = 0; r < k; r++)
              { double *Mr = &(M[r*n]);
                double s = rn_dot(n, Mr, Mk);
                rn_mix_in(n, -s, Mr, Mk);
              }
            /* Normalize and accept if not too small: */
            double len = rn_dir(n, Mk, Mk);
            if (len >= 1.0e-4) { break; }
          }
      }
  }

void rmxn_throw_directions(int32_t m, int32_t n, double U[])
  {
    for (int32_t i = 0; i < m; i++)
      { double *Ui = &(U[i*n]); 
        if (i < n)
          { rn_axis(n, i, Ui); }
        else
          { double minMaxCos = +INF;
            double u[n];
            int32_t maxTry = 10;
            int32_t nTry = 0;
            while(TRUE)
              { nTry++;
                /* Generate a random direction {u}: */
                rn_throw_dir (n, u);
                /* Check whether {u} is far enough from previous directions: */
                double maxCos = 0.0; /* Max {fabs(cos(u,Uk))} for previous rows {Uk}. */
                for (int32_t k = 0; k < i; k++)
                  { double *Uk = &(U[k*n]); 
                    double cos = fabs(rn_dot(n, u, Uk));
                    if (cos > maxCos) { maxCos = cos; }
                  }
                if (maxCos < minMaxCos)
                  { rn_copy(n, u, Ui); 
                    minMaxCos = maxCos;
                    if ((nTry >= maxTry) || (minMaxCos < cos(M_PI/6))) { break; }
                  }
              }
          }
      }
  }

void rmxn_spin_rows(int32_t m, int32_t n, double A[], double M[])
  { /* Generate a random orthonormal {n×n} matrix {N}: */
    double N[n*n];
    rmxn_throw_ortho(n, N);
    /* Map each row of {A} by {N} (beware of aliasing between {A} and {M}): */
    double v[n];
    for (int32_t i = 0; i < m; i++)
      { rmxn_map_row(n, n, &(A[i*n]), N, v);
        rn_copy(n, v, &(M[i*n]));
      }
  }

void rmxn_spin_cols(int32_t m, int32_t n, double A[], double M[])
  { /* Generate a random orthonormal {m×m} matrix {N}: */
    double N[m*m];
    rmxn_throw_ortho(m, N);
    /* Map each col of {A} by {N} (beware of aliasing between {A} and {M}): */
    double a[m], v[m];
    for (int32_t j = 0; j < n; j++)
      { rmxn_get_col(m, n, A, j, a);
        rmxn_map_col(m, m, a, N, v);
        rmxn_set_col(m, n, M, j, v);
      }
  }
  
void rmxn_shift_rows(int32_t m, int32_t n, double A[], double v[], double M[])
  { for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++) 
          { int32_t ij = i*n + j; M[ij] = A[ij] + v[j]; }
      }
  }
  
void rmxn_shift_cols(int32_t m, int32_t n, double v[], double A[], double M[])
  { for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++) 
          { int32_t ij = i*n + j; M[ij] = A[ij] + v[i]; }
      }
  }

