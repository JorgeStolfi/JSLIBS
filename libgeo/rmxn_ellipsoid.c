/* See {minn_constr.h}. */
/* Last edited on 2024-12-05 11:48:30 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
 
#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <sym_eigen.h>
#include <affirm.h>
#include <rn.h>
#include <rmxn.h>
#include <rmxn_throw.h>

#include <rmxn_ellipsoid.h>

void rmxn_ellipsoid_pierce
  ( uint32_t n,
    double arad[],
    uint32_t d, 
    double S[],
    double U[], 
    double urad[]
  )
  {
    bool_t debug = FALSE;

    if (debug) { fprintf(stderr, "piercing the ellipsoid...\n"); }

    demand(n >= 0, "invalid {n}");
    demand((d >= 0) && (d <= n), "invalid {d}");
    if ((d > 0) && (n > 0))
      { 
        if (debug) { fprintf(stderr, "  computing the metric matrix {M} for {\\RF} ...\n"); }
        double M[d*d];
        for (uint32_t r = 0; r < d; r++)
          { double *Sr = &(S[r*n]);
            for (uint32_t s = 0;  s <= r; s++)
              { double *Ss = &(S[s*n]);
                double sum = 0.0;
                for (uint32_t i = 0;  i < n; i++)
                  { double Sri = Sr[i];
                    double Ssi = Ss[i];
                    double ri = arad[i];
                    if (ri > 0)
                      { sum += Sri*Ssi/(ri*ri); }
                    else
                      { demand (fabs(Sri) + fabs(Ssi) < 1.0e-12, "bad constraints"); }
                  }
                M[r*d + s] = sum;
                M[s*d + r] = sum; /* Diag i
                s assigned twice, but OK. */
              }
          }
        if (debug) { fprintf(stderr, "  computing the eigen decomp {S,d} of {M} ...\n"); }
        double Q[d*d]; /* Eigenvector matrix. */
        double e[d]; /* Eigenvalues. */
        uint32_t nev;
        sym_eigen(d, M, e, Q, &nev);
        demand(nev == d, "failed to determine eigenvalues of {M}");
        /* Check that all eigenvalues are positive and
          compute the search radii {urad[0..d-1]} of {\RF}: */
        for (uint32_t k = 0;  k < d; k++)
          { demand(e[k] > 0.0, "non-positive eigenvalue");
            urad[k] = 1.0/sqrt(e[k]);
          }

        if (debug) { fprintf(stderr, "  computing the basis {U = Q S} aligned with axes of {\\RF} ...\n"); }
        for (uint32_t k = 0;  k < d; k++)
          { double *Qk = &(Q[k*d]);
            double *Uk = &(U[k*n]);
            for (uint32_t i = 0;  i < n; i++)
              { double sum = 0.0;
                for (uint32_t s = 0;  s < d; s++)
                  { double *Ss = &(S[s*n]);
                    double Qrs = Qk[s];
                    double Ssi = Ss[i];
                    sum += Qrs*Ssi;
                  }
                Uk[i] = sum;
              }
            if (debug) 
              { fprintf(stderr, "U[%d,*] = ", k);
                rn_print(stderr, n, Uk);
                fprintf(stderr, " rad = %12.7f\n", urad[k]);
              }
          }
      }
    if (debug) { fprintf(stderr, "done piercing the ellipsoid.\n"); }
  }

void rmxn_ellipsoid_cut
  ( uint32_t n,
    double arad[],
    uint32_t m, 
    double C[],
    uint32_t d,
    double U[], 
    double urad[]
  )
  { 
    bool_t debug = FALSE;
    
    if (debug) { fprintf(stderr, "cutting the ellipsoid...\n"); }
    demand(m >= 0, "invalid {m}");
    demand(d >= 0, "invalid {d}");
    demand(n == m+d, "invalid {m,d,n}");
    if ((d > 0) && (n > 0))
      { 
        if (debug) { fprintf(stderr, "  computing the complementary orthogonal basis {S}...\n"); }
        double S[d*n];
        rmxn_throw_ortho_complement(n, m, C, d, S);
        rmxn_ellipsoid_pierce(n, arad, d, S, U, urad);
      }
    if (debug) { fprintf(stderr, "done cutting the ellipsoid.\n"); }
  }

void rmxn_ellipsoid_normalize_constraints
  ( uint32_t n,
    double arad[],
    uint32_t q,
    double A[],
    bool_t verbose, 
    uint32_t *m_P,
    double **C_P
  )
  { 
    bool_t debug = FALSE;
    
    if (debug) { fprintf(stderr, "  normalizing the constraints...\n"); }
    
    double *C = rmxn_alloc(n,n); /* The orthonormalized constraints. */
    uint32_t m = 0; /* Number of independent constraints found. */
    
    double v[n]; /* Work vector. */

    auto void add_constraint(double a[]);
      /* Makes the vector {a[0..n-1]} orthogonal to rows {0..m-1} of {C}.
        If the result is sufficiently distinct from the zero vector,
        appends it as row {m} of {C}, and increments {C}. */
    
    /* Add the given constraints: */    
    for (uint32_t i = 0;  i < q; i++)
      { if (verbose) { fprintf(stderr, "    adding given constraint {v*A[%d,*]' == 0} ...", i); }
        double *Ai = &(A[i*n]);
        add_constraint(Ai);
      }

    /* Add the constraints from zero base radii: */
    double e[n]; /* Canonical basis vector. */
    for (int32_t j = 0;  j < n; j++) { e[j] = 0.0; }
    for (int32_t i = 0;  i < n; i++)
      { demand(arad[i] >= 0, "invalid {arad}");
        if (arad[i] < 1.0e-100)
          { if (verbose) { fprintf(stderr, "    adding axial constraint {v[%d] == 0} ...", i); }
            e[i] = 1.0;
            add_constraint(e);
            e[i] = 0.0; /* Restore {e} to all zeros. */
          }
      }
    if (verbose) { fprintf(stderr, "  got %d non-redundant constraints\n", m); }
    if (m < n) { C = realloc(C, m*n); }
    (*m_P) = m;
    (*C_P) = C;
    if (debug) { fprintf(stderr, "  done normalizing the constraints.\n"); }
    return;
     
     /* INTERNAL IMPLS */
     
     void add_constraint(double a[])
       { double am = rn_norm(n, a);
         if (am > 1.0e-100)
           { rn_copy(n, a, v);
             for (uint32_t k = 0;  k < m; k++)
               { double *Ck = &(C[k*n]);
                 double d = rn_dot(n, v, Ck);
                 if (fabs(d) > 1.0e-200)
                   { rn_mix_in(n, -d, Ck, v); }
               }
             double vm = rn_dir(n, v, v);
             if (vm > 1.0e-8*am)
               { affirm(m < n, "orthogonalization failure"); 
                 if (verbose) { fprintf(stderr, " accepted as constraint %d\n", m); }
                 double *Cp = &(C[m*n]);
                 rn_copy(n, v, Cp);
                 m++;
               }
             else
               { if (verbose) { fprintf(stderr, " seems redundant\n"); } }
           }
         else
           { if (verbose) { fprintf(stderr, " null\n"); } }
       }
  }

