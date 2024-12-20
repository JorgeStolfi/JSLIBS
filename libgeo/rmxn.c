/* See rmxn.h. */
/* Last edited on 2024-12-01 00:18:14 by stolfi */

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <rn.h>
#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <affirm.h>
#include <gausol_triang.h>
#include <gausol_solve.h>

#include <rmxn.h>

void rmxn_zero(uint32_t m, uint32_t n, double *M)
  { int32_t t = 0;
    for (uint32_t i = 0;  i < m; i++)
      for (uint32_t j = 0;  j < n; j++)
        { M[t] = 0.0; t++; }
  }

void rmxn_copy(uint32_t m, uint32_t n, double *A, double *M)
  { uint32_t mn = m*n;
    for (uint32_t ij = 0;  ij < mn; ij++) { M[ij] = A[ij]; }
  }

void rmxn_get_row(uint32_t m, uint32_t n, double *A, uint32_t i, double *r)
  { double *Ai = &(A[i*n]);
    for (uint32_t j = 0;  j < n; j++) { r[j] = Ai[j]; }
  }
  
void rmxn_set_row(uint32_t m, uint32_t n, double *A, uint32_t i, double *r)
  { double *Ai = &(A[i*n]);
    for (uint32_t j = 0;  j < n; j++) { Ai[j] = r[j]; }
  }

void rmxn_get_col(uint32_t m, uint32_t n, double *A, uint32_t j, double *r)
  { double *Aij = &(A[j]);
    for (uint32_t i = 0;  i < m; i++) { r[i] = (*Aij); Aij += n; }
  }

void rmxn_set_col(uint32_t m, uint32_t n, double *A, uint32_t j, double *r)
  { double *Aij = &(A[j]);
    for (uint32_t i = 0;  i < m; i++) { (*Aij) = r[i]; Aij += n; }
  }

void rmxn_ident(uint32_t m, uint32_t n, double *M)
  { int32_t t = 0;
    for (uint32_t i = 0;  i < m; i++)
      for (uint32_t j = 0;  j < n; j++)
        { M[t] = (i == j ? 1.0 : 0.0); t++; }
  }

void rmxn_add(uint32_t m, uint32_t n, double *A, double *B, double *M)
  { rn_add(m*n, A, B, M); }
    
void rmxn_sub(uint32_t m, uint32_t n, double *A, double *B, double *M)
  { rn_sub(m*n, A, B, M); }

void rmxn_map_row (uint32_t m, uint32_t n, double *x, double *A, double *r)
  { for (uint32_t j = 0;  j < n; j++)
      { double sum = 0.0, corr = 0.0; 
        uint32_t t = j;
        for (uint32_t i = 0;  i < m; i++) 
          { double term = x[i] * A[t];
            /* Kahan's summation: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
            t += n;
          }
        r[j] = sum;
      }
  }

void rmxn_map_col (uint32_t m, uint32_t n, double *A, double *x, double *r)
  { int32_t t = 0;
    for (uint32_t i = 0;  i < m; i++)
      { double sum = 0.0, corr = 0.0; 
        for (uint32_t j = 0;  j < n; j++)
          { double term = A[t] * x[j]; 
            /* Kahan's summation: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
            t++;
          }
        r[i] = sum;
      }
  }

void rmxn_mul (uint32_t m, uint32_t p, uint32_t n, double *A, double *B, double *M)
  { uint32_t r = 0, v = 0;
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++)
          { double sum = 0.0, corr = 0.0;
            uint32_t t = j;
            for (uint32_t k = 0;  k < p; k++)
              { double term = A[r+k]*B[t];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
                t += n;
              }
            M[v] = sum; v++;
          }
        r += p;
      }
  }

void rmxn_mul_tr (uint32_t m, uint32_t n, uint32_t p, double *A, double *B, double *M)
  { uint32_t v = 0, r = 0;
    for (uint32_t i = 0;  i < m; i++)
      { uint32_t s = 0;
        for (uint32_t j = 0;  j < n; j++)
          { double sum = 0.0, corr = 0.0;
            for (uint32_t k = 0;  k < p; k++) 
              { double term = A[r+k]*B[s+k];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            M[v] = sum; v++;
            s += p;
          }
        r += p;
      }
  }

void rmxn_tr_mul (uint32_t p, uint32_t m, uint32_t n, double *A, double *B, double *M)
  { uint32_t v = 0;
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++)
          { double sum = 0.0, corr = 0.0;
            uint32_t r = i, s = j;
            for (uint32_t k = 0;  k < p; k++) 
              { double term = A[r]*B[s];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
                r += m; s += n;
              }
            M[v] = sum; v++;
          }
      }
  }

double rmxn_det (uint32_t n, double *A)
  { uint32_t n2 = n*n;
    double *C = rmxn_alloc(n,n);
    for (uint32_t t = 0;  t < n2; t++) { C[t] = A[t]; }
    /* !!! Using full pivoting. Would partial pivoting be better? !!! */
    uint32_t prow[n], pcol[n]; /* Permutation matrices for Gaussian eleimination. */
    double det_tri;
    uint32_t rank;
    double tiny = 1.0e-180;
    gausol_triang_reduce(n, prow, n, pcol, C, 0, NULL, tiny, &det_tri, &rank);
    assert(! isnan(det_tri));
    double det_A = (rank < n ? 0.0 : det_tri); 
    free(C);
    return det_A;
  }

double rmxn_det_by_enum(uint32_t m, uint32_t n, double A[], uint32_t q)
  { if ((q > m) || (q > n)) { return 0.0; }
    demand(q <= rmxn_det_by_enum_SIZE_MAX, "too many determinant terms to enumerate");
    
    uint32_t perm[q];

    for (uint32_t i = 0;  i < q; i++) { perm[i] = i; }

    double det = 0.0, corr = 0.0;

    auto void add_terms(uint32_t k, double product);
      /* Generates all permutations of {perm[k..q-1]}, and adds to
        {det} the partial {product} times {A[r,perm[r]} for {r} in {k..q-1}.
        Upon return, {perm[k..q-1]} are back to their initial
        state. */

    add_terms(0, 1.0);
    return det;

    void add_terms(uint32_t k, double product)
      { if (k >= q)
          { /* The {product} is a complete determinant term */
            double term = product; 
            /* Kahan's detmation: */
            double tcorr = term - corr;
            double newDet = det + tcorr;
            corr = (newDet - det) - tcorr;
            det = newDet;
          }
        else
          { for (uint32_t j = k; j < q; j++)
              { if (j != k) { uint32_t t = perm[k]; perm[k] = perm[j]; perm[j] = t; }
                uint32_t row = k, col = perm[k];
                double prod_more = product * (k == j ? 1.0 : -1.0) * A[row*n + col];
                add_terms(k+1, prod_more);
                if (j != k) { uint32_t t = perm[k]; perm[k] = perm[j]; perm[j] = t; }
              }
          }
      }
  }

double rmxn_inv (uint32_t n, double *A, double *M)
  { bool_t debug = FALSE;
    double *D = rmxn_alloc(n,n);
    /* Copy {A} into {C}, fill {D} with the identity: */
    for (uint32_t i = 0;  i < n; i++) 
      { double *Dij = &(D[n*i]); 
        for (uint32_t j = 0;  j < n; j++) { (*Dij) = (j == i ? 1.0 : 0.0); Dij++; }
      }
    uint32_t rank;
    double det_A;
    double pivot_rows = TRUE;
    double pivot_cols = TRUE;
    gausol_solve(n, n, A, n, D, M, pivot_rows, pivot_cols, 0.0, &det_A, &rank);
    if (debug) { fprintf(stderr, "n = %d  rank = %d  det = %24.16e\n", n, rank, det_A); }
    assert(isfinite(det_A));
    assert((det_A != 0.0) == (rank == n));
    free(D);
    return det_A;
  }

void rmxn_scale(uint32_t m, uint32_t n, double s, double *A, double *M) 
  { uint32_t k = 0;
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++)
          { M[k] = s*A[k]; k++; }
      }
  }

void rmxn_mix(uint32_t m, uint32_t n, double s, double *A, double t, double *B, double *M)
  { uint32_t k = 0;
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++)
          { M[k] = s*A[k] + t*B[k]; k++; }
      }
  }

void rmxn_rel_diff(uint32_t m, uint32_t n, double *A, double *B, double *M)
  { uint32_t k = 0;
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++)
          { M[k] = rel_diff(A[k], B[k]); k++; }
      }
  }

double rmxn_norm_sqr(uint32_t m, uint32_t n, double *A)
  { /* !!! Should use Kahan's summation !!! */
    double s = 0;
    uint32_t k = 0;
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++)
          { double Aij = A[k];
            s += Aij*Aij; k++; 
          }
      }
    return s;
  }

double rmxn_norm(uint32_t m, uint32_t n, double *A)
  { double s = rmxn_norm_sqr(m, n, A);
    return sqrt(s);
  }

double rmxn_normalize(uint32_t m, uint32_t n, double *A)
  { 
    double w = rmxn_norm(m, n, A);
    if (w != 0)
      { uint32_t k = 0;
        for (uint32_t i = 0;  i < m; i++)
          { for (uint32_t j = 0;  j < n; j++)
             { A[k] /= w; k++; }
          }
      }
    else
      { uint32_t k = 0;
        for (uint32_t i = 0;  i < m; i++)
          { for (uint32_t j = 0;  j < n; j++)
             { A[k] = NAN; k++; }
          }
      }
    return w;
  }

double rmxn_mod_norm_sqr(uint32_t n, double *A)
  {/* !!! Should use Kahan's summation !!! */
    double s = 0.0;
    uint32_t k = 0;
    for (uint32_t i = 0;  i < n; i++)
      for (uint32_t j = 0;  j < n; j++)
        { double Aij = A[k]; k++; 
          double Dij = Aij - (i == j ? 1 : 0);
          s += Dij*Dij;
        }
    return s; 
  }

double rmxn_max_abs_elem(uint32_t m, uint32_t n, double *A)
  { double emax = 0.0;
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++)
          { double Mij = fabs(A[i*n + j]);
            if (Mij > emax) { emax = Mij; }
          }
      }
    return emax;
  }

double rmxn_max_abs_elem_in_col(uint32_t m, uint32_t n, double M[], uint32_t j)
  { double emax = 0.0;
    for (uint32_t i = 0; i < m; i++)
      { double Mij = fabs(M[i*n + j]);
        if (Mij > emax) { emax = Mij; }
      }
    return emax;
  }

double rmxn_max_abs_elem_in_row(uint32_t m, uint32_t n, double M[], uint32_t i)
  { double emax = 0.0;
    for (uint32_t j = 0; j < n; j++)
      { double Mij = fabs(M[i*n + j]);
        if (Mij > emax) { emax = Mij; }
      }
    return emax;
  }

void rmxn_LT_inv_map_row(uint32_t n, double *y, double *L, double *r)
  { uint32_t j = n;
    while (j > 0)
      { j--; /* Safe since {j > 0}. */
        double sum = y[j], corr = 0.0;
        for (uint32_t i = j+1; i < n; i++)
          { double term = -r[i]*L[n*i + j];
            /* Kahan's summation: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
          }
        r[j] = sum/L[n*j + j];
      }
  }

void rmxn_LT_inv_map_col(uint32_t m, double *L, double *y, double *r)
  { for (uint32_t i = 0;  i < m; i++)
      { double sum = y[i], corr = 0.0;
        for (uint32_t j = 0;  j < i; j++)
          { double term = -L[m*i + j]*r[j];
            /* Kahan's summation: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
          }
        r[i] = sum/L[m*i + i];
      }
  }

void rmxn_LT_pos_div(uint32_t m, uint32_t n, double *A, double *L, double *M)
  { for (uint32_t k = 0;  k < m; k++)
      { uint32_t j = n;
        while (j > 0)
          { j--; /* Safe here since {j > 0}. */
            double sum = A[n*k + j], corr = 0.0;
            for (uint32_t i = j+1; i < n; i++)
              { double term = -M[n*k + i]*L[n*i + j];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            M[n*k + j] = sum/L[n*j + j];
          }
      }
  }

void rmxn_LT_pre_div(uint32_t m, uint32_t n, double *L, double *A, double *M)
  { for (uint32_t k = 0;  k < n; k++)
      { for (uint32_t i = 0;  i < m; i++)
          { double sum = A[n*i + k], corr = 0.0;
            for (uint32_t j = 0;  j < i; j++)
              { double term = -L[m*i + j]*M[n*j + k];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            M[n*i + k] = sum/L[m*i + i];
          }
      }
  }

void rmxn_cholesky(uint32_t n, double *A, double *L)
  { 
    /* Andre-Louis Cholesky (spelled with a 'y'), born in France in
      1875, was a geodesist in the French military. He developed his
      computational procedure to solve geodetic problems, among which
      was, in the words of his eulogy, the "problem of levelling in
      Morocco." He died in battle in 1918. Benoit published the
      computational method in "Bulletin geodesique" in 1924. 
        -- David Pattison, Washington, DC

      Although the name is originally Polish or Russian, Cholesky
      probably pronounced it himself in the French fashion -- namely,
      stressed on the final syllable, and with an "sh" sound, as in
      "Chopin" and "Chostakovich" (and not with a "tsh" sound, which in
      French would be spelled "tch", as in "Tchaikovski" or
      "Tchekov"). */
  
    for (uint32_t i = 0;  i < n; i++)
      { /* Compute {L[i,j] = X[i,j] - SUM{ L[i,k]*A[j,k] : k = 0..j-1 })/L[j,j]} */
        /* where X[i,j] = (j <= i ? A[i,j] : 0) */
        for (uint32_t j = 0;  j < i; j++)
          { double sum = 0.0, corr = 0.0;
            for (uint32_t k = 0;  k < j; k++)
              { double term = L[n*i + k]*L[n*j + k];
                /* Kahan's summation formula: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            double Ljj = L[n*j + j];
            affirm(Ljj != 0.0, "zero element in diagonal");
            L[n*i + j] = (A[n*i + j] - sum)/Ljj;
          }
        
        /* Compute {L[i,i] = sqrt(A[i,i] - SUM{ L[i,j]^2 : j = 0..i-1 })} */
        double sum = 0.0, corr = 0.0;
        for (uint32_t j = 0;  j < i; j++)
          { double w = L[n*i+j], term = w*w; 
            /* Kahan's summation formula: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
          }
        double Lii = A[n*i+i] - sum;
        affirm (Lii >= 0.0, "matrix is not positive definite?");
        L[n*i + i] = sqrt(Lii);
        for (uint32_t j = i+1; j < n; j++) { L[n*i + j] = 0.0; }
      }
  }

void rmxn_perturb_unif(uint32_t m, uint32_t n, double pabs, double prel, double *A)
  {
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++) 
          { double *Aij = &(A[i*n + j]);
            double mag = pabs + prel*fabs(*Aij);
            double d = dabrandom(-mag, +mag);
            (*Aij) += d; 
          }
      }
  }
  
void rmxn_cleanup(uint32_t m, uint32_t n, double *A, double tiny)
  {
    if (tiny <= 0.0) { return; }
    for (uint32_t i = 0; i < m; i++)
      { for (uint32_t j = 0; j < n; j++) 
          { double *Aij = &(A[i*n + j]);
            if (fabs(*Aij) < tiny) { (*Aij) = 0.0; }
          }
      }
  }

void rmxn_print (FILE *f, uint32_t m, uint32_t n, double *A)
  { rmxn_gen_print(f, m, n, A, NULL, NULL, NULL, NULL, NULL, NULL, NULL); }

void rmxn_gen_print
  ( FILE *f, 
    uint32_t m, 
    uint32_t n, double *A,
    char *fmt, 
    char *olp, char *osep, char *orp,
    char *ilp, char *isep, char *irp
  )
  {
    rmxn_gen_print3 
      ( f, m, 
        n,  A, 
        0, NULL, 
        0, NULL, 
        fmt,
        olp, osep, orp,
        ilp, isep, irp,
        NULL
      );
  }  

void rmxn_gen_print2 
  ( FILE *f, 
    uint32_t m,
    uint32_t n1, double *A1,
    uint32_t n2, double *A2,
    char *fmt, 
    char *olp, char *osep, char *orp,     /* Outer delimiters. */
    char *ilp, char *isep, char *irp,     /* Inner delimiters. */
    char *msep                            /* Matrix separator. */
  )
  {
    rmxn_gen_print3 
      ( f, m, 
        n1, A1, 
        n2, A2, 
        0, NULL, 
        fmt,
        olp, osep, orp,
        ilp, isep, irp,
        msep
      );
  }

void rmxn_gen_print3 
  ( FILE *f, 
    uint32_t m,
    uint32_t n1, double *A1,
    uint32_t n2, double *A2,
    uint32_t n3, double *A3,
    char *fmt, 
    char *olp, char *osep, char *orp,     /* Outer delimiters. */
    char *ilp, char *isep, char *irp,     /* Inner delimiters. */
    char *msep                            /* Matrix separator. */
  )
  {
    if (olp == NULL) { olp = "(\n"; }
    if (osep == NULL) { osep = "\n"; }
    if (orp == NULL) { orp = "\n)"; }
    if (ilp == NULL) { ilp = "  ("; }
    if (isep == NULL) { isep = " "; }
    if (msep == NULL) { msep = "  "; }
    if (irp == NULL) { irp = ")"; }
    if (fmt == NULL) { fmt = "%16.8e"; }
    
    fputs(olp, f);
    for (uint32_t i = 0;  i < m; i++)
      { if (i > 0) { fputs(osep, f); }
        for (uint32_t k = 1;  k <= 3; k++)
          { uint32_t nk = (k == 1 ? n1 : (k == 2 ? n2 : n3));
            double *Ak = (k == 1 ? A1 : (k == 2 ? A2 : A3));
            if (Ak != NULL)
              { double *Aki = (nk == 0 ? NULL : &(Ak[i*nk]));
                if (k > 1) { fputs(msep, f); }
                fputs(ilp, f);
                for (uint32_t j = 0;  j < nk; j++) 
                  { if (j > 0) { fputs(isep, f); }
                    fprintf(f, fmt, Aki[j]);
                  }
                fputs(irp, f);
              }
          }
      }
    fputs(orp, f);
    fflush(f);
  }  

double *rmxn_alloc(uint32_t m, uint32_t n)
  { double *p = talloc(m*n, double);
    return p;
  }
