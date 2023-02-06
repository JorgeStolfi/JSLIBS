/* See rmxn.h. */
/* Last edited on 2023-02-03 05:40:23 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <rn.h>
#include <bool.h>
#include <jsmath.h>
#include <affirm.h>
#include <gauss_elim.h>

#include <rmxn.h>

void rmxn_zero(int32_t m, int32_t n, double *M)
  { int32_t i, j, t;
    t = 0;
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        { M[t] = 0.0; t++; }
  }

void rmxn_copy(int32_t m, int32_t n, double *A, double *M)
  { int32_t mn = m*n;
    int32_t ij;
    for (ij = 0; ij < mn; ij++) { M[ij] = A[ij]; }
  }

void rmxn_get_row(int32_t m, int32_t n, double *A, int32_t i, double *r)
  { double *Ai = &(A[i*n]);
    int32_t j;
    for (j = 0; j < n; j++) { r[j] = Ai[j]; }
  }
  
void rmxn_set_row(int32_t m, int32_t n, double *A, int32_t i, double *r)
  { double *Ai = &(A[i*n]);
    int32_t j;
    for (j = 0; j < n; j++) { Ai[j] = r[j]; }
  }

void rmxn_get_col(int32_t m, int32_t n, double *A, int32_t j, double *r)
  { double *Aij = &(A[j]);
    int32_t i;
    for (i = 0; i < m; i++) { r[i] = (*Aij); Aij += n; }
  }

void rmxn_set_col(int32_t m, int32_t n, double *A, int32_t j, double *r)
  { double *Aij = &(A[j]);
    int32_t i;
    for (i = 0; i < m; i++) { (*Aij) = r[i]; Aij += n; }
  }

void rmxn_ident(int32_t m, int32_t n, double *M)
  { int32_t i, j, t;
    t = 0;
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        { M[t] = (i == j ? 1.0 : 0.0); t++; }
  }

void rmxn_map_row (int32_t m, int32_t n, double *x, double *A, double *r)
  { int32_t i, j;
    for (j = 0; j < n; j++)
      { double sum = 0.0, corr = 0.0; 
        int32_t t = j;
        for (i = 0; i < m; i++) 
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

void rmxn_map_col (int32_t m, int32_t n, double *A, double *x, double *r)
  { int32_t i, j, t = 0;
    for (i = 0; i < m; i++)
      { double sum = 0.0, corr = 0.0; 
        for (j = 0; j < n; j++)
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

void rmxn_mul (int32_t m, int32_t p, int32_t n, double *A, double *B, double *M)
  { int32_t i, j, k, r = 0, v = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { double sum = 0.0, corr = 0.0;
            int32_t t = j;
            for (k = 0; k < p; k++)
              { double term = A[r+k]*B[t];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
                t+= n;
              }
            M[v] = sum; v++;
          }
        r += p;
      }
  }

void rmxn_mul_tr (int32_t m, int32_t n, int32_t p, double *A, double *B, double *M)
  { int32_t i, j, k;
    int32_t v = 0;
    int32_t r = 0;
    for (i = 0; i < m; i++)
      { int32_t s = 0;
        for (j = 0; j < n; j++)
          { double sum = 0.0, corr = 0.0;
            for (k = 0; k < p; k++) 
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

void rmxn_tr_mul (int32_t p, int32_t m, int32_t n, double *A, double *B, double *M)
  { int32_t i, j, k;
    int32_t v = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { double sum = 0.0, corr = 0.0;
            int32_t r = i, s = j;
            for (k = 0; k < p; k++) 
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

double rmxn_det (int32_t n, double *A)
  { int32_t n2 = n*n, t;
    double *C = (double *)notnull(malloc(n2*sizeof(double)), "no mem for C");
    double det = 1.0;
    for (t = 0; t < n2; t++) { C[t] = A[t]; }
    gsel_triangularize(n, n, C, FALSE, 0.0);
    for (t = 0; t < n2; t += n+1) { det *= C[t]; }
    free(C);
    return det;
  }

double rmxn_inv (int32_t n, double *A, double *M)
  { int32_t i, j;
    int32_t nC = 2*n;
    int32_t nA = n;
    int32_t nM = n;
    double *C = (double *)notnull(malloc(n*nC*sizeof(double)), "no mem for C");
    /* Copy {A} into the left half of {C}, fill the right half with the identity: */
    for (i = 0; i < n; i++) 
      { double *Cij = &(C[nC*i]); 
        double *Aij = &(A[nA*i]);
        for (j = 0; j < nA; j++) { (*Cij) = (*Aij); Cij++; Aij++; }
        for (j = 0; j < nM; j++) { (*Cij) = (j == i ? 1.0 : 0.0); Cij++; }
      }
    gsel_triangularize(n, nC, C, FALSE, 0.0);
    double det = 1.0;
    for (i = 0; i < n; i++) { det *= C[nC*i + i]; }
    gsel_diagonalize(n, nC, C, 0);
    gsel_normalize(n, nC, C, 0);
    for (i = 0; i < n; i++) 
      { double *Cij = &(C[nC*i + n]); 
        double *Mij = &(M[nM*i]);
        for (j = 0; j < nM; j++) { (*Mij) = (*Cij); Mij++; Cij++; }
      }
    free(C);
    return det;
  }
  
double rmxn_inv_full (int32_t n, double *A, double *M)
  { 
    int32_t i, j, t;
    /* Copy {A} into {M}: */
    t = 0;
    for (i = 0; i < n; i++) { for (j = 0; j < n; j++) { M[t] = A[t]; t++; } }
 
    double det = 1.0; /* Accumulates the determinant. */
 
    int32_t prow[n], pcol[n]; 
    int32_t k = 0;
    
    while (k < n)
      { /* Process the remaining {n-k}x{n-k} submatrix starting at {A[k][k]}. */

        /* Find the element with largest abs value in the submatrix. */
        /* Store its value in {biga}, its indices in {prow[k],pcol[k]}: */
        prow[k] = k;
        pcol[k] = k;
        double biga = M[k*n+k];
        for (i = k; i < n; i++)
          { for (j = k; j < n; j++)
              { if (fabs (M[i*n+j]) > fabs (biga))
                  { biga = M[i*n+j]; prow[k] = i; pcol[k] = j; }
              }
          }
 
        /* Accumulate the determinant: */
        det *= biga; /* Product of pivots */

        /* If the matrix is all zeros, we break here: */
        if (biga == 0.0) { break; }
        
        /* Swap rows {k} and {prow[k]}: */
        i = prow[k];
        if (i > k)
          { for (j = 0; j < n; j++)
              { double hold = -M[k*n+j];
                M[k*n+j] = M[i*n+j];
                M[i*n+j] = hold;
              }
          }
 
        /* Swap columns {k} and {pcol[k]}: */
        j = pcol[k];
        if (j > k)
          { for (i = 0; i < n; i++)
              { double hold = -M[i*n+k];
                M[i*n+k] = M[i*n+j];
                M[i*n+j] = hold;
              }
          }
 
        /* Negate column {k} and divide by minus the pivot {biga}. */
        for (i = 0; i < n; i++) { if (i != k) { M[i*n+k] /= -biga; } }
 
        /* Reduce the matrix: */
        for (i = 0; i < n; i++)
          { if (i != k)
              { double hold = M[i*n+k];
                for (j = 0; j < n; j++)
                  { if (j != k) { M[i*n+j] += hold * M[k*n+j]; } }
              }
          }
 
        /* Divide row {k} by pivot {biga}: */
        for (j = 0; j < n; j++) { if (j != k) { M[k*n+j] /= biga; } }
 
        /* Set the diagonal element to its reciprocal: */
        M[k*n+k] = 1.0 / biga;
         
        /* Shrink: */
        k++;
      }
 
    /* Undo the row and column interchanges, transposed: */
    while(k > 0)
      { 
        /* Back up in {k}: */
        k--;
      
        /* Swap column {k} with column {prow[k]}: */
        i = prow[k];
        if (i > k)
          { for (j = 0; j < n; j++) 
              { double hold = M[j*n+k];
                M[j*n+k] = -M[j*n+i];
                M[j*n+i] = hold;
              }
          }
 
        /* Swap row {k} with row {pcol[k]}: */
        j = pcol[k];
        if (j > k)
          { for (i = 0; i < n; i++)
              { double hold = M[k*n+i];
                M[k*n+i] = -M[j*n+i];
                M[j*n+i] = hold;
              }
          }
      }

    /* Return the determinant: */
    return det;
  }


void rmxn_scale(int32_t m, int32_t n, double s, double *A, double *M) 
  { int32_t i, j, k = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { M[k] = s*A[k]; k++; }
      }
  }

void rmxn_mix(int32_t m, int32_t n, double s, double *A, double t, double *B, double *M)
  { int32_t i, j, k = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { M[k] = s*A[k] + t*B[k]; k++; }
      }
  }

void rmxn_rel_diff(int32_t m, int32_t n, double *A, double *B, double *M)
  { int32_t i, j, k = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { M[k] = rel_diff(A[k], B[k]); k++; }
      }
  }

double rmxn_norm_sqr(int32_t m, int32_t n, double *A)
  { double s = 0;
    int32_t i, j, k = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { double Aij = A[k]; k++; 
            s += Aij*Aij;
          }
      }
    return s;
  }

double rmxn_norm(int32_t m, int32_t n, double *A)
  { double s = 0;
    int32_t i, j, k = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { double Aij = A[k]; k++; 
            s += Aij*Aij;
          }
      }
    return sqrt(s);
  }

double rmxn_mod_norm_sqr(int32_t n, double *A)
  {
    double s = 0.0;
    int32_t i, j;
    int32_t k = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
        { double Aij = A[k]; k++; 
          double Dij = Aij - (i == j ? 1 : 0);
          s += Dij*Dij;
        }
    return s; 
  }

double rmxn_max_abs_elem(int32_t m, int32_t n, double *A)
  { int32_t i, j;
    double emax = 0.0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { double Mij = fabs(A[i*n + j]);
            if (Mij > emax) { emax = Mij; }
          }
      }
    return emax;
  }

void rmxn_LT_inv_map_row(int32_t n, double *y, double *L, double *r)
  { int32_t i, j;
    for (j = n-1; j >= 0; j--)
      { double sum = y[j], corr = 0.0;
        for (i = j+1; i < n; i++)
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

void rmxn_LT_inv_map_col(int32_t m, double *L, double *y, double *r)
  { int32_t i, j;
    for (i = 0; i < m; i++)
      { double sum = y[i], corr = 0.0;
        for (j = 0; j < i; j++)
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

void rmxn_LT_pos_div(int32_t m, int32_t n, double *A, double *L, double *M)
  { int32_t i, j, k;
    for (k = 0; k < m; k++)
      { for (j = n-1; j >= 0; j--)
          { double sum = A[n*k + j], corr = 0.0;
            for (i = j+1; i < n; i++)
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

void rmxn_LT_pre_div(int32_t m, int32_t n, double *L, double *A, double *M)
  { int32_t i, j, k;
    for (k = 0; k < n; k++)
      { for (i = 0; i < m; i++)
          { double sum = A[n*i + k], corr = 0.0;
            for (j = 0; j < i; j++)
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

void rmxn_cholesky(int32_t n, double *A, double *L)
  { 
    /* Andre-Louis Cholesky (spelled with a 'y'), born in France in
      1875, was a geodesist in the French military. He developed his
      computational procedure to solve geodetic problems, among which
      was, in the words of his eulogy, the "problem of levelling in
      Morocco." He died in battle in 1918. Benoit published the
      computational method in "Bulletin geodesique" in 1924. 
        -- David Pattison, Washington, DC

      Although the name is originally Polish or Russian, Cholesky 
      probably pronounced it himself in the French fashion, i.e.
      with an "sh" sound, as in "Chopin" and "Chostakovich"; and not
      with a "tsh" sound --- which in French would be spelled "Tch",
      as in "Tchaikovski" or "Tchekov". */
  
    int32_t i, j, k;

    for (i = 0; i < n; i++)
      { double Lii;
        /* Compute {L[i,j] = X[i,j] - SUM{ L[i,k]*A[j,k] : k = 0..j-1 })/L[j,j]} */
        /* where X[i,j] = (j <= i ? A[i,j] : 0) */
        for (j = 0; j < i; j++)
          { double sum = 0.0, corr = 0.0;
            for (k = 0; k < j; k++)
              { double term = L[n*i + k]*L[n*j + k];
                /* Kahan's summation formula: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            { double Ljj = L[n*j + j];
              affirm(Ljj != 0.0, "zero element in diagonal");
              L[n*i + j] = (A[n*i + j] - sum)/Ljj;
            } 
          }
        
        /* Compute {L[i,i] = sqrt(A[i,i] - SUM{ L[i,j]^2 : j = 0..i-1 })} */
        { double sum = 0.0, corr = 0.0;
          for (j = 0; j < i; j++)
            { double w = L[n*i+j], term = w*w; 
              /* Kahan's summation formula: */
              double tcorr = term - corr;
              double newSum = sum + tcorr;
              corr = (newSum - sum) - tcorr;
              sum = newSum;
            }
          Lii = A[n*i+i] - sum;
        }
        affirm (Lii >= 0.0, "matrix is not positive definite?");
        L[n*i + i] = sqrt(Lii);
        for (j = i+1; j < n; j++) { L[n*i + j] = 0.0; }
      }
  }

void rmxn_print (FILE *f, int32_t m, int32_t n, double *A)
  { rmxn_gen_print(f, m, n, A, NULL, NULL, NULL, NULL, NULL, NULL, NULL); }

void rmxn_gen_print
  ( FILE *f, int32_t m, int32_t n, double *A,
    char *fmt, 
    char *olp, char *osep, char *orp,
    char *ilp, char *isep, char *irp
  )
  {
    rmxn_gen_print3 
      ( f, m, 
        n, A, -1, NULL, -1, NULL, 
        fmt,
        olp, osep, orp,
        ilp, isep, irp,
        NULL
      );
  }  

void rmxn_gen_print2 
  ( FILE *f, int32_t m,
    int32_t n1, double *A1,
    int32_t n2, double *A2,
    char *fmt, 
    char *olp, char *osep, char *orp,     /* Outer delimiters. */
    char *ilp, char *isep, char *irp,     /* Inner delimiters. */
    char *msep                            /* Matrix separator. */
  )
  {
    rmxn_gen_print3 
      ( f, m, 
        n1, A1, n2, A2, -1, NULL, 
        fmt,
        olp, osep, orp,
        ilp, isep, irp,
        msep
      );
  }

void rmxn_gen_print3 
  ( FILE *f, int32_t m,
    int32_t n1, double *A1,
    int32_t n2, double *A2,
    int32_t n3, double *A3,
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
    for (int32_t i = 0; i < m; i++)
      { if (i > 0) { fputs(osep, f); }
        for (int32_t k = 1; k <= 3; k++)
          { int32_t nk = (k == 1 ? n1 : (k == 2 ? n2 : n3));
            double *Ak = (k == 1 ? A1 : (k == 2 ? A2 : A3));
            if (nk >= 0)
              { double *Aki = (nk == 0 ? NULL : &(Ak[i*nk]));
                if (k > 1) { fputs(msep, f); }
                fputs(ilp, f);
                for (int32_t j = 0; j < nk; j++) 
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

double *rmxn_alloc(int32_t m, int32_t n)
  { void *p = malloc(m*n*sizeof(double));
    affirm(p != NULL, "no memory for rmxn_t");
    return (double *)p;
  }

