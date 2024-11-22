/* See rmxn.h. */
/* Last edited on 2024-11-22 05:33:18 by stolfi */

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
#include <gauss_elim.h>

#include <rmxn.h>

void rmxn_zero(uint32_t m, uint32_t n, double *M)
  { int32_t t = 0;
    for (int32_t i = 0; i < m; i++)
      for (int32_t j = 0; j < n; j++)
        { M[t] = 0.0; t++; }
  }

void rmxn_copy(uint32_t m, uint32_t n, double *A, double *M)
  { uint32_t mn = m*n;
    for (int32_t ij = 0; ij < mn; ij++) { M[ij] = A[ij]; }
  }

void rmxn_get_row(uint32_t m, uint32_t n, double *A, uint32_t i, double *r)
  { double *Ai = &(A[i*n]);
    for (int32_t j = 0; j < n; j++) { r[j] = Ai[j]; }
  }
  
void rmxn_set_row(uint32_t m, uint32_t n, double *A, uint32_t i, double *r)
  { double *Ai = &(A[i*n]);
    for (int32_t j = 0; j < n; j++) { Ai[j] = r[j]; }
  }

void rmxn_get_col(uint32_t m, uint32_t n, double *A, uint32_t j, double *r)
  { double *Aij = &(A[j]);
    for (int32_t i = 0; i < m; i++) { r[i] = (*Aij); Aij += n; }
  }

void rmxn_set_col(uint32_t m, uint32_t n, double *A, uint32_t j, double *r)
  { double *Aij = &(A[j]);
    for (int32_t i = 0; i < m; i++) { (*Aij) = r[i]; Aij += n; }
  }

void rmxn_ident(uint32_t m, uint32_t n, double *M)
  { int32_t t = 0;
    for (int32_t i = 0; i < m; i++)
      for (int32_t j = 0; j < n; j++)
        { M[t] = (i == j ? 1.0 : 0.0); t++; }
  }

void rmxn_add(uint32_t m, uint32_t n, double *A, double *B, double *M)
  { rn_add(m*n, A, B, M); }
    
void rmxn_sub(uint32_t m, uint32_t n, double *A, double *B, double *M)
  { rn_sub(m*n, A, B, M); }

void rmxn_map_row (uint32_t m, uint32_t n, double *x, double *A, double *r)
  { for (int32_t j = 0; j < n; j++)
      { double sum = 0.0, corr = 0.0; 
        int32_t t = j;
        for (int32_t i = 0; i < m; i++) 
          { double term = x[i] * A[t];
            /* Kahan's summation: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
            t += (int32_t)n;
          }
        r[j] = sum;
      }
  }

void rmxn_map_col (uint32_t m, uint32_t n, double *A, double *x, double *r)
  { int32_t t = 0;
    for (int32_t i = 0; i < m; i++)
      { double sum = 0.0, corr = 0.0; 
        for (int32_t j = 0; j < n; j++)
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
  { int32_t r = 0, v = 0;
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { double sum = 0.0, corr = 0.0;
            int32_t t = j;
            for (int32_t k = 0; k < p; k++)
              { double term = A[r+k]*B[t];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
                t += (int32_t)n;
              }
            M[v] = sum; v++;
          }
        r += (int32_t)p;
      }
  }

void rmxn_mul_tr (uint32_t m, uint32_t n, uint32_t p, double *A, double *B, double *M)
  { int32_t v = 0, r = 0;
    for (int32_t i = 0; i < m; i++)
      { int32_t s = 0;
        for (int32_t j = 0; j < n; j++)
          { double sum = 0.0, corr = 0.0;
            for (int32_t k = 0; k < p; k++) 
              { double term = A[r+k]*B[s+k];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            M[v] = sum; v++;
            s += (int32_t)p;
          }
        r += (int32_t)p;
      }
  }

void rmxn_tr_mul (uint32_t p, uint32_t m, uint32_t n, double *A, double *B, double *M)
  { int32_t v = 0;
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { double sum = 0.0, corr = 0.0;
            int32_t r = i, s = j;
            for (int32_t k = 0; k < p; k++) 
              { double term = A[r]*B[s];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
                r += (int32_t)m; s += (int32_t)n;
              }
            M[v] = sum; v++;
          }
      }
  }

double rmxn_det (uint32_t n, double *A)
  { uint32_t n2 = n*n;
    double *C = rmxn_alloc(n,n);
    double det = 1.0;
    for (int32_t t = 0; t < n2; t++) { C[t] = A[t]; }
    gsel_triangularize(n, n, C, FALSE, 0.0);
    for (int32_t t = 0; t < n2; t += (int32_t)n+1) { det *= C[t]; }
    free(C);
    return det;
  }

double rmxn_inv (uint32_t n, double *A, double *M)
  { uint32_t nC = 2*n;
    uint32_t nA = n;
    uint32_t nM = n;
    double *C = rmxn_alloc(n,nC);
    /* Copy {A} into the left half of {C}, fill the right half with the identity: */
    for (int32_t i = 0; i < n; i++) 
      { double *Cij = &(C[(int32_t)nC*i]); 
        double *Aij = &(A[(int32_t)nA*i]);
        for (int32_t j = 0; j < nA; j++) { (*Cij) = (*Aij); Cij++; Aij++; }
        for (int32_t j = 0; j < nM; j++) { (*Cij) = (j == i ? 1.0 : 0.0); Cij++; }
      }
    gsel_triangularize(n, nC, C, TRUE, 0.0);
    double det = 1.0;
    for (int32_t i = 0; i < n; i++) { det *= C[(int32_t)nC*i + i]; }
    if (det == 0.0)
      { /* Singular matrix: */
        double *Mij = M;
        for (int32_t i = 0; i < n; i++) 
          { for (int32_t j = 0; j < n; j++)
              { (*Mij) = NAN; Mij++; }
          }
      }
    else
      { gsel_diagonalize(n, nC, C);
        gsel_normalize(n, nC, C);
        for (int32_t i = 0; i < n; i++) 
          { double *Cij = &(C[(int32_t)nC*i + (int32_t)n]); 
            double *Mij = &(M[(int32_t)nM*i]);
            for (int32_t j = 0; j < nM; j++) { (*Mij) = (*Cij); Mij++; Cij++; }
          }
      }
    free(C);
    return det;
  }
  
double rmxn_inv_full (uint32_t n, double *A, double *M)
  { 
    /* Copy {A} into {M}: */
    int32_t t = 0;
    for (int32_t i = 0; i < n; i++) 
      { for (int32_t j = 0; j < n; j++) 
          { M[t] = A[t]; t++; }
      }
 
    double det = 1.0; /* Accumulates the determinant. */
 
    int32_t prow[n], pcol[n]; 
    int32_t k = 0;
    while (k < n)
      { /* Process the remaining {n-k}x{n-k} submatrix starting at {A[k][k]}. */

        /* Find the element with largest abs value in the submatrix. */
        /* Store its value in {biga}, its indices in {prow[k],pcol[k]}: */
        prow[k] = k;
        pcol[k] = k;
        double biga = M[k*(int32_t)n+k];
        for (int32_t i = k; i < n; i++)
          { for (int32_t j = k; j < n; j++)
              { if (fabs (M[i*(int32_t)n+j]) > fabs (biga))
                  { biga = M[i*(int32_t)n+j]; prow[k] = i; pcol[k] = j; }
              }
          }
 
        /* Accumulate the determinant: */
        det *= biga; /* Product of pivots */

        /* If the matrix is all zeros, we break here: */
        if (biga == 0.0) { break; }
        
        /* Swap rows {k} and {prow[k]}: */
        { int32_t i = prow[k];
          if (i > k)
            { for (int32_t j = 0; j < n; j++)
                { double hold = -M[k*(int32_t)n+j];
                  M[k*(int32_t)n+j] = M[i*(int32_t)n+j];
                  M[i*(int32_t)n+j] = hold;
                }
            }
        }
 
        /* Swap columns {k} and {pcol[k]}: */
        { int32_t j = pcol[k];
          if (j > k)
            { for (int32_t i = 0; i < n; i++)
                { double hold = -M[i*(int32_t)n+k];
                  M[i*(int32_t)n+k] = M[i*(int32_t)n+j];
                  M[i*(int32_t)n+j] = hold;
                }
            }
        }
 
        /* Negate column {k} and divide by minus the pivot {biga}. */
        for (int32_t i = 0; i < n; i++) { if (i != k) { M[i*(int32_t)n+k] /= -biga; } }
 
        /* Reduce the matrix: */
        for (int32_t i = 0; i < n; i++)
          { if (i != k)
              { double hold = M[i*(int32_t)n+k];
                for (int32_t j = 0; j < n; j++)
                  { if (j != k) { M[i*(int32_t)n+j] += hold * M[k*(int32_t)n+j]; } }
              }
          }
 
        /* Divide row {k} by pivot {biga}: */
        for (int32_t j = 0; j < n; j++) { if (j != k) { M[k*(int32_t)n+j] /= biga; } }
 
        /* Set the diagonal element to its reciprocal: */
        M[k*(int32_t)n+k] = 1.0 / biga;
         
        /* Shrink: */
        k++;
      }
 
    /* Undo the row and column interchanges, transposed: */
    while(k > 0)
      { 
        /* Back up in {k}: */
        k--;
      
        /* Swap column {k} with column {prow[k]}: */
        { int32_t i = prow[k];
          if (i > k)
            { for (int32_t j = 0; j < n; j++) 
                { double hold = M[j*(int32_t)n+k];
                  M[j*(int32_t)n+k] = -M[j*(int32_t)n+i];
                  M[j*(int32_t)n+i] = hold;
                }
            }
        }
 
        /* Swap row {k} with row {pcol[k]}: */
        { int32_t j = pcol[k];
          if (j > k)
            { for (int32_t i = 0; i < n; i++)
                { double hold = M[k*(int32_t)n+i];
                  M[k*(int32_t)n+i] = -M[j*(int32_t)n+i];
                  M[j*(int32_t)n+i] = hold;
                }
            }
        }
      }

    /* Return the determinant: */
    return det;
  }


void rmxn_scale(uint32_t m, uint32_t n, double s, double *A, double *M) 
  { int32_t k = 0;
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { M[k] = s*A[k]; k++; }
      }
  }

void rmxn_mix(uint32_t m, uint32_t n, double s, double *A, double t, double *B, double *M)
  { int32_t k = 0;
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { M[k] = s*A[k] + t*B[k]; k++; }
      }
  }

void rmxn_rel_diff(uint32_t m, uint32_t n, double *A, double *B, double *M)
  { int32_t k = 0;
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { M[k] = rel_diff(A[k], B[k]); k++; }
      }
  }

double rmxn_norm_sqr(uint32_t m, uint32_t n, double *A)
  { /* !!! Should use Kahan's summation !!! */
    double s = 0;
    int32_t k = 0;
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { double Aij = A[k]; k++; 
            s += Aij*Aij;
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
      { int32_t k = 0;
        for (int32_t i = 0; i < m; i++)
          { for (int32_t j = 0; j < n; j++)
             { A[k] /= w; k++; }
          }
      }
    else
      { int32_t k = 0;
        for (int32_t i = 0; i < m; i++)
          { for (int32_t j = 0; j < n; j++)
             { A[k] = NAN; k++; }
          }
      }
    return w;
  }

double rmxn_mod_norm_sqr(uint32_t n, double *A)
  {/* !!! Should use Kahan's summation !!! */
    double s = 0.0;
    int32_t k = 0;
    for (int32_t i = 0; i < n; i++)
      for (int32_t j = 0; j < n; j++)
        { double Aij = A[k]; k++; 
          double Dij = Aij - (i == j ? 1 : 0);
          s += Dij*Dij;
        }
    return s; 
  }

double rmxn_max_abs_elem(uint32_t m, uint32_t n, double *A)
  { double emax = 0.0;
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { double Mij = fabs(A[i*(int32_t)n + j]);
            if (Mij > emax) { emax = Mij; }
          }
      }
    return emax;
  }

void rmxn_LT_inv_map_row(uint32_t n, double *y, double *L, double *r)
  { /* Note {int} rather than {uint} because loop ends at {-1}: */
    for (int32_t j = (int32_t)n-1; j >= 0; j--)
      { double sum = y[j], corr = 0.0;
        for (int32_t i = j+1; i < n; i++)
          { double term = -r[i]*L[(int32_t)n*i + j];
            /* Kahan's summation: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
          }
        r[j] = sum/L[(int32_t)n*j + j];
      }
  }

void rmxn_LT_inv_map_col(uint32_t m, double *L, double *y, double *r)
  { for (int32_t i = 0; i < m; i++)
      { double sum = y[i], corr = 0.0;
        for (int32_t j = 0; j < i; j++)
          { double term = -L[(int32_t)m*i + j]*r[j];
            /* Kahan's summation: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
          }
        r[i] = sum/L[(int32_t)m*i + i];
      }
  }

void rmxn_LT_pos_div(uint32_t m, uint32_t n, double *A, double *L, double *M)
  { for (int32_t k = 0; k < m; k++)
      { /* Note hack for decreasing loop with uint counter. */
        for (int32_t jp1 = (int32_t)n; jp1 > 0; jp1--)
          { int32_t j = jp1-1;
            double sum = A[(int32_t)n*k + j], corr = 0.0;
            for (int32_t i = j+1; i < n; i++)
              { double term = -M[(int32_t)n*k + i]*L[(int32_t)n*i + j];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            M[(int32_t)n*k + j] = sum/L[(int32_t)n*j + j];
          }
      }
  }

void rmxn_LT_pre_div(uint32_t m, uint32_t n, double *L, double *A, double *M)
  { for (int32_t k = 0; k < n; k++)
      { for (int32_t i = 0; i < m; i++)
          { double sum = A[(int32_t)n*i + k], corr = 0.0;
            for (int32_t j = 0; j < i; j++)
              { double term = -L[(int32_t)m*i + j]*M[(int32_t)n*j + k];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            M[(int32_t)n*i + k] = sum/L[(int32_t)m*i + i];
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
      probably pronounced it himself in the French fashion, i.e.
      with an "sh" sound, as in "Chopin" and "Chostakovich"; and not
      with a "tsh" sound --- which in French would be spelled "Tch",
      as in "Tchaikovski" or "Tchekov". */
  
    for (int32_t i = 0; i < n; i++)
      { /* Compute {L[i,j] = X[i,j] - SUM{ L[i,k]*A[j,k] : k = 0..j-1 })/L[j,j]} */
        /* where X[i,j] = (j <= i ? A[i,j] : 0) */
        for (int32_t j = 0; j < i; j++)
          { double sum = 0.0, corr = 0.0;
            for (int32_t k = 0; k < j; k++)
              { double term = L[(int32_t)n*i + k]*L[(int32_t)n*j + k];
                /* Kahan's summation formula: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            double Ljj = L[(int32_t)n*j + j];
            affirm(Ljj != 0.0, "zero element in diagonal");
            L[(int32_t)n*i + j] = (A[(int32_t)n*i + j] - sum)/Ljj;
          }
        
        /* Compute {L[i,i] = sqrt(A[i,i] - SUM{ L[i,j]^2 : j = 0..i-1 })} */
        double sum = 0.0, corr = 0.0;
        for (int32_t j = 0; j < i; j++)
          { double w = L[(int32_t)n*i+j], term = w*w; 
            /* Kahan's summation formula: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
          }
        double Lii = A[(int32_t)n*i+i] - sum;
        affirm (Lii >= 0.0, "matrix is not positive definite?");
        L[(int32_t)n*i + i] = sqrt(Lii);
        for (int32_t j = i+1; j < n; j++) { L[(int32_t)n*i + j] = 0.0; }
      }
  }

void rmxn_perturb_unif(uint32_t m, uint32_t n, double pabs, double prel, double *A)
  {
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++) 
          { double *Aij = &(A[i*(int32_t)n + j]);
            double mag = pabs + prel*fabs(*Aij);
            double d = dabrandom(-mag, +mag);
            (*Aij) += d; 
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
    for (int32_t i = 0; i < m; i++)
      { if (i > 0) { fputs(osep, f); }
        for (int32_t k = 1; k <= 3; k++)
          { uint32_t nk = (k == 1 ? n1 : (k == 2 ? n2 : n3));
            double *Ak = (k == 1 ? A1 : (k == 2 ? A2 : A3));
            if (Ak != NULL)
              { double *Aki = (nk == 0 ? NULL : &(Ak[i*(int32_t)nk]));
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

double *rmxn_alloc(uint32_t m, uint32_t n)
  { double *p = talloc(m*n, double);
    return p;
  }
