/* See rfmxn.h. */
/* Last edited on 2021-08-17 05:54:53 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <rn.h>
#include <rmxn.h>
#include <bool.h>
#include <jsmath.h>
#include <affirm.h>
#include <gauss_elim.h>

#include <rfn.h>
#include <rfmxn.h>

void rfmxn_zero(int32_t m, int32_t n, float *M)
  { int32_t i, j, t;
    t = 0;
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        { M[t] = 0.0; t++; }
  }

void rfmxn_copy(int32_t m, int32_t n, float *A, float *M)
  { int32_t mn = m*n;
    int32_t ij;
    for (ij = 0; ij < mn; ij++) { M[ij] = A[ij]; }
  }

void rfmxn_get_row(int32_t m, int32_t n, float *A, int32_t i, float *r)
  { float *Ai = &(A[i*n]);
    int32_t j;
    for (j = 0; j < n; j++) { r[j] = Ai[j]; }
  }
  
void rfmxn_set_row(int32_t m, int32_t n, float *A, int32_t i, float *r)
  { float *Ai = &(A[i*n]);
    int32_t j;
    for (j = 0; j < n; j++) { Ai[j] = r[j]; }
  }

void rfmxn_get_col(int32_t m, int32_t n, float *A, int32_t j, float *r)
  { float *Aij = &(A[j]);
    int32_t i;
    for (i = 0; i < m; i++) { r[i] = (*Aij); Aij += n; }
  }

void rfmxn_set_col(int32_t m, int32_t n, float *A, int32_t j, float *r)
  { float *Aij = &(A[j]);
    int32_t i;
    for (i = 0; i < m; i++) { (*Aij) = r[i]; Aij += n; }
  }

void rfmxn_ident(int32_t m, int32_t n, float *M)
  { int32_t i, j, t;
    t = 0;
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        { M[t] = (i == j ? 1.0 : 0.0); t++; }
  }

void rfmxn_map_row (int32_t m, int32_t n, float *x, float *A, float *r)
  { int32_t i, j;
    for (j = 0; j < n; j++)
      { double sum = 0.0, corr = 0.0; 
        int32_t t = j;
        for (i = 0; i < m; i++) 
          { double term = (double)(x[i]) * A[t];
            /* Kahan's summation: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
            t += n;
          }
        r[j] = (float)sum;
      }
  }

void rfmxn_map_col (int32_t m, int32_t n, float *A, float *x, float *r)
  { int32_t i, j, t = 0;
    for (i = 0; i < m; i++)
      { double sum = 0.0, corr = 0.0; 
        for (j = 0; j < n; j++)
          { double term = (double)(A[t]) * x[j]; 
            /* Kahan's summation: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
            t++;
          }
        r[i] = (float)sum;
      }
  }

void rfmxn_mul (int32_t m, int32_t p, int32_t n, float *A, float *B, float *M)
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
                t += n;
              }
            M[v] = (float)sum; v++;
          }
        r += p;
      }
  }

void rfmxn_mul_tr (int32_t m, int32_t n, int32_t p, float *A, float *B, float *M)
  { int32_t i, j, k;
    int32_t v = 0;
    int32_t r = 0;
    for (i = 0; i < m; i++)
      { int32_t s = 0;
        for (j = 0; j < n; j++)
          { double sum = 0.0, corr = 0.0;
            for (k = 0; k < p; k++) 
              { double term = (double)(A[r+k])*B[s+k];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            M[v] = (float)sum; v++;
            s += p;
          }
        r += p;
      }
  }

void rfmxn_tr_mul (int32_t p, int32_t m, int32_t n, float *A, float *B, float *M)
  { int32_t i, j, k;
    int32_t v = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { double sum = 0.0, corr = 0.0;
            int32_t r = i, s = j;
            for (k = 0; k < p; k++) 
              { double term = (double)(A[r])*B[s];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
                r += m; s += n;
              }
            M[v] = (float)sum; v++;
          }
      }
  }

double rfmxn_det (int32_t n, float *A)
  { int32_t n2 = n*n;
    double *C = (double *)notnull(malloc(n2*sizeof(double)), "no mem for C");
    for (int32_t t = 0; t < n2; t++) { C[t] = (double)A[t]; }
    gsel_triangularize(n, n, C, FALSE, 0.0);
    double det = 1.0;
    for (int32_t t = 0; t < n2; t += n+1) { det *= C[t]; }
    free(C);
    return det;
  }

double rfmxn_inv (int32_t n, float *A, float *M)
  { int32_t n2 = n*n;
    double *C = (double *)notnull(malloc(n2*sizeof(double)), "no mem for C");
    for (int32_t t = 0; t < n2; t++) { C[t] = (double)A[t]; }
    double det = rmxn_inv(n, C, C);
    for (int32_t t = 0; t < n2; t++) { M[t] = (float)C[t]; }
    free(C);
    return det;
  }
  
double rfmxn_inv_full (int32_t n, float *A, float *M)
  { int32_t n2 = n*n;
    double *C = (double *)notnull(malloc(n2*sizeof(double)), "no mem for C");
    for (int32_t t = 0; t < n2; t++) { C[t] = (double)A[t]; }
    double det = rmxn_inv_full(n, C, C);
    for (int32_t t = 0; t < n2; t++) { M[t] = (float)C[t]; }
    free(C);
    return det;
  }


void rfmxn_scale(int32_t m, int32_t n, float s, float *A, float *M) 
  { int32_t i, j, k = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { M[k] = (float)(s*A[k]); k++; }
      }
  }

void rfmxn_mix(int32_t m, int32_t n, double s, float *A, double t, float *B, float *M)
  { int32_t i, j, k = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { M[k] = (float)(s*A[k] + t*B[k]); k++; }
      }
  }

void rfmxn_rel_diff(int32_t m, int32_t n, float *A, float *B, float *M)
  { int32_t i, j, k = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { M[k] = (float)rel_diff((double)A[k], (double)B[k]); k++; }
      }
  }

double rfmxn_norm_sqr(int32_t m, int32_t n, float *A)
  { double s = 0;
    int32_t i, j, k = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { double Aij = (double)A[k]; k++; 
            s += Aij*Aij;
          }
      }
    return s;
  }

double rfmxn_norm(int32_t m, int32_t n, float *A)
  { double s = 0;
    int32_t i, j, k = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { double Aij = (double)A[k]; k++; 
            s += Aij*Aij;
          }
      }
    return sqrt(s);
  }

double rfmxn_mod_norm_sqr(int32_t n, float *A)
  {
    double s = 0.0;
    int32_t i, j;
    int32_t k = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
        { double Aij = (double)A[k]; k++; 
          double Dij = Aij - (i == j ? 1 : 0);
          s += Dij*Dij;
        }
    return s; 
  }

double rfmxn_max_abs_elem(int32_t m, int32_t n, float *A)
  { int32_t i, j;
    float emax = 0.0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { float Mij = fabsf(A[i*n + j]);
            if (Mij > emax) { emax = Mij; }
          }
      }
    return emax;
  }

void rfmxn_LT_inv_map_row(int32_t n, float *y, float *L, float *r)
  { int32_t i, j;
    for (j = n-1; j >= 0; j--)
      { double sum = (double)y[j], corr = 0.0;
        for (i = j+1; i < n; i++)
          { double term = -(double)(r[i])*L[n*i + j];
            /* Kahan's summation: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
          }
        r[j] = (float)(sum/L[n*j + j]);
      }
  }

void rfmxn_LT_inv_map_col(int32_t m, float *L, float *y, float *r)
  { int32_t i, j;
    for (i = 0; i < m; i++)
      { double sum = (double)y[i], corr = 0.0;
        for (j = 0; j < i; j++)
          { double term = -(double)(L[m*i + j])*r[j];
            /* Kahan's summation: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
          }
        r[i] = (float)(sum/L[m*i + i]);
      }
  }

void rfmxn_LT_pos_div(int32_t m, int32_t n, float *A, float *L, float *M)
  { int32_t i, j, k;
    for (k = 0; k < m; k++)
      { for (j = n-1; j >= 0; j--)
          { double sum = (double)A[n*k + j], corr = 0.0;
            for (i = j+1; i < n; i++)
              { double term = -(double)(M[n*k + i])*L[n*i + j];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            M[n*k + j] = (float)(sum/L[n*j + j]);
          }
      }
  }

void rfmxn_LT_pre_div(int32_t m, int32_t n, float *L, float *A, float *M)
  { int32_t i, j, k;
    for (k = 0; k < n; k++)
      { for (i = 0; i < m; i++)
          { double sum = (double)A[n*i + k], corr = 0.0;
            for (j = 0; j < i; j++)
              { double term = -(double)(L[m*i + j])*M[n*j + k];
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            M[n*i + k] = (float)(sum/L[m*i + i]);
          }
      }
  }

void rfmxn_cholesky(int32_t n, float *A, float *L)
  { 
    int32_t n2 = n*n;
    double *C = (double *)notnull(malloc(n2*sizeof(double)), "no mem for C");
    for (int32_t t = 0; t < n2; t++) { C[t] = (double)A[t]; }
    rmxn_cholesky(n, C, C);
    for (int32_t t = 0; t < n2; t++) { L[t] = (float)C[t]; }
    free(C);
  }

void rfmxn_print (FILE *f, int32_t m, int32_t n, float *A)
  { rfmxn_gen_print(f, m, n, A, NULL, NULL, NULL, NULL, NULL, NULL, NULL); }

void rfmxn_gen_print
  ( FILE *f, int32_t m, int32_t n, float *A,
    char *fmt, 
    char *olp, char *osep, char *orp,
    char *ilp, char *isep, char *irp
  )
  {
    int32_t i,j, t;
    if (olp == NULL) { olp = "(\n"; }
    if (osep == NULL) { osep = "\n"; }
    if (orp == NULL) { orp = "\n)"; }
    if (ilp == NULL) { ilp = "  ("; }
    if (isep == NULL) { isep = " "; }
    if (irp == NULL) { irp = ")"; }
    if (fmt == NULL) { fmt = "%13.6e"; }
    fputs(olp, f);
    t = 0;
    for (i = 0; i < m; i++)
      {
        if (i > 0) { fputs(osep, f); }
        fputs(ilp, f);
        for (j = 0; j < n; j++) 
          { if (j > 0) { fputs(isep, f); }
            fprintf(f, fmt, A[t]); t++;
          }
        fputs(irp, f);
      }
    fputs(orp, f);
    fflush(f);
  }  


float *rfmxn_alloc(int32_t m, int32_t n)
  { void *p = malloc(m*n*sizeof(float));
    affirm(p != NULL, "no memory for rfmxn_t");
    return (float *)p;
  }

