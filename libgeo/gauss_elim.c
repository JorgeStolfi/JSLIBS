/* See gauss_elim.h */
/* Last edited on 2021-06-09 20:14:57 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <math.h>
#include <limits.h>

#include <assert.h>
#include <affirm.h>
#include <bool.h>
#include <rmxn.h>

#include <gauss_elim.h>

/* IMPLEMENTATIONS */

int32_t gsel_solve(int32_t m, int32_t n, double A[], int32_t p, double B[], double X[], double tiny)
  { bool_t debug = FALSE;
  
    /* Work array: */
    int32_t np = n+p; /* Total columns in {A} and {B}. */
    double AB[m*np]; /* Matrices {A} and {B} side by side. */
    
    /* Copy system matrices into work array: */
    int32_t i, j, k;
    for (i = 0; i < m; i ++)
      { for (j = 0; j < n; j++) { AB[i*np + j] = A[i*n + j]; }
        for (k = 0; k < p; k++) { AB[i*np + n + k] = B[i*p + k]; }
      }
    
    /* Solve system: */
    if (debug) { gsel_print_array(stderr, "%9.5f", "original:",  m, np, AB, ""); }
    
    gsel_triangularize(m, np, AB, TRUE, tiny);
    if (debug) { gsel_print_array(stderr, "%9.5f", "triangularized:",  m, np, AB, ""); }

    gsel_diagonalize(m, np, AB, TRUE);
    if (debug) { gsel_print_array(stderr, "%9.5f", "diagonalized:",  m, np, AB, ""); }

    gsel_normalize(m, np, AB, TRUE);
    if (debug) { gsel_print_array(stderr, "%9.5f", "normalized:", m, np, AB, ""); }

    int32_t r = gsel_extract_solution(m, np, AB, p, X, TRUE);
    if (debug) 
      { if (r < m) { fprintf(stderr, "there may be %d unsatisfied solutions\n", m - r); }
        if (r < np-p) { fprintf(stderr, "there are %d degrees of indeterminacy\n", np - p - r); }
        gsel_print_array(stderr, "%9.5f", "solution:", n, p, X, "");
      }
    
    return r;
  }
  
double gsel_determinant(int32_t m, int32_t n, double A[], int32_t q)
  { if ((q > m) || (q > n)) { return 0.0; }
    /* Make a work copy of the first {q} rows and columns of {A}: */
    double M[q*q];
    int32_t i, j;
    for (i = 0; i < q; i ++)
      { for (j = 0; j < q; j++) { M[i*q + j] = A[i*n + j]; } }
    
    /* Triangularize and get determinant: */
    gsel_triangularize(q, q, M, FALSE, 0.0);
    return gsel_triangular_det(q, q, M, q);
  }

void gsel_print_system
  ( FILE *wr, 
    char *fmt, 
    char *head, 
    int32_t m, 
    int32_t n, 
    double A[], 
    int32_t p, 
    double B[], 
    char *foot
  )
  { if (head != NULL) { fprintf(wr, "%s\n", head); }
    int32_t i, j, k;
    for (i = 0; i < m; i++)
      { fprintf(wr, "  ");
        for (j = 0; j < n; j++)
          { fprintf(wr, " "); fprintf(wr, fmt, A[i*n + j]); }
        fprintf(wr, " | "); 
        for (k = 0; k < p; k++)
          { fprintf(wr, " "); fprintf(wr, fmt, B[i*p + k]); }
        fprintf(wr, "\n");
      }
    if (foot != NULL) { fprintf(wr, "%s\n", foot); }
  }

void gsel_print_array(FILE *wr, char *fmt, char *head, int32_t m, int32_t n, double M[], char *foot)
  { if (head != NULL) { fprintf(wr, "%s\n", head); }
    int32_t i, j;
    for (i = 0; i < m; i++)
      { fprintf(wr, "  ");
        for (j = 0; j < n; j++)
          { fprintf(wr, " "); fprintf(wr, fmt, M[i*n + j]); }
        fprintf(wr, "\n");
      }
    if (foot != NULL) { fprintf(wr, "%s\n", foot); }
  }

void gsel_triangularize(int32_t m, int32_t n, double M[], int32_t total, double tiny)
  { int32_t i = 0;
    int32_t j = 0;
    double *Mij = &(M[0]);
    while ((i < m-1) && (j < n))
      { /* Elements {M[r][s]} with {r >= i} and {s < j} are all zero. */
        /* Clear elements {M[k][j]} with {k > i} */
        int32_t k; double *Mkj;
        for (k = i+1, Mkj = Mij + n; k < m; k++, Mkj += n)
          { /* Swap rows {i} and {k} if needed, negating one of them: */
            if (fabs(*Mkj) > fabs(*Mij))
              { int32_t r; double *Mkr; double *Mir;
                for (r = j, Mkr = Mkj, Mir = Mij; r < n; r++, Mkr++, Mir++)
                  { double t = (*Mkr); (*Mkr) = - (*Mir); (*Mir) = t; }
              }
            /* Subtract from row {k} a multiple of row {i} that cancels {M[k,j]}: */
            if ((*Mkj) != 0.0)
              { int32_t r; double *Mkr; double *Mir;
                double s = (*Mkj)/(*Mij);
                (*Mkj) = 0.0;
                for (r = j+1, Mkr = Mkj+1, Mir = Mij+1; r < n; r++, Mkr++, Mir++)
                  { double old = (*Mkr);
                    (*Mkr) = old - s*(*Mir);
                    /* Feeble attempt to clean out entries that are just roundoff error: */
                    if (1.0e12*fabs(*Mkr) < fabs(old)) { (*Mkr) = 0.0; }
                  }
              }
          }
        if ((! total) || ((*Mij) != 0.0)) { i++; Mij += n; }
        j++; Mij++;
      }
  }

void gsel_diagonalize(int32_t m, int32_t n, double M[], int32_t total)
  { int32_t i = 0;
    int32_t j = 0;
    double *Mij = &(M[0]);
    while ((i < m) && (j < n))
      { if ((*Mij) != 0.0)
          { /* Clear elements {M[k][j]} with {k < i} */
            int32_t k; double *Mkj;
            for (k = i-1, Mkj = Mij-n; k >= 0; k--, Mkj -= n)
              { /* Sub from row {k} a multiple of row {i} that cancels {M[k,j]}: */
                if ((*Mkj) != 0.0)
                  { int32_t r; double *Mkr; double *Mir;
                    double s = (*Mkj)/(*Mij);
                    (*Mkj) = 0.0;
                    for (r = j+1, Mkr = Mkj+1, Mir = Mij+1; r < n; r++, Mkr++, Mir++)
                      { (*Mkr) -= s*(*Mir); }
                  }
              }
            i++; Mij += n;
          }
        else if (! total)
          { i++; Mij += n; }
        j++; Mij++;
      }
  }

void gsel_normalize(int32_t m, int32_t n, double M[], int32_t total)
  { int32_t i = 0;
    int32_t j = 0;
    double *Mij = &(M[0]);
    while ((i < m) && (j < n))
      { if ((*Mij) != 0.0)
          { /* Scale row {i} by {1/M[i,j]}: */
            int32_t r; double *Mir;
            double s = (*Mij);
            (*Mij) = 1.0;
            for (r = j+1, Mir = Mij+1; r < n; r++, Mir++) { (*Mir) /= s; }
            i++; Mij += n;
          }
        else if (! total)
          { i++; Mij += n; }
        j++; Mij++;
      }
  }

int32_t gsel_extract_solution(int32_t m, int32_t n, double M[], int32_t p, double X[], int32_t total)
  { affirm(n >= p, "bad array dimensions");
    int32_t q = n-p;
    int32_t i, j, k;
    /* Scan unknowns and set them: */
    i = 0; j = 0;
    double *Aij = &(M[0]);
    int32_t r = 0; /* Number of equations actually used. */
    while (j < q)
      { /* Set unknowns {X[j,0..p-1]}: */
        double *Xjk = &(X[j*p]);
        double piv = (i < m ? (*Aij) : 0.0); /* Pivot value. */
        if (piv != 0.0)
          { /* Equation {i} defines {X[j,0..p-1]}: */
            double *Bik = &(M[i*n + q]);
            for (k = 0; k < p; k++, Xjk++, Bik++) { (*Xjk) = (*Bik)/piv; }
            r++;
            if (i < m) { i++; Aij += n; }
          }
        else
          { /* Equation {i} skips variables {X[j,0..p-1]}, so set them to zero: */
            for (k = 0; k < p; k++, Xjk++) { (*Xjk) = 0.0; }
            if ((i < m) && (! total)) { i++; Aij += n; }
          }
        j++; Aij++;
      }
    return r;
  }

double gsel_triangular_det(int32_t m, int32_t n, double M[], int32_t q)
  { if ((q > m) || (q > n)) { return 0.0; }
    double det = 1.0;
    int32_t i;
    for (i = 0; i < q; i++) { det *= M[i*n + i];}
    return det;
  }

void gsel_residual(int32_t m, int32_t n, double A[], int32_t p, double B[], double X[], double R[])
  { /* Compute {A X - B} with care: */
    int32_t i, j, k;
    int32_t in0 = 0; /* {== i*n}. */
    int32_t ipj = 0; /* {== i*p+j}. */
    for (i = 0; i < m; i++)
      { for (j = 0; j < p; j++)
          { double sum = 0.0, corr = 0.0;
            int32_t kpj = j; /* {== k*p+j}. */
            for (k = 0; k <= n; k++)
              { double term = (k == n ? -B[ipj] : A[in0+k]*X[kpj]);
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
                kpj += p;
              }
            R[ipj] = sum; ipj++;
          }
        in0 += n;
      }
  }
