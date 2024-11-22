/* See gauss_elim.h */
/* Last edited on 2024-11-22 02:16:21 by stolfi */

#include <stdint.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include <assert.h>
#include <affirm.h>
#include <bool.h>
#include <rmxn.h>
#include <jsmath.h>

#include <gauss_elim.h>

/* IMPLEMENTATIONS */

uint32_t gsel_solve(uint32_t m, uint32_t n, double A[], uint32_t p, double B[], double X[], double tiny)
  { bool_t debug = FALSE;
  
    /* Work array: */
    uint32_t np = n+p; /* Total columns in {A} and {B}. */
    double AB[m*np]; /* Matrices {A} and {B} side by side. */
    
    /* Copy system matrices into work array: */
    
    for (int32_t i = 0; i < m; i ++)
      { for (int32_t j = 0; j < n; j++) { AB[i*np + j] = A[i*n + j]; }
        for (int32_t k = 0; k < p; k++) { AB[i*np + n + k] = B[i*p + k]; }
      }
    
    /* Solve system: */
    if (debug) { gsel_print_array(stderr, 4, "%9.5f", "original:",  m, np,"AB",AB, ""); }
    
    gsel_triangularize(m, np, AB, TRUE, tiny);
    if (debug) { gsel_print_array(stderr, 4, "%9.5f", "triangularized:",  m, np,"AB",AB, ""); }
    
    gsel_diagonalize(m, np, AB);
    if (debug) { gsel_print_array(stderr, 4, "%9.5f", "diagonalized:",  m, np,"AB",AB, ""); }

    gsel_normalize(m, np, AB);
    if (debug) { gsel_print_array(stderr, 4, "%9.5f", "normalized:", m, np,"AB",AB, ""); }

    uint32_t rank_ext = gsel_extract_solution(m, np, AB, p, X);
    if (debug) 
      { if (rank_ext < m) { fprintf(stderr, "there may be %d unsatisfied solutions\n", m - rank_ext); }
        if (rank_ext < np-p) { fprintf(stderr, "there are %d degrees of indeterminacy\n", np - p - rank_ext); }
        gsel_print_array(stderr, 4, "%9.5f", "solution:", n, p,"X",X, "");
      }
    
    return rank_ext;
  }
  
double gsel_determinant(uint32_t m, uint32_t n, double A[], uint32_t q)
  { if ((q > m) || (q > n)) { return 0.0; }
    /* Make a work copy of the first {q} rows and columns of {A}: */
    double M[q*q];
    
    for (int32_t i = 0; i < q; i ++)
      { for (int32_t j = 0; j < q; j++) { M[i*q + j] = A[i*n + j]; } }
    
    /* Triangularize and get determinant: */
    gsel_triangularize(q, q, M, FALSE, 0.0);
    return gsel_triangular_det(q, q, M, q);
  }

void gsel_triangularize(uint32_t m, uint32_t n, double M[], bool_t total, double tiny)
  { int32_t i = 0;
    int32_t j = 0;
    double *Mij = &(M[0]);
    while ((i < m) && (j < n))
      { double *Mkj;
        /* Elements {M[r][s]} with {r >= i} and {s < j} are all zero. */
        assert(isfinite(*Mij));
        if (i < m-1)
          { /* Pivoting: */
            /* Find row {kmax} in {i..m-1} that maximizes {|M[kmax,j]|} */
            int32_t kmax = i;
            double Mmax = (*Mij);
            Mkj = Mij + n;
            for (int32_t k = i+1; k < m; k++)
              { assert(isfinite(*Mkj));
                if (fabs(*Mkj) > fabs(Mmax)) { kmax = k; Mmax = (*Mkj); } 
                Mkj += n;
              }
            if (kmax != i)
              { /* Swap rows {i} and {kmax} negating one of them: */
                Mkj = &(M[kmax*n + j]);
                double *Mkr = Mkj; 
                double *Mir = Mij;
                for (int32_t r = j; r < n; r++)
                  { double t = (*Mkr); (*Mkr) = -(*Mir); (*Mir) = t;
                    Mkr++; Mir++;
                  }
              }
          }
        if ((*Mij) != 0.0)
          { /* Clear elements {M[k][j]} with {k > i} */
            Mkj = Mij + n;
            for (int32_t k = i+1; k < m; k++)
              { if ((*Mkj) != 0.0)
                  { /* Subtract from row {k} a multiple of row {i} that cancels {M[k,j]}: */
                    double s = (*Mkj)/(*Mij);
                    assert(isfinite(s));
                    (*Mkj) = 0.0;
                    double *Mkr = Mkj+1;
                    double *Mir = Mij+1;
                    for (int32_t r = j+1; r < n; r++)
                      { double old = (*Mkr);
                        (*Mkr) = old - s*(*Mir);
                        /* Feeble attempt to clean out entries that are just roundoff error: */
                        if (fabs(*Mkr) < tiny*fabs(old)) { (*Mkr) = 0.0; }
                        Mkr++; Mir++;
                      }
                  }
                Mkj += n;
              }
            /* Advance to next row: */
            i++; Mij += n;
          }
        else
          { if (! total) { /* Advance to next row anyway: */ i++; Mij += n; } }
        j++; Mij++;
      }
  }

void gsel_diagonalize(uint32_t m, uint32_t n, double M[])
  { int32_t i = 0;
    int32_t j = 0;
    double *Mij = &(M[0]);
    while ((i < m) && (j < n))
      { if ((*Mij) != 0.0)
          { /* Clear elements {M[k][j]} with {k < i} */
            double *Mkj = Mij-n;
            for (int32_t k = (int32_t)i-1; k >= 0; k--)
              { /* Sub from row {k} a multiple of row {i} that cancels {M[k,j]}: */
                if ((*Mkj) != 0.0)
                  { double s = (*Mkj)/(*Mij);
                    (*Mkj) = 0.0;
                    double *Mkr = Mkj+1; 
                    double *Mir = Mij+1;
                    for (int32_t r = j+1; r < n; r++)
                      { (*Mkr) -= s*(*Mir); Mkr++; Mir++; }
                  }
                Mkj -= n;
              }
            i++; Mij += n;
          }
        j++; Mij++;
      }
  }

void gsel_normalize(uint32_t m, uint32_t n, double M[])
  { int32_t i = 0;
    int32_t j = 0;
    double *Mij = &(M[0]);
    while ((i < m) && (j < n))
      { if ((*Mij) != 0.0)
          { /* Scale row {i} by {1/M[i,j]}: */
            double s = (*Mij);
            (*Mij) = 1.0;
            double *Mir = Mij+1;
            for (int32_t r = j+1; r < n; r++) { (*Mir) /= s;  Mir++; }
            i++; Mij += n;
          }
        j++; Mij++;
      }
  }

uint32_t gsel_extract_solution(uint32_t m, uint32_t n, double M[], uint32_t p, double X[])
  { affirm(n >= p, "bad array dimensions");
    uint32_t q = n-p;
    
    /* Scan unknowns and set them: */
    int32_t i = 0;
    int32_t j = 0; /* Elements of {M[i,k]} with {k < j} should be all zero. */
    double *Aij = &(M[0]);
    uint32_t neq = 0; /* Number of equations actually used. */
    while (j < q)
      { /* Set unknowns {X[j,0..p-1]}: */
        double *Xjk = &(X[j*p]);
        double piv = (i < m ? (*Aij) : 0.0); /* Pivot value. */
        demand(isfinite(piv), "invalid element in matrix");
        if (piv != 0.0)
          { /* Equation {i} defines {X[j,0..p-1]}: */
            double *Bik = &(M[i*n + q]);
            for (int32_t k = 0; k < p; k++, Xjk++, Bik++) { (*Xjk) = (*Bik)/piv; }
            neq++;
            if (i < m) { i++; Aij += n; }
          }
        else
          { /* Equation {i} skips variables {X[j,0..p-1]}, so set them to zero: */
            for (int32_t k = 0; k < p; k++, Xjk++) { (*Xjk) = 0.0; }
          }
        j++; Aij++;
      }
    return neq;
  }

double gsel_triangular_det(uint32_t m, uint32_t n, double M[], uint32_t q)
  { if ((q > m) || (q > n)) { return 0.0; }
    double det = 1.0;
    
    for (int32_t i = 0; i < q; i++) { det *= M[i*n + i];}
    return det;
  }

void gsel_residual(uint32_t m, uint32_t n, double A[], uint32_t p, double B[], double X[], double R[])
  { /* Compute {A X - B} with care: */
    
    int32_t in0 = 0; /* {== i*n}. */
    int32_t ipj = 0; /* {== i*p+j}. */
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < p; j++)
          { double sum = 0.0, corr = 0.0;
            int32_t kpj = j; /* {== k*p+j}. */
            for (int32_t k = 0; k <= n; k++)
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

void gsel_print_array
  ( FILE *wr,
    uint32_t indent,
    char *fmt,
    char *head,
    uint32_t m,
    uint32_t n,
    char *Mname,
    double M[],
    char *foot
  )
  { 
    gsel_print_system(wr, indent, fmt, head, m, n,Mname,M, 0,NULL,NULL, 0,NULL,NULL, foot);
  }

void gsel_print_system
  ( FILE *wr, 
    uint32_t indent,
    char *fmt, 
    char *head, 
    uint32_t m, 
    uint32_t n, 
    char *Aname,
    double A[], 
    uint32_t p, 
    char *Bname,
    double B[], 
    uint32_t q,
    char *Cname,
    double C[], 
    char *foot
  )
  { if (Aname == NULL) { Aname = ""; }
    char *Aeq = (strlen(Aname) > 0 ? " = " : "");
    
    if (Bname == NULL) { Bname = ""; }
    char *Beq = (strlen(Bname) > 0 ? " = " : "");
    
    if (Cname == NULL) { Cname = ""; }
    char *Ceq = (strlen(Cname) > 0 ? " = " : "");
    
    auto void print_row(uint32_t i, char *name, char *eq, uint32_t r, double M[]);
      /* Prints a row of {M} on the current line, without newline. */
    
    if (head != NULL) { fprintf(wr, "%*s%s\n", indent, "", head); }
    for (int32_t i = 0; i < m; i++)
      { if (indent > 0) { fprintf(wr, "%*s", indent, ""); }
        if (A != NULL) { print_row(i, Aname, Aeq, n, A); }
        if (B != NULL) { print_row(i, Bname, Beq, p, B); }
        if (C != NULL) { print_row(i, Cname, Ceq, q, C); }
        fprintf(wr, "\n");
      }
    if (foot != NULL) { fprintf(wr, "%*s%s\n", indent, "", foot); }
    return;
    
    void print_row(uint32_t i, char *name, char *eq, uint32_t r, double M[])
      { int32_t skip = (int32_t)(strlen(name) + strlen(eq));
        if (i == (m-1)/2)
          { fprintf(wr, "  %s%s[ ", name, eq); }
        else
          { fprintf(wr, "  %*s[ ", skip, ""); }
        for (int32_t k = 0; k < r; k++)
          { fprintf(wr, " "); fprintf(wr, fmt, M[i*r + k]); }
        fprintf(wr, " ]");
      }
  }
