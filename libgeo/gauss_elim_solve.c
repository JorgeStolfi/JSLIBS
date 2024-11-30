/* See gauss_elim_solve.h */
/* Last edited on 2024-11-25 04:14:30 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>

#include <gauss_elim.h>
#include <gauss_elim_triangularize.h>
#include <gauss_elim_diagonalize.h>
#include <gauss_elim_normalize.h>

#include <gauss_elim_solve.h>

uint32_t gauss_elim_solve
  ( uint32_t m,
    uint32_t n,
    double A[],
    uint32_t p,
    double B[],
    double X[],
    double tiny
  )
  { /* Work array: */
    uint32_t np = n+p; /* Total columns in {A} and {B}. */
    double AB[m*np]; /* Matrices {A} and {B} side by side. */
    
    /* Copy system matrices into work array: */
    
    for (int32_t i = 0; i < m; i ++)
      { for (int32_t j = 0;  j < n; j++) { AB[i*(int32_t)np + j] = A[i*(int32_t)n + j]; }
        for (int32_t k = 0;  k < p; k++) { AB[i*(int32_t)np + (int32_t)n + k] = B[i*(int32_t)p + k]; }
      }
    uint32_t rank;
    gauss_elim_solve_packed(m, n, p, AB, X, tiny, &rank, NULL);
    return rank;
  }
    
void gauss_elim_solve_packed
  ( uint32_t m,
    uint32_t n,
    uint32_t p,
    double AB[],
    double X[],
    double tiny,
    uint32_t *rank_P,
    double *det_P
  )
  {
    bool_t debug = FALSE;
  
    uint32_t np = n+p; /* Num of columns in {AB}. */
    
    /* Solve system: */
    if (debug) { gauss_elim_print_array(stderr, 4, "%9.5f", "original:",  m, np,"AB",AB, ""); }
    
    gauss_elim_triangularize(m, np, AB, TRUE, tiny);
    if (debug) { gauss_elim_print_array(stderr, 4, "%9.5f", "triangularized:",  m, np,"AB",AB, ""); }
    
    gauss_elim_diagonalize(m, np, AB);
    if (debug) { gauss_elim_print_array(stderr, 4, "%9.5f", "diagonalized:",  m, np,"AB",AB, ""); }
    
    double det;
    if (m != n)
      { det = 0.0; }
    else
      { det = 1.0;
        for (uint32_t i = 0; (i < n) && (det != 0.0); i++)
          { det *= AB[i*np + i]; }
      }
    if (debug) { fprintf(stderr, "determinant = %24.16e\n", det); }

    gauss_elim_normalize(m, np, AB);
    if (debug) { gauss_elim_print_array(stderr, 4, "%9.5f", "normalized:", m, np,"AB",AB, ""); }

    uint32_t rank = gauss_elim_extract_solution(m, np, AB, p, X);
    if (debug) 
      { if (rank < m) { fprintf(stderr, "there may be %d unsatisfied solutions\n", m - rank); }
        if (rank < np-p) { fprintf(stderr, "there are %d degrees of indeterminacy\n", np - p - rank); }
        gauss_elim_print_array(stderr, 4, "%9.5f", "solution:", n, p,"X",X, "");
      }
    
    if (rank_P != NULL) (*rank_P) = rank;
    if (det_P != NULL) (*det_P) = det;
  }
 
uint32_t gauss_elim_extract_solution
  ( uint32_t m,
    uint32_t n,
    double M[],
    uint32_t p,
    double X[]
  )
  { demand(n >= p, "bad array dimensions");
    uint32_t q = n-p;
    
    /* Scan unknowns and set them: */
    int32_t i = 0;
    int32_t j = 0; /* Elements of {M[i,k]} with {k < j} should be all zero. */
    double *Aij = &(M[0]);
    uint32_t neq = 0; /* Number of equations actually used. */
    while (j < q)
      { /* Set unknowns {X[j,0..p-1]}: */
        double *Xjk = &(X[j*(int32_t)p]);
        double piv = (i < m ? (*Aij) : 0.0); /* Pivot value. */
        demand(isfinite(piv), "invalid element in matrix");
        if (piv != 0.0)
          { /* Equation {i} defines {X[j,0..p-1]}: */
            double *Bik = &(M[i*(int32_t)n + (int32_t)q]);
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
