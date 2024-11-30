/* See gausol_solve.h */
/* Last edited on 2024-11-30 05:30:44 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>

#include <gausol_triang.h>
#include <gausol_print.h>

#include <gausol_solve.h>

/* INTERNAL PROTOPYPES */

/* IMPLEMENTTAIONS */

void gausol_solve
  ( uint32_t m,
    uint32_t n,
    double A[],
    uint32_t p,
    double B[],
    double X[],
    bool_t pivot_rows,
    bool_t pivot_cols,
    double tiny,
    uint32_t *rank_P,
    double *det_P
  )
  { 
    bool_t debug = TRUE;
  
    if (debug) { gausol_print_system(stderr, 4, "%10.4f", "original:",  m,NULL,0, n,NULL,0,"A",A, p,"B",B, 0,NULL,NULL, ""); }
    
    uint32_t *prow = (pivot_rows ? talloc(m, uint32_t) : NULL);
    uint32_t *pcol = (pivot_cols ? talloc(n, uint32_t) : NULL);
    double det;
    uint32_t rank;
    gausol_triang_reduce(m, prow, n, pcol, A, p, B, tiny, &rank, &det);
    if (debug) { gausol_print_system(stderr, 4, "%10.4f", "triangularized:",  m,prow,rank, n,pcol,rank,"A",A, p,"B",B, 0,NULL,NULL, ""); }

    if (det_P != NULL)
      { if (debug) { fprintf(stderr, "determinant = %24.16e\n", det); }
        (*det_P) = det;
      }
    
    gausol_triang_diagonalize(m, prow, n, pcol, A, p, B, rank, tiny);
    if (debug) { gausol_print_system(stderr, 4, "%10.4f", "diagonalized:",  m,prow,rank, n,pcol,rank,"A",A, p,"B",B, 0,NULL,NULL, ""); }

    gausol_triang_normalize(m, prow, n, pcol, A, p, B, rank, tiny);
    if (debug) { gausol_print_system(stderr, 4, "%10.4f", "normalized:", m,prow,rank, n,pcol,rank,"A",A, p,"B",B, 0,NULL,NULL, ""); }

    auto void extract_solution(void);
      /* This procedure tries to set {X} to a solution of the system {A X = B},
        assuming {A} and {B} have been processed by {gausol_triang_reduce}, 
        {gausol_triang_diagonalize}, and {gausol_triang_normalize};
        so that submatrix {P10} is zeros, and {P00} into the {rank×rank} identity. Therefore the system to be
        solved is {P11 Y1 = Q1} and {Y0 = Q0 - P01 Y1}.

        Then, if {m = n = rank} then {P01}, {P11} and {Q1} are empty and
        {Y0 = Y}, so the system has a single {Y = Y0 = Q0}.

        Otherwise the system may be under- or over-determined, and 
        may have infinitely may solutions or no solution.  In that case
        the procedure will set {Y1} to all zeros, ignoring {P11} and {Q1},
        and then set {Y0} to {Q0}; so that at least the first {rank} equations
        of {P Y = Q} will be satisfied. */

    extract_solution();
    if (debug) 
      { if (rank < m) { fprintf(stderr, "there may be %d unsatisfied solutions\n", m - rank); }
        if (rank < n) { fprintf(stderr, "there are %d degrees of indeterminacy\n", n - rank); }
        gausol_print_array(stderr, 4, "%10.4f", "solution:", n,pcol,rank, p,NULL,0, "X",X, "");
      }
      
    (*rank_P) = rank;
    
    if (prow != NULL) { free(prow); }
    if (pcol != NULL) { free(pcol); }

    return;

    void extract_solution(void)
      { /* Fill {Y1} with zeros: */
        for (uint32_t j = rank; j < n; j++)
          { uint32_t pcol_j = (pcol == NULL ? j : pcol[j]); assert(pcol_j < n);
            for (uint32_t k = 0; k < p; k++) { X[pcol_j*p + k] = 0.0; }
          }

        /* Get {Y0} from {P00 Y0 = Q0 - P01 Y1}: */
        /* Since {P00} is identity and {Y1} is zero, it is just {Y0 = Q0}: */
        for (uint32_t t = 0; t < rank; t++)
          { uint32_t prow_t = (prow == NULL ? t : prow[t]); assert(prow_t < m);
            uint32_t pcol_t = (pcol == NULL ? t : pcol[t]); assert(pcol_t < n);
            for (uint32_t k = 0; k < p; k++) { X[pcol_t*p + k] = B[prow_t*p + k]; }
          }
      }
    
  }
