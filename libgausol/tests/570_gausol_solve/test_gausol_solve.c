/* test_gausol_solve.c --- test program for gausol_solve.h  */
/* Last edited on 2024-11-30 05:26:21 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include <jsmath.h>
#include <jsrandom.h>
#include <bool.h>
#include <affirm.h>

#include <gausol_print.h>
#include <gausol_test_tools.h>

#include <gausol_solve.h>

/* GENERAL PARAMETERS */

#define MAX_RUNS 100
  /* Max number of trials per test. */

#define MAX_ROWS 10
  /* Max number of rows in main matrices of generated linear systems. */

#define MAX_COLS 10
  /* Max number of columns in main matrices of generated linear systems. */

#define MAX_PRBS 3
  /* Max number of columns in right-hand sides of generated linear systems. */

/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);

void test_gausol_solve(uint32_t trial, bool_t verbose);
  /* Tests {gausol_solve}. */

void test_gausol_solve_packed(uint32_t trial, bool_t verbose);
  /* Tests {gausol_solve_packed}. */

void check_solve
  ( uint32_t m, uint32_t n, double A[],
    uint32_t p, double B[],
    double X[], uint32_t rank,
    double Bmax[]
  );
  /* Checks the outcome of {gausol_solve}. The parameter {X} should be
    the solution computed by {gausol_solve}, and {rank} should be its
    return value. */

void check_residual
  ( uint32_t m, uint32_t n, double A[],
    uint32_t p, double B[],
    double X[],
    double R[],
    double Bmax[]
  );
  /* Checks the outcome of {gausol_residual}. The parameter {R} should be
    the residual computed by {gausol_residual} for the given {A}, {B},
    and {X}. */

void check_solution_with_reference
  ( uint32_t m, uint32_t n, uint32_t p,
    double X[], double X_ref[],
    double Bmax[]
  );
  /* Checks a putative solution {X} of a system {A X = B} against the
    `true' solution {X_ref}; where {A} is {m × n}, {X} and {X_ref} are
    {n × p}, and {B} is {m × p}. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    for (uint32_t i = 0;  i < MAX_RUNS; i++) 
      { test_gausol_solve(i, i < 5); }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_gausol_solve (uint32_t trial, bool_t verbose)
  {
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);

    uint32_t m, n, p;
    double *A, *B, *X_ref;
    uint32_t m_max = MAX_ROWS;
    uint32_t n_max = MAX_COLS;
    uint32_t p_max = MAX_PRBS;
    
    double tiny = ((trial ^ 1) == 0 ? 1.0e-180 : 1.0e-15);
     
    gausol_test_tools_choose_system
      ( trial, m_max, n_max, p_max,
        &m, &n, &A, &p, &B, &X_ref, 
        tiny, verbose
      );
   
    fprintf(stderr, "testing with m = %d  n = %d  p = %d  tiny = %24.16e\n", m, n, p, tiny);

    /* For determinant checking: */
    double det_ref = NAN;
    if ((m == n) && (m <= gausol_test_tools_det_by_enum_SIZE_MAX))
      { det_ref = gausol_test_tools_det_by_enum(m, n, A, m);
        if (verbose) { fprintf(stderr, "det(A) by enum = %24.16e\n", det_ref); }
      }
    double rms_A = gausol_test_tools_elem_RMS(m, n, A);
    if (verbose) { fprintf(stderr, "rms_A = %24.16e\n", rms_A); }

    /* Save copies of {A,B}: */
    double *AT = talloc(m*n, double);     /* Scratch copy of {A}. */
    double *BT = talloc(m*p, double);     /* Scratch copy of {B}. */
    double *X_cmp = talloc(n*p, double);  /* Computed solution. */
    double err_X_best = +INF;
    bool_t pvrows_best, pvcols_best;
    for (bool_t pvrows = FALSE; pvrows <= TRUE; pvrows++)
      { for (bool_t pvcols = FALSE; pvcols <= TRUE; pvcols++)
          {
            fprintf(stderr, "----------------------------------------------------------------------\n");
            fprintf(stderr, "pivot_rows = %c  pivot_cols = %c\n", "FT"[pvrows], "FT"[pvcols]);

            /* Make scratch copies {AT,BT} of original {A,B}: */
            for (uint32_t ij = 0; ij < m*n; ij++) { AT[ij] = A[ij]; }
            for (uint32_t ik = 0; ik < m*p; ik++) { BT[ik] = B[ik]; }

            /* Call procedure: */
            uint32_t rank;
            double det_cmp;
            gausol_solve(m, n, AT, p, BT, X_cmp, pvrows, pvcols, tiny, &rank, &det_cmp);
            if (verbose) 
              { fprintf(stderr, "  obtained rank = %d\n", rank);
                if (rank < m) { fprintf(stderr, "  rows deficit = %d\n", m - rank); }
                if (rank < n) { fprintf(stderr, "  cols deficit = %d\n", n - rank); }
                fprintf(stderr, "\n");
                fprintf(stderr, "  obtained det(A) = %24.16e\n", det_cmp);
                gausol_print_array(stderr, 4, "%12.6f", "solution {X}:", n,NULL,0, p,NULL,0, "X",X_cmp, "");
              }

            /* Check result: */
            double err_X = gausol_test_tools_check_solve(m, n, A, p, B, X_ref, X_cmp, rank, verbose);
            if (! isnan(err_X))
              { if (verbose) 
                  { fprintf(stderr, "    pivot %c %c", "FT"[pvrows], "FT"[pvcols]);
                    fprintf(stderr, "  RMS error = %24.16e\n", err_X); 
                  }
                if (err_X < err_X_best)
                  { err_X_best = err_X; pvrows_best = pvrows; pvcols_best = pvcols; }
              }

            if (! isnan(det_ref))
              { gausol_test_tools_compare_determinants(m, n, rank, rms_A, det_ref, det_cmp, verbose); }

            fprintf(stderr, "----------------------------------------------------------------------\n");
          }
      }
    if (isfinite(err_X_best))
      { fprintf(stderr, "best solution: pivot rows = %c cols = %c", "FT"[pvrows_best], "FT"[pvcols_best]);
        fprintf(stderr, "  RMS error = %24.16e\n", err_X_best);
      }
    else
      { fprintf(stderr, "no unique solutions found with any pivoting\n"); }
    
    free(A); free(B); free(AT); free(BT);
    free(X_ref); free(X_cmp);
    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }
