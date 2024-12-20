/* test_gausol_solve.c --- test program for gausol_solve.h  */
/* Last edited on 2024-11-30 23:59:28 by stolfi */

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

void test_gausol_solve__gausol_solve_in_place(uint32_t trial, bool_t in_place, bool_t verbose);

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    for (uint32_t i = 0;  i < MAX_RUNS; i++) 
      { test_gausol_solve__gausol_solve_in_place(i, TRUE, i < 5); }
    for (uint32_t i = 0;  i < 100; i++) 
      { test_gausol_solve__gausol_solve_in_place(i, FALSE, ((i >= 10) &&(i <= 12))); }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_gausol_solve__gausol_solve_in_place (uint32_t trial, bool_t in_place, bool_t verbose)
  {
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);

    uint32_t m, n, p;
    double *A0, *B0, *X_ref;
    uint32_t m_max = MAX_ROWS;
    uint32_t n_max = MAX_COLS;
    uint32_t p_max = MAX_PRBS;
    
    double tiny = ((trial ^ 1) == 0 ? 1.0e-180 : 1.0e-15);
     
    gausol_test_tools_choose_system
      ( trial, m_max, n_max, p_max,
        &m, &n, &A0, &p, &B0, &X_ref, 
        tiny, verbose
      );
   
    fprintf(stderr, "testing with m = %d  n = %d  p = %d  tiny = %24.16e\n", m, n, p, tiny);

    /* Save copies of {A,B}: */
    double *A = talloc(m*n, double);      /* Scratch copy of {A}. */
    double *B = talloc(m*p, double);      /* Scratch copy of {B}. */
    double *X_cmp = talloc(n*p, double);  /* Computed solution. */
    double err_X_best = +INF;
    bool_t pvrows_best, pvcols_best;
    for (bool_t pvrows = FALSE; pvrows <= TRUE; pvrows++)
      { for (bool_t pvcols = FALSE; pvcols <= TRUE; pvcols++)
          {
            fprintf(stderr, "----------------------------------------------------------------------\n");
            fprintf(stderr, "in_place = %c", "FT"[in_place]);
            fprintf(stderr, "  pivot_rows = %c  pivot_cols = %c\n", "FT"[pvrows], "FT"[pvcols]);

            /* Make scratch copies {A,B} of original {A,B}: */
            for (uint32_t ij = 0; ij < m*n; ij++) { A[ij] = A0[ij]; }
            for (uint32_t ik = 0; ik < m*p; ik++) { B[ik] = B0[ik]; }

            /* Call procedure: */
            uint32_t rank;
            double det_cmp;
            if (in_place)
              { gausol_solve_in_place(m, n, A, p, B, X_cmp, pvrows, pvcols, tiny, &det_cmp, &rank); }
            else
              { gausol_solve(m, n, A, p, B, X_cmp, pvrows, pvcols, tiny, &det_cmp, &rank);
                for (uint32_t ij = 0; ij < m*n; ij++) { demand(A[ij] == A0[ij], "array {A} changed"); }
                for (uint32_t ik = 0; ik < m*p; ik++) { demand(B[ik] == B0[ik], "array {B} changed"); }
              }

            if (verbose) 
              { fprintf(stderr, "  obtained rank = %d\n", rank);
                if (rank < m) { fprintf(stderr, "  rows deficit = %d\n", m - rank); }
                if (rank < n) { fprintf(stderr, "  cols deficit = %d\n", n - rank); }
                fprintf(stderr, "\n");
                fprintf(stderr, "  obtained det(A) = %24.16e\n", det_cmp);
                gausol_print_array(stderr, 4, "%12.6f", "solution {X}:", n,NULL,0, p,NULL,0, "X",X_cmp, "");
              }

            /* Check result: */
            bool_t pivoted = (pvrows || pvcols);
            double err_X = gausol_test_tools_check_solve(m, pivoted, n, A0, p, B0, X_ref, X_cmp, rank, det_cmp, verbose);
            if (! isnan(err_X))
              { if (verbose) 
                  { fprintf(stderr, "    pivot %c %c", "FT"[pvrows], "FT"[pvcols]);
                    fprintf(stderr, "  RMS error = %24.16e\n", err_X); 
                  }
                if (err_X < err_X_best)
                  { err_X_best = err_X; pvrows_best = pvrows; pvcols_best = pvcols; }
              }

            fprintf(stderr, "----------------------------------------------------------------------\n");
          }
      }
    if (isfinite(err_X_best))
      { fprintf(stderr, "best solution: pivot rows = %c cols = %c", "FT"[pvrows_best], "FT"[pvcols_best]);
        fprintf(stderr, "  RMS error = %24.16e\n", err_X_best);
      }
    else
      { fprintf(stderr, "no unique solutions found with any pivoting\n"); }
      
    free(A0); free(B0); free(A); free(B);
    free(X_ref); free(X_cmp);
    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }
