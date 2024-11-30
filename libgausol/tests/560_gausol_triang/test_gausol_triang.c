/* test_gausol_tri-diag_.c --- test program for gausol_triang.h  */
/* Last edited on 2024-11-30 04:31:43 by stolfi */

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

#include <gausol_triang.h>

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

void test_gausol_triang (uint32_t trial, bool_t verbose);
  /* Tests {gausol_triang_reduce}, {gausol_triangular_det},
    {gausol_triang_diagonalize}. */

/* IMPLEMENTATIONS */
  
/* The do { .. } while is a hak to allow semicolon after it. */
#define BAD(fmt, ...) \
  do { ok = FALSE; fprintf(stderr, "** " fmt, ##__VA_ARGS__); } while (FALSE)

int32_t main (int32_t argc, char **argv)
  {
    for (uint32_t trial = 0;  trial < MAX_RUNS; trial++) 
      { bool_t verbose = (trial < 20);
        test_gausol_triang(trial, verbose);
      }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_gausol_triang (uint32_t trial, bool_t verbose)
  {
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s  trial = %d\n", __FUNCTION__, trial);

    srand(20230225 + 2*trial);
    srandom(20230225 + 2*trial);
      
    /* Choose a 'tiny number' threshold: */
    double tiny = (drandom() < 0.5 ? 1.0e-16 : 1.0e-180);
    
    if (verbose) { fprintf(stderr, "  throwing the system ...\n"); }
    uint32_t m, n, p;
    double *A = NULL, *B = NULL, *X = NULL;
    gausol_test_tools_choose_system
      ( trial,
        /*m_max*/ MAX_ROWS, /*n_max*/ MAX_COLS, /*p_max*/ MAX_PRBS,
        /*m_P*/ &m, /*n_P*/ &n, /*A_P*/ &A,
        /*p_P*/ &p, /*B_P*/ &B, 
        /*X_P*/ &X,
        tiny, verbose
      );
    
    /* Choose the pivoting options: */
    bool_t pivot_rows = (drandom() < 0.5);
    bool_t pivot_cols = (drandom() < 0.5);

    fprintf(stderr, "  m = %2d  n = %2d  p = %2d", m, n, p);
    fprintf(stderr, "  pivot_rows = %c  pivot_cols = %c  tiny = %24.16e\n", "FT"[pivot_rows], "FT"[pivot_cols], tiny);

    if ((m >= 1) && (n >= 1) && (drandom() < 0.3))
      { /* Introduce some linear dependencies: */
        uint32_t i0 = uint32_abrandom(0, m-1);
        if (verbose) { fprintf(stderr, "  making row %d of {A} dependent ...\n", i0); }
        gausol_test_tools_make_row_dependent(i0, m, n, A);
        if (m >= 2)
          { uint32_t i1 = uint32_abrandom(0, m-2);
            if (i1 >= i0) { i1++; }
            if (verbose) { fprintf(stderr, "  making row %d of {A} dependent ...\n", i1); }
            gausol_test_tools_make_row_dependent(i1, m, n, A);
          }
        if (verbose) { fprintf(stderr, "  recomputing RHS {B} ...\n"); }
        gausol_test_tools_multiply(m, n, p, A, X, B, tiny);
        if (verbose) 
          { gausol_print_system
              ( stderr, 6, "%12.6f", "modified system", 
                m,NULL,0,  n,NULL,0, "A",A,  p,"B",B,  0,NULL,NULL, ""
              ); 
          }
      } 

    /* Some procedures below need this: */
    double rms_A = gausol_test_tools_elem_RMS(m, n, A);

    double det_ref = NAN;  /* Original determinant of {A}, o rank {NAN} if not computed. */
    if ((m == n) && (m < gausol_test_tools_det_by_enum_SIZE_MAX))
      { if (verbose) { fprintf(stderr, "  computing 'true' determinant of {A} ...\n"); }
        det_ref = gausol_test_tools_det_by_enum(m, m, A, m);
        if (verbose) { fprintf(stderr, "  determinant = %22.14e\n\n", det_ref); }
        demand(isfinite(det_ref), "overflow in determinant computation");
      }
   
    if (verbose) { fprintf(stderr, "  calling {gausol_triang_reduce} ...\n"); }
    uint32_t *prow = (pivot_rows ? talloc(m, uint32_t) : NULL);
    uint32_t *pcol = (pivot_cols ? talloc(n, uint32_t) : NULL);
    double det_cmp;
    uint32_t rank;
    gausol_triang_reduce(m, prow, n, pcol, A, p, B, tiny, &rank, &det_cmp);
    
    if (verbose) { fprintf(stderr, "  returned rank rank = %d  det = %24.16e\n", rank, det_cmp); }
    gausol_test_tools_check_triang_reduce
      ( m,prow,  n,pcol,A,  p,B, tiny, rank, det_ref,det_cmp, rms_A, verbose);
      
    if (verbose) { fprintf(stderr, "  calling {gausol_triang_diagonalize} ...\n"); }
    gausol_triang_diagonalize(m,prow, n,pcol,A, p,B, rank, tiny);
    gausol_test_tools_check_diagonalize(m,prow, n,pcol,A, p,B, rank, verbose);

    free(A); free(B); free(X);
    if (prow != NULL) { free(prow); }
    if (pcol != NULL) { free(pcol); }
    if (verbose) { fprintf(stderr, "done.\n"); }
  }
