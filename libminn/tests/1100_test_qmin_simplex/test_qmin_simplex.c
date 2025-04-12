/* test_qmin_simplex.c --- test program for qmin_simplex.h  */
/* Last edited on 2025-03-19 13:58:10 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include <rn.h>
#include <rmxn.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <jsprintf.h>
#include <bool.h>
#include <affirm.h>

#include <qmin_simplex.h>

#include <qmin_simplex.h>

/* INTERNAL PROTOTYPES */

#define MAX_RUNS 100

int32_t main (int32_t argc, char **argv);

void test_qmin_simplex(uint32_t trial, bool_t verbose);
  /* Tests {gausol_quadratic_min}. */

void throw_quadratic_fn(uint32_t n, double A[], double b[], bool_t verbose);
  /* Generates a random {n × n} coefficient matrix {A}, positive semidefinite,
    and a random {n}-vector {b}. */

void check_quadratic_min(uint32_t n, double A[], double b[], double x[], double bmax);
  /* Checks the putative solution of the constrained
    problem of minimizing {Q(x) = x' A x - 2 x'b + c} subject to
    {x[i] >= 0} for all {i}. */

double max_abs_elem(uint32_t m, uint32_t n, double M[]);
  /* Returns the maximum absolute value of the elements in the {m × n} matrix {M}. */

double max_abs_col_elem(uint32_t m, uint32_t n, double M[], uint32_t j);
  /* Returns the maximum absolute value of the elements in column {j} of
     the {m × n} matrix {M}. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    for (uint32_t i = 0;  i < MAX_RUNS; i++) { test_triangularize(i, i < 10); }
    for (uint32_t i = 0;  i < MAX_RUNS; i++) { test_gausol(i, i < 10); }
    for (uint32_t i = 0;  i < MAX_RUNS; i++) { test_solve(i, i < 5); }
    for (uint32_t i = 0;  i < MAX_RUNS; i++) { test_qmin_simplex(i, i < 30); }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_qmin_simplex (uint32_t trial, bool_t verbose)
  {
    srand(1665 + 418*trial);
    srandom(1665 + 418*trial);
    uint32_t n = uint32_abrandom(1, MAX_COLS); /* Size of matrix. */

    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);
    fprintf(stderr, "testing with  n = %d ...\n", n);

    double A[n*n];   /* Main systems matrix. */
    double b[n];     /* Right-hand-side matrix. */

    if (verbose) { fprintf(stderr, "  generating system...\n\n"); }
    throw_quadratic_fn(n, A, b, verbose);

    /* Some procedures below need this: */
    double bmax = max_abs_col_elem(n, 1, b, 0);

    /* Call procedures and check results: */
    if (verbose) { fprintf(stderr, "  running tests...\n\n"); }
    double x[n];      /* Computed solution. */

    qms_quadratic_min(n, A, b, x);
    if (verbose) { gausol_print_array(stderr, 4, "%12.6f", "gausol_quadratic_min:", n, 1,"x",x, ""); }
    check_quadratic_min(n, A, b, x, bmax);

    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void check_quadratic_min(uint32_t n, double A[], double b[], double x[], double bmax)
  {
    double tol = max_residual_roundoff(n, n, bmax);

    for (uint32_t i = 0;  i < n; i++)
      { /* Compute the derivative {yi = (A x - b)[i]} of {Q} w.r.t {x[i]}: */
        double yi = 0.0;
        for (uint32_t j = 0;  j < n; j++) { yi += A[i*n + j]*x[j]; }
        yi -= b[i];
        /* Check conditions for {x[i]} and {(A x - b)[i]}: */
        if (x[i] < 0.0)
          { fprintf(stderr, "x[%d] = %24.16e\n", i, x[i]);
            demand(FALSE, "** this element should be non-negative");
          }
        else if (x[i] > 0.0)
          { if (fabs(yi) > tol)
              { fprintf(stderr, "x[%d]         = %24.16e\n", i, x[i]);
                fprintf(stderr, "(A x - b)[%d] = %24.16e  tol = %24.16e\n", i, yi, tol);
                demand(FALSE, "** this residual should be zero");
              }
          }
        else
          { if (yi < -tol)
              { fprintf(stderr, "x[%d]         = %24.16e\n", i, x[i]);
                fprintf(stderr, "(A x - b)[%d] = %24.16e  tol = %24.16e\n", i, yi, tol);
                demand(FALSE, "** this residual should be positive");
              }
          }
      }
  }

void throw_matrix(uint32_t m, uint32_t n, double A[], bool_t verbose)
  {
    /* Generate power-of-ten scale factors: */
    double Ascale = pow(10.0, int32_abrandom(0, 3));
    if (verbose) { fprintf(stderr, "  scale for A = %8.1e\n", Ascale); }

    /* Generate a random coefficient matrix {A}: */
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++)
         { A[i*n + j] = Ascale * dabrandom(-0.5*m*n, +0.5*m*n); }
      }

   if (verbose) { gausol_print_array(stderr, 4, "%12.6f", "  original matrix:", m, n,"A",A, ""); }
  }

#define MAX_PHIS MAX_COLS
  /* Max number of sample points in least squares basis. */

void throw_quadratic_fn(uint32_t n, double A[], double b[], bool_t verbose)
  {
    /* Generate power-of-ten scale factors: */
    double Fscale = pow(10.0, int32_abrandom(0, 2));
    double xscale = pow(10.0, int32_abrandom(0, 3));
    if (verbose) { fprintf(stderr, "  scales: F = %8.1e x = %8.1e\n", Fscale, xscale); }


    /* Generate a random goal point {x}: */
    double x[n];
    for (uint32_t j = 0;  j < n; j++)
      { x[j] = xscale * dabrandom(-0.5*n, +0.5*n); }
    if (verbose) { gausol_print_array(stderr, 4, "%12.6f", "original solution:", n, 1,"x",x, ""); }

    /* Generate a random {q × n} basis matrix {F}: */
    uint32_t q = (n >= MAX_PHIS ? n : n + uint32_abrandom(0, MAX_PHIS - n));
    if (verbose) { fprintf(stderr, "  basis width q = %d\n", q); }
    double F[q*n]; /* Least squares basis. */
    double Fmax = n/2.718281828;
    for (uint32_t k = 0;  k < q; k++)
      { for (uint32_t i = 0;  i < n; i++)
          { F[k*n + i] = Fscale * dabrandom(-Fmax, +Fmax); }
      }
    if (verbose) { gausol_print_array(stderr, 4, "%12.6f", "basis matrix:", q, n, "F",F, ""); }

    /* Compute {A = F' F: */
    for (uint32_t i = 0;  i < n; i++)
      { for (uint32_t j = 0;  j < n; j++)
          { double s = 0;
            for (uint32_t k = 0;  k < q; k++) { s += F[k*n + i]*F[k*n + j]; }
            A[i*n + j] = s;
          }
      }

    /* Compute the right-hand side {b = A x}: */
    for (uint32_t i = 0;  i < n; i++)
      { double s = 0;
        for (uint32_t j = 0;  j < n; j++) { s += A[i*n + j]*x[j]; }
        b[i] = s;
      }
    if (verbose) { gausol_print_system(stderr, 4, "%12.6f", "original system:", n, n,"A",A, 1,"b",b, 0,NULL,NULL, ""); }
  }

double max_abs_elem(uint32_t m, uint32_t n, double M[])
  {
    double emax = 0.0;
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++)
          { double Mij = fabs(M[i*n + j]);
            if (Mij > emax) { emax = Mij; }
          }
      }
    return emax;
  }

double max_abs_col_elem(uint32_t m, uint32_t n, double M[], uint32_t j)
  {
    double emax = 0.0;
    for (uint32_t i = 0;  i < m; i++)
      { double Mij = fabs(M[i*n + j]);
        if (Mij > emax) { emax = Mij; }
      }
    return emax;
  }
