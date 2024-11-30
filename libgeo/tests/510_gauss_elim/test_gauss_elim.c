/* test_gauss_elim.c --- test program for gauss_elim.h  */
/* Last edited on 2024-11-25 03:23:48 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <rn.h>
#include <rmxn.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <bool.h>
#include <affirm.h>

#include <qmin_simplex.h>
#include <gauss_elim.h>

int32_t main (int32_t argc, char **argv);

void test_gauss_elim_print_array (uint32_t trial, bool_t verbose);
  /* Tests {gauss_elim_triangularize} and {gauss_elim_triangular_det}. */

void test_gauss_elim_print_system (uint32_t trial, bool_t verbose);
  /* Tests {gauss_elim_triangularize} and {gauss_elim_triangular_det}. */

void throw_quadratic_fn(uint32_t n, double A[], double b[], bool_t verbose);
  /* Generates a random {n × n} coefficient matrix {A}, positive semidefinite,
    and a random {n}-vector {b}. */

/* CHECKING GAUSSIAN ELIMINATION STEPS

  In the following procedures, the parameters {Bmax[k]}, for {k} in
  {0..p-1}, should be the largest absolute value of any element on
  column {k} of the {B} array. */

void check_satisfaction
  ( uint32_t m,
    uint32_t n,
    double M[],
    uint32_t p,
    double X_ref[],
    double Bmax[]
  );
  /* Checks whether the {(n-p) × p} matrix {X} is still a solution of the
    system {A X = B}  where {A} is the first {n-p} columns of {M},
    and {B} is the last {p} columns. */

void check_extract_solution
  ( uint32_t m,
    uint32_t n,
    double M[],
    uint32_t p,
    double X[],
    uint32_t rank_ext,
    uint32_t lead[],
    double X_ref[],
    double Bmax[]
  );
  /* Checks the outcome of {gauss_elim_extract_solution}. The parameter {X}
    should be the solution computed by {gauss_elim_extract_solution},
    {X_ref} should be the `true' solution, and {r} should be its
    return value. */

void check_solve
  ( uint32_t m, uint32_t n, double A[],
    uint32_t p, double B[],
    double X[], uint32_t r,
    double Bmax[]
  );
  /* Checks the outcome of {gauss_elim_solve}. The parameter {X} should be
    the solution computed by {gauss_elim_solve}, and {r} should be its
    return value. */

void check_residual
  ( uint32_t m, uint32_t n, double A[],
    uint32_t p, double B[],
    double X[],
    double R[],
    double Bmax[]
  );
  /* Checks the outcome of {gauss_elim_residual}. The parameter {R} should be
    the residual computed by {gauss_elim_residual} for the given {A}, {B},
    and {X}. */

void check_solution_with_reference
  ( uint32_t m, uint32_t n, uint32_t p,
    double X[], double X_ref[],
    double Bmax[]
  );
  /* Checks a putative solution {X} of a system {A X = B} against the
    `true' solution {X_ref}; where {A} is {m × n}, {X} and {X_ref} are
    {n × p}, and {B} is {m × p}. */

void check_quadratic_min(uint32_t n, double A[], double b[], double x[], double bmax);
  /* Checks the putative solution of the constrained
    problem of minimizing {Q(x) = x' A x - 2 x'b + c} subject to
    {x[i] >= 0} for all {i}. */

double max_residual_roundoff(uint32_t m, uint32_t n, double bmax);
  /* Estimates the maximum roundoff error in the computed residual
    {A x - b} with {m} equations and {n} unknowns, given the
    absolute magnitude {bmax} of the right-hand-side elements. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    for (uint32_t i = 0;  i < MAX_RUNS; i++) { test_triangularize(i, i < 10); }
    for (uint32_t i = 0;  i < MAX_RUNS; i++) { test_gauss_elim(i, i < 10); }
    for (uint32_t i = 0;  i < MAX_RUNS; i++) { test_solve(i, i < 5); }
    for (uint32_t i = 0;  i < MAX_RUNS; i++) { test_quadratic_min(i, i < 30); }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_triangularize (uint32_t trial, bool_t verbose)
  {
    srand(20230225 + 2*trial);
    srandom(20230225 + 2*trial);
    uint32_t m = uint32_abrandom(1, MAX_ROWS); /* Rows (equations). */
    uint32_t n = uint32_abrandom(1, MAX_COLS); /* Main columns (unknowns). */
    bool_t total = (rand()/8) % 2 == 0; /* Do half the tests with {total=TRUE}. */

    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);
    fprintf(stderr, "testing with m = %d  n = %d total = %d ...\n", m, n, total);

    double A[m*n]; /* Main systems matrix. */

    if (verbose) { fprintf(stderr, "  generating matrix...\n\n"); }
    throw_matrix(m, n, A, verbose);

    /* Compute the determinant of the first {min(m,n)} columns of {A}, recursively: */
    uint32_t q = (m < n ? m : n);
    if (verbose) { fprintf(stderr, "  computing 'true' determinant...\n"); }
    double detA_ref = rmxn_det_by_enum(m, n, A, q);
    if (verbose) { fprintf(stderr, "  determinant = %22.14e\n\n", detA_ref); }

    /* Some procedures below need this: */
    double Amax = rmxn_max_abs_elem(m, n, A);

    /* Call procedures and check results: */
    if (verbose) { fprintf(stderr, "  running tests...\n\n"); }

    gauss_elim_triangularize(m, n, A, total, 0.0);
    if (verbose)
      { gauss_elim_print_array(stderr, 4, "%12.6f", "gauss_elim_triangularize:", m, n,"A",A, ""); }
    uint32_t lead[m];
    check_triangularize(m, n, A, total, lead);

    /* Up to this point, the determinant of the first {n} cols should not have changed: */
    double detA = gauss_elim_triangular_det(m, n, A, q);
    if (verbose) { fprintf(stderr, "  gauss_elim_triangular_det: result = %24.16e\n\n", detA); }
    check_determinant(q, Amax, detA, detA_ref);

    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void test_gauss_elim (uint32_t trial, bool_t verbose)
  {
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);

    uint32_t m, n, p;
    double *A, *B, *X_ref;
    choose_system(trial, &m, &n, &A, &p, &B, &X_ref, verbose);

    fprintf(stderr, "testing with m = %d  n = %d  p = %d ...\n", m, n, p);

    /* Compute the determinant of the first {min(m,n)} columns of {A}, recursively: */
    uint32_t q = (m < n ? m : n);
    if (verbose) { fprintf(stderr, "  computing 'true' determinant...\n"); }
    double detA_ref = rmxn_det_by_enum(m, n, A, q);
    if (verbose) { fprintf(stderr, "  determinant = %22.14e\n\n", detA_ref); }

    /* Some procedures below need this: */
    double Amax = rmxn_max_abs_elem(m, n, A);
    double Bmax[p];
    for (uint32_t k = 0;  k < p; k++) 
      { Bmax[k] = rmxn_max_abs_elem_in_col(n, p, B, k); }

    /* Pack arrays into a single array: */
    uint32_t np = n+p;       /* Total columns of {A} and {B}. */
    double AB[m*np];    /* Combined {A} and {B} matrices, side by side. */
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++) { AB[i*np + j] = A[i*n + j]; }
        for (uint32_t k = 0;  k < p; k++) { AB[i*np + n + k] = B[i*p + k]; }
      }
    if (verbose) { gauss_elim_print_array(stderr, 4, "%12.6f", "original packed array:", m, np,"AB",AB, ""); }

    gauss_elim_triangularize(m, np, AB, TRUE, 0.0);
    if (verbose)
      { gauss_elim_print_array(stderr, 4, "%12.6f", "gauss_elim_triangularize:", m, np,"AB",AB, ""); }
    uint32_t lead[m]; /* {lead[i]} will be {Lead{M,i)}, or {np} if {+oo}. */
    check_triangularize(m, np, AB, TRUE, lead);
    /* Find if there are all-zeros equations: */
    if (verbose)
      { for (uint32_t i = 0; i < m; i++)
         { uint32_t j = lead[i];
           if (j >= n)
             { fprintf(stderr, "Lead(AB,%d) = %d%s\n", i, j, (j >= n ? " (equation is all zeros)" : "")); }
         }
      }
    /* The equations should be equivalent, so {X_ref} should still be a solution: */
    check_satisfaction(m, np, AB, p, X_ref, Bmax);

    gauss_elim_diagonalize(m, np, AB);
    if (verbose) { gauss_elim_print_array(stderr, 4, "%12.6f", "gauss_elim_diagonalize:", m, np,"AB",AB, ""); }
    /* The triangularization conditions should still hold: */
    check_triangularize(m, np, AB, TRUE, lead);
    /* Check the diagonalization conditions: */
    check_diagonalize(m, np, AB, lead);

    /* Up to this point, the determinant of the first {n} cols should not have changed: */
    double detA = gauss_elim_triangular_det(m, np, AB, q);
    if (verbose) { fprintf(stderr, "  gauss_elim_triangular_det: result = %24.16e\n\n", detA); }
    check_determinant(q, Amax, detA, detA_ref);

    /* We should be able to normalize and solve: */
    gauss_elim_normalize(m, np, AB);
    if (verbose) { gauss_elim_print_array(stderr, 4, "%12.6f", "gauss_elim_normalize:", m, np,"AB",AB, ""); }
    /* The triangularization conditions should still hold: */
     check_triangularize(m, np, AB, TRUE, lead);
    /* The diagonalization conditions should still hold: */
    check_diagonalize(m, np, AB, lead);
    /* Check the normalization condiitons: */
    check_normalize(m, np, AB, p, lead);

    /* Call procedures and check results: */
    if (verbose) { fprintf(stderr, "  running tests...\n\n"); }
    double X[n*p];      /* Computed solution. */

    uint32_t rank_ext = gauss_elim_extract_solution(m, np, AB, p, X);
    if (verbose) { gauss_elim_print_array(stderr, 4, "%12.6f", "gauss_elim_extract_solution:", n, p,"X",X, ""); }
    if (verbose) { fprintf(stderr, "  used %d out of %d equations\n", rank_ext, m); }
    check_extract_solution(m, np, AB, p, X, rank_ext, lead, X_ref, Bmax);
    if ((rank_ext == m) && (rank_ext >= n))
      { check_solution_with_reference(m, n, p, X, X_ref, Bmax); }

    free(A); free(B); free(X_ref);
    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void test_solve (uint32_t trial, bool_t verbose)
  {
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);

    uint32_t m, n, p;
    double *A, *B, *X_ref;
    choose_system(trial, &m, &n, &A, &p, &B, &X_ref, verbose);

    fprintf(stderr, "testing with m = %d  n = %d  p = %d ...\n", m, n, p);

    /* Compute the determinant of the first {max(m,n)} columns of {A}, recursively: */
    uint32_t q = (m > n ? m : n);
    if (verbose) { fprintf(stderr, "  computing 'true' determinant...\n"); }
    double detA_ref = rmxn_det_by_enum(m, n, A, q);
    if (verbose) { fprintf(stderr, "  determinant = %22.14e\n\n", detA_ref); }

    /* Some procedures below need this: */
    double Amax = rmxn_max_abs_elem(m, n, A);
    double Bmax[p];
    for (uint32_t k = 0;  k < p; k++) 
      { Bmax[k] = rmxn_max_abs_elem_in_col(n, p, B, k); }

    /* Call procedures and check results: */
    if (verbose) { fprintf(stderr, "  running tests...\n\n"); }

    /* Check the determinant: */
    double detA = gauss_elim_determinant(m, n, A, q);
    if (verbose) { fprintf(stderr, "  gauss_elim_determinant: result = %24.16e\n\n", detA); }
    check_determinant(q, Amax, detA, detA_ref);

    double X[n*p];      /* Computed solution. */
    uint32_t r = gauss_elim_solve(m, n, A, p, B, X, 0.0);
    if (verbose) { gauss_elim_print_array(stderr, 4, "%12.6f", "solution {X}:", n, p,"X",X, ""); }
    if (verbose) { fprintf(stderr, "  used %d out of %d equations\n", r, m); }
    check_solve(m, n, A, p, B, X, r, Bmax);

    if ((r == m) && (r >= n))
      { check_solution_with_reference(m, n, p, X, X_ref, Bmax); }

    double R[m*p];      /* Residual {A X - B}. */
    gauss_elim_residual(m, n, A, p, B, X, R);
    if (verbose) { gauss_elim_print_array(stderr, 4, "%12.6f", "gauss_elim_residual:", m, p,"R",R, ""); }
    check_residual(m, n, A, p, B, X, R, Bmax);

    free(A); free(B); free(X_ref);
    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void test_quadratic_min (uint32_t trial, bool_t verbose)
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
    double bmax = rmxn_max_abs_elem_in_col(n, 1, b, 0);

    /* Call procedures and check results: */
    if (verbose) { fprintf(stderr, "  running tests...\n\n"); }
    double x[n];      /* Computed solution. */

    qms_quadratic_min(n, A, b, x);
    if (verbose) { gauss_elim_print_array(stderr, 4, "%12.6f", "gauss_elim_quadratic_min:", n, 1,"x",x, ""); }
    check_quadratic_min(n, A, b, x, bmax);

    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void check_satisfaction
  ( uint32_t m,
    uint32_t n,
    double M[],
    uint32_t p,
    double X[],
    double Bmax[]
  )
  { uint32_t q = n - p;

    for (uint32_t k = 0;  k < p; k++)
      { double tol = max_residual_roundoff(m, n, Bmax[k]);
        for (uint32_t i = 0;  i < m; i++)
          { double s = 0.0;
            for (uint32_t j = 0;  j < q; j++) { s += M[i*n + j]*X[j*p + k]; }
            s -= M[i*n + q + k];
            if (fabs(s) > tol)
              { fprintf(stderr, "(A X - B)[%d,%d] = %24.16e  tol = %24.16e\n", i, k, s, tol);
                demand(FALSE, "** system is no longer satisfied");
              }
          }
      }
  }

void check_extract_solution
  ( uint32_t m,
    uint32_t np,
    double AB[],
    uint32_t p,
    double X[],
    uint32_t rank_ext,
    uint32_t lead[],
    double X_ref[],
    double Bmax[]
  )
  { /* The normalization conditions should still hold: */
    check_normalize(m, np, AB, p, X_ref, Bmax);

    /* Check the {gauss_elim_extract_solution} post-conditions: */
    double tol[p]; /* Max roundoff error for each column of {X} and {B}. */
    for (uint32_t k = 0;  k < p; k++) { tol[k] = max_residual_roundoff(m, np, Bmax[k]); }

    uint32_t n = np - p;
    uint32_t neq = 0; /* Number of non-zero equations. */
    uint32_t j = 0; /* Next candidate for pivot. */
    for (uint32_t i = 0;  i < m; i++)
      { while ((j < n) && (j < lead[i]))
          { for (uint32_t k = 0; k < p; k++)
              { double Xjk = X[j*p + k];
                if ((! isfinite(Xjk)) || (Xjk != 0.0))
                  { fprintf(stderr, "**element {X[%d,%d]} = %24.16e should be zero\n", i, j, Xij);
                    ok = FALSE;
                  }
              }
          }
        if (j < n)
          { /* {AB[i,j]} is the pivot on line {i}. */
            assert(j == lead[i]);
            double ABij = AB[i*np + jr]; /* Pivot element. */
            assert(isfinite(ABij));
            if ((! isfinite(ABij)) || (ABij == 0.0))
              { fprintf(stderr, "** pivot {A[%d,%d]}} = %24.16e is invalid or zero\n", i, j, ABij);
                ok = FALSE;
              }
            else
              { /* We have one more effective equation: */
                for (uint32_t k = 0; k < p; k++)
                  { double Aij = AB[i*np + j];
                    double Bik = AB[i*np + n + k];
                    double Xjk_ref = X_ref[j*p + k]; /* Intended solution. */
                    double Bik_exp = Aij * Xjk_ref;  /* Expected {B}. */
                    double Xjk_cmp = X[j*p + k];     /* Solution by library. */
                    double Bik_cmp = Aij * Xjk_cmp;  /* Its {B}. */
                    /* The residual must be zero or nearly so: */
                    if (fabs(Bik_exp - Bik_cmp) > tol[k])
                      { fprintf(stderr, "** A[%d,%d] = %24.16e", i, j, Aij);
                        fprintf(stderr, " B[%d,%d] = %24.16e", i, k, Bik);
                        fprintf(stderr, " (A X)[%d,%d] = %24.16e B[%d,%d] = %24.16e dif = %24.16e tol = %24.16e\n", i, k, Yik, tol[k]);
                        fprintf(stderr, "** (A X)[%d,%d] = %24.16e B[%d,%d] = %24.16e dif = %24.16e tol = %24.16e\n", i, k, Yik, tol[k]);
                        demand(FALSE, "** this residual should be zero");
                  }
              }
          }
        else
          { /* Equation was ignored, nothing to check: */ }

        /* In any case, the next row's pivot must be at least one col to the right: */
        if (jr < UINT32_MAX) { jr++; }
      }

    /* Check number of equations used. */
    if (neq != rank_ext)
      { fprintf(stderr, " number of useful equations = %d  returned by {extract_solution} = %d\n", neq, rank_ext);
        demand(FALSE, "** counts don't match");
      }
  }

void check_solution_with_reference
  ( uint32_t m, uint32_t n, uint32_t p,
    double X[], double X_ref[],
    double Bmax[]
  )
  {
    for (uint32_t k = 0;  k < p; k++)
      { double tol = max_residual_roundoff(m, n, Bmax[k]);

        for (uint32_t j = 0;  j < n; j++)
          { double xj = X[j*p + k];
            double xj_ref = X_ref[j*p + k];
            if (fabs(xj - xj_ref) > tol)
              { fprintf(stderr, "X[%d,%d] = %24.16e\n", j, k, X[j*p + k]);
                fprintf(stderr, "X_ref[%d,%d] = %24.16e\n", j, k, X_ref[j*p + k]);
                fprintf(stderr, "difference = %24.16e  tol = %24.16e\n", fabs(xj - xj_ref), tol);
                demand(FALSE, "** computed solution does not match reference sol");
              }
          }
      }
  }


void check_solve(uint32_t m, uint32_t n, double A[], uint32_t p, double B[], double X[], uint32_t r, double Bmax[])
  { if (r < m)
      { fprintf(stderr, "  excess equations, solution not checked.\n"); }
    else
      { for (uint32_t k = 0;  k < p; k++)
          { double tol = max_residual_roundoff(m, n, Bmax[k]);
            for (uint32_t i = 0;  i < m; i++)
              { /* Compute the residual {Yik = (A X - B)[i,k]}: */
                double Yik = 0.0;
                for (uint32_t j = 0;  j < n; j++) { Yik += A[i*n + j]*X[j*p + k]; }
                Yik -= B[i*p + k];
                if (fabs(Yik) > tol)
                  { fprintf(stderr, "(A X - B)[%d,%d] = %24.16e  tol = %24.16e\n", i, k, Yik, tol);
                    demand(FALSE, "** this residual should be zero");
                  }
              }
          }
      }
  }

void check_residual
  ( uint32_t m, uint32_t n, double A[],
    uint32_t p, double B[],
    double X[],
    double R[],
    double Bmax[]
  )
  {
    for (uint32_t k = 0;  k < p; k++)
      { /* The residual must be more accurate than the solution, right? */
        double tol = 0.1 * max_residual_roundoff(m, n, Bmax[k]);
        for (uint32_t i = 0;  i < m; i++)
          { /* Compute the residual {Yik = (A X - B)[i,k]}, carefully: */
            double sum = 0.0, corr = 0.0;
            for (uint32_t j = 0;  j <= n; j++)
              { double term = ( j == n ? -B[i*p + k] : A[i*n + j]*X[j*p + k] );
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            double Yik = sum;
            /* Compare with residual from library routine: */
            if (fabs(Yik - R[i*p + k]) > tol)
              { fprintf(stderr, "(A X - B)[%d,%d] = %24.16e\n", i, k, Yik);
                fprintf(stderr, "R[%d,%d] =         %24.16e\n", i, k, R[i*p + k]);
                fprintf(stderr, "tol =            %24.16e\n", tol);
                demand(FALSE, "** residual does not match");
              }
          }
      }
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
    if (verbose) { gauss_elim_print_array(stderr, 4, "%12.6f", "original solution:", n, 1,"x",x, ""); }

    /* Generate a random {q × n} basis matrix {F}: */
    uint32_t q = (n >= MAX_PHIS ? n : n + uint32_abrandom(0, MAX_PHIS - n));
    if (verbose) { fprintf(stderr, "  basis width q = %d\n", q); }
    double F[q*n]; /* Least squares basis. */
    double Fmax = n/2.718281828;
    for (uint32_t k = 0;  k < q; k++)
      { for (uint32_t i = 0;  i < n; i++)
          { F[k*n + i] = Fscale * dabrandom(-Fmax, +Fmax); }
      }
    if (verbose) { gauss_elim_print_array(stderr, 4, "%12.6f", "basis matrix:", q, n, "F",F, ""); }

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
    if (verbose) { gauss_elim_print_system(stderr, 4, "%12.6f", "original system:", n, n,"A",A, 1,"b",b, 0,NULL,NULL, ""); }
  }

double max_residual_roundoff(uint32_t m, uint32_t n, double bmax)
  { /* Assumes that the magnitude of the products {A[i,j]*x[j]}
      is comparable to the magnitude {bmax} of the right-hand-side
      vector. Also assumes that roundoff errors are independent: */
    return 1.0e-14 * bmax * sqrt(n);
  }
