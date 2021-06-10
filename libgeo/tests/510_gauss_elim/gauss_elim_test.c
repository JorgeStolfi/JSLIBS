/* gauss_elim_test --- test program for gauss_elim.h  */
/* Last edited on 2021-06-09 23:38:42 by jstolfi */

#define _GNU_SOURCE
#include <gauss_elim.h>
#include <qmin_simplex.h>

#include <affirm.h>
#include <bool.h>
#include <flt.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <rmxn.h>
#include <rn.h>

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

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

void test_gauss_elim(int32_t trial, bool_t verbose);
  /* Tests the low-level routines (except {gsel_solve} and {gsel_quadratic_min}). */

void test_solve(int32_t trial, bool_t verbose);
  /* Tests {gsel_solve}. */

void test_quadratic_min(int32_t trial, bool_t verbose);
  /* Tests {gsel_quadratic_min}. */

void throw_system(int32_t m, int32_t n, double A[], int32_t p, double B[], double X[], bool_t verbose);
  /* Generates a random {m × n} coefficient matrix {A} and a random {m × p}
    solution matrix {X}, then computes the {n × p} right-hand-side matrix 
    {B = A X}. */

void throw_quadratic_fn(int32_t n, double A[], double b[], bool_t verbose);
  /* Generates a random {n × n} coefficient matrix {A}, positive semidefinite,
    and a random {n}-vector {b}. */

double determinant(int32_t m, int32_t n, double M[], int32_t q);
  /* Determinant of the first {q} rows and columns of {M}, computed by the
    basic definition (sum of {q!} terms) Returns zero if {q > m} or {q > n}. */

/* CHECKING GAUSSIAN ELIMINATION STEPS

  In the following procedures, the parameters {Bmax[k]}, for {k} in
  {0..p-1}, should be the largest absolute value of any element on
  column {k} of the {B} array. */

void check_satisfaction(int32_t m, int32_t n, double M[], int32_t p, double X_ref[], double Bmax[]);
  /* Checks whether the {(n-p) × p} matrix {X} is still a solution of the 
    system {A X = B}  where {A} is the first {n-p} columns of {M},
    and {B} is the last {p} columns. */

void check_triangularize
  ( int32_t m, int32_t n, double M[], 
    int32_t p, double X_ref[], 
    bool_t total, 
    double Bmax[]
  );
  /* Checks the outcome of {gsel_triangularize}. The tests include
    {check_satisfaction(m,n,M,p,X_ref,Bmax)}. */

void check_diagonalize
  ( int32_t m, int32_t n, double M[], 
    int32_t p, double X_ref[], 
    bool_t total, 
    double Bmax[]
  );
  /* Checks the outcome of {gsel_diagonalize}. The tests include
    {check_satisfaction(m,n,M,p,X_ref,Bmax)}. */

void check_normalize
  ( int32_t m, int32_t n, double M[], 
    int32_t p, double X_ref[], 
    bool_t total, 
    double Bmax[]
  );
  /* Checks the outcome of {gsel_normalize}. The tests include
    {check_satisfaction(m,n,M,p,X_ref,Bmax)}. */

void check_extract_solution
  ( int32_t m, int32_t n, double M[], 
    int32_t p, double X[], int32_t r, 
    double X_ref[],
    bool_t total, double Bmax[]
  );
  /* Checks the outcome of {gsel_extract_solution}. The parameter {X}
    should be the solution computed by {gsel_extract_solution},
    {X_ref} should be the `true' solution, and {r} should be its
    return value. */

void check_solve
  ( int32_t m, int32_t n, double A[], 
    int32_t p, double B[], 
    double X[], int32_t r,
    double Bmax[]
  );
  /* Checks the outcome of {gsel_solve}. The parameter {X} should be
    the solution computed by {gsel_solve}, and {r} should be its
    return value. */

void check_residual
  ( int32_t m, int32_t n, double A[], 
    int32_t p, double B[], 
    double X[],
    double R[],
    double Bmax[]
  );
  /* Checks the outcome of {gsel_residual}. The parameter {R} should be
    the residual computed by {gsel_residual} for the given {A}, {B},
    and {X}. */

void check_solution_with_reference
  ( int32_t m, int32_t n, int32_t p, 
    double X[], double X_ref[], 
    double Bmax[]
  );    
  /* Checks a putative solution {X} of a system {A X = B} against the
    `true' solution {X_ref}; where {A} is {m × n}, {X} and {X_ref} are
    {n × p}, and {B} is {m × p}. */

void check_quadratic_min(int32_t n, double A[], double b[], double x[], double bmax);
  /* Checks the putative solution of the constrained
    problem of minimizing {Q(x) = x' A x - 2 x'b + c} subject to
    {x[i] >= 0} for all {i}. */

void check_determinant(int32_t m, double Amax, double detA, double detR);
  /* Checks whether the determinant {detA} of an {m × m} matrix is equal
    to {detR} except for roundoff.  Assumes that the matrix entries are
    limited to {Amax} in absolute value. */

double max_abs_elem(int32_t m, int32_t n, double M[]);
  /* Returns the maximum absolute value of the elements in the {m × n} matrix {M}. */

double max_abs_col_elem(int32_t m, int32_t n, double M[], int32_t j);
  /* Returns the maximum absolute value of the elements in column {j} of
     the {m × n} matrix {M}. */

double max_residual_roundoff(int32_t m, int32_t n, double bmax);
  /* Estimates the maximum roundoff error in the computed residual
    {A x - b} with {m} equations and {n} unknowns, given the 
    absolute magnitude {bmax} of the right-hand-side elements. */

double max_determinant_roundoff(int32_t m, double Amax);
  /* Estimates the maximum roundoff error in the computation of the
    determinant of an {m × m} matrix, given the absolute magnitude
    {bmax} of the matrix elements. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { int32_t i;
    for (i = 0; i < MAX_RUNS; i++) { test_gauss_elim(i, i < 5); }
    for (i = 0; i < MAX_RUNS; i++) { test_solve(i, i < 5); }
    for (i = 0; i < MAX_RUNS; i++) { test_quadratic_min(i, i < 30); }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_gauss_elim (int32_t trial, bool_t verbose)
  { 
    srand(1665 + 2*trial);
    srandom(1665 + 2*trial);
    int32_t m = rand()/(RAND_MAX/MAX_ROWS) + 1; /* Rows (equations). */
    int32_t n = rand()/(RAND_MAX/MAX_COLS) + 1; /* Main columns (unknowns). */
    int32_t p = rand()/(RAND_MAX/MAX_PRBS) + 1; /* RHS columns (problems). */
    bool_t total = (rand()/8) % 2 == 0; /* Do half the tests with {total=TRUE}. */

    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);
    fprintf(stderr, "testing with m = %d  n = %d  p = %d  total = %d ...\n", m, n, p, total);
    int32_t i, j, k;
    
    double A[m*n]; /* Main systems matrix. */
    double B[m*p]; /* Right-hand-side matrix. */
    double X_ref[n*p];  /* True (mostly) solution. */

    if (verbose) { fprintf(stderr, "  generating system...\n\n"); }
    throw_system(m, n, A, p, B, X_ref, verbose);

    /* Compute the determinant of the first {max(m,n)} columns of {A}, recursively: */
    int32_t q = (m > n ? m : n);
    if (verbose) { fprintf(stderr, "  computing 'true' determinant...\n"); }
    double detA_ref = determinant(m, n, A, q);
    if (verbose) { fprintf(stderr, "  determinant = %22.14e\n\n", detA_ref); }
    
    /* Some procedures below need this: */
    double Amax = max_abs_elem(m, n, A);
    double Bmax[p];
    for (k = 0; k < p; k++) { Bmax[k] = max_abs_col_elem(n, p, B, k); }

    /* Pack arrays into a single array: */
    int32_t np = n+p;       /* Total columns of {A} and {B}. */
    double AB[m*np];    /* Combined {A} and {B} matrices, side by side. */
    for (i = 0; i < m; i++) 
      { for (j = 0; j < n; j++) { AB[i*np + j] = A[i*n + j]; }
        for (k = 0; k < p; k++) { AB[i*np + n + k] = B[i*p + k]; }
      }
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  original packed array:", m, np, AB, ""); }
    
    /* Call procedures and check results: */
    if (verbose) { fprintf(stderr, "  running tests...\n\n"); }
    double X[n*p];      /* Computed solution. */

    gsel_triangularize(m, np, AB, total, 0.0);
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  gsel_triangularize:", m, np, AB, ""); }
    check_triangularize(m, np, AB, p, X_ref, total, Bmax);

    gsel_diagonalize(m, np, AB, total);
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  gsel_diagonalize:", m, np, AB, ""); }
    check_diagonalize(m, np, AB, p, X_ref, total, Bmax);

    /* Up to this point, the determinant of the first {n} cols should not have changed: */
    double detA = gsel_triangular_det(m, np, AB, q);
    if (verbose) { fprintf(stderr, "  gsel_triangular_det: result = %24.16e\n\n", detA); }
    check_determinant(q, Amax, detA, detA_ref);

    /* We should be able to normalize and solve: */
    gsel_normalize(m, np, AB, total);
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  gsel_normalize:", m, np, AB, ""); }
    check_normalize(m, np, AB, p, X_ref, total, Bmax);

    int32_t r = gsel_extract_solution(m, np, AB, p, X, total);
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  gsel_extract_solution:", n, p, X, ""); }
    if (verbose) { fprintf(stderr, "  used %d out of %d equations\n", r, m); }
    check_extract_solution(m, np, AB, p, X, r, X_ref, total, Bmax);

    if ((r == m) && (r >= n)) { check_solution_with_reference(m, n, p, X, X_ref, Bmax); }

    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void test_solve (int32_t trial, bool_t verbose)
  { 
    srand(1665 + 12*trial);
    srandom(1665 + 12*trial);
    int32_t m = rand()/(RAND_MAX/MAX_ROWS) + 1; /* Rows (equations). */
    int32_t n = rand()/(RAND_MAX/MAX_COLS) + 1; /* Main columns (unknowns). */
    int32_t p = rand()/(RAND_MAX/MAX_PRBS) + 1; /* RHS columns (problems). */
    
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);
    fprintf(stderr, "testing with m = %d  n = %d  p = %d ...\n", m, n, p);
    int32_t k;
    
    double A[m*n];      /* Main systems matrix. */
    double B[m*p];      /* Right-hand-side matrix. */
    double X_ref[n*p];  /* True (mostly) solution. */
    double R[m*p];      /* Residual {A X - B}. */

    if (verbose) { fprintf(stderr, "  generating system...\n\n"); }
    throw_system(m, n, A, p, B, X_ref, verbose);

    /* Compute the determinant of the first {max(m,n)} columns of {A}, recursively: */
    int32_t q = (m > n ? m : n);
    if (verbose) { fprintf(stderr, "  computing 'true' determinant...\n"); }
    double detA_ref = determinant(m, n, A, q);
    if (verbose) { fprintf(stderr, "  determinant = %22.14e\n\n", detA_ref); }
    
    /* Some procedures below need this: */
    double Amax = max_abs_elem(m, n, A);
    double Bmax[p];
    for (k = 0; k < p; k++) { Bmax[k] = max_abs_col_elem(n, p, B, k); }

    /* Call procedures and check results: */
    if (verbose) { fprintf(stderr, "  running tests...\n\n"); }
    double X[n*p];      /* Computed solution. */

    /* Check the determinant: */
    double detA = gsel_determinant(m, n, A, q);
    if (verbose) { fprintf(stderr, "  gsel_determinant: result = %24.16e\n\n", detA); }
    check_determinant(q, Amax, detA, detA_ref);

    int32_t r = gsel_solve(m, n, A, p, B, X, 0.0);
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  gsel_solve:", n, p, X, ""); }
    if (verbose) { fprintf(stderr, "  used %d out of %d equations\n", r, m); }
    check_solve(m, n, A, p, B, X, r, Bmax);

    if ((r == m) && (r >= n)) 
      { check_solution_with_reference(m, n, p, X, X_ref, Bmax); }

    gsel_residual(m, n, A, p, B, X, R);
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  gsel_residual:", m, p, R, ""); }
    check_residual(m, n, A, p, B, X, R, Bmax);

    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void test_quadratic_min (int32_t trial, bool_t verbose)
  { 
    srand(1665 + 418*trial);
    srandom(1665 + 418*trial);
    int32_t n = rand()/(RAND_MAX/MAX_COLS) + 1; /* Size of matrix. */

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
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  gsel_quadratic_min:", n, 1, x, ""); }
    check_quadratic_min(n, A, b, x, bmax);

    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void check_satisfaction(int32_t m, int32_t n, double M[], int32_t p, double X[], double Bmax[])
  { int32_t q = n - p;
    int32_t i, j, k;
    for (k = 0; k < p; k++) 
      { double tol = max_residual_roundoff(m, n, Bmax[k]);
        for (i = 0; i < m; i++)
          { double s = 0.0;
            for (j = 0; j < q; j++) { s += M[i*n + j]*X[j*p + k]; }
            s -= M[i*n + q + k];
            if (fabs(s) > tol)
              { fprintf(stderr, "(A X - B)[%d,%d] = %24.16e  tol = %24.16e\n", i, k, s, tol);
                demand(FALSE, "** system is no longer satisfied");
              }
          }
      }
  }

void check_triangularize
  ( int32_t m, int32_t n, double M[], 
    int32_t p, double X_ref[], 
    bool_t total, 
    double Bmax[]
  )
  { /* The equations should be equivalent, so {X} should still be a solution: */
    check_satisfaction(m, n, M, p, X_ref, Bmax);
    /* Check the {gsel_triangularize} post-conditions: */
    int32_t i, j;
    int32_t jr = -1; /* Leading non-zero column in previous row */
    for (i = 0; i < m; i++)
      { /* Find the leading nonzero column {j} in row {i}: */
        j = 0; 
        while ((j < n) && (M[i*n + j] == 0.0)) { j++; }
        if (j >= n) { j = INT32_MAX; }
        /* Check triangulation condition: */
        bool_t ok; char *msg; 
        if (total)
          { ok = (j > jr) || ((jr == INT32_MAX) && (j == INT32_MAX));
            msg = "** leading element in wrong column"; 
          }
        else
          { ok = j >= i;
            msg = "** nonzero element below diagonal";
          }
        if (!ok) 
          { fprintf(stderr, "M[%d,%d] = %24.16e\n", i, j, M[i*n + j]);
            demand(FALSE, msg);
          }
        jr = j;
      }
  }

void check_diagonalize
  ( int32_t m, int32_t n, double M[], 
    int32_t p, double X_ref[], 
    bool_t total, 
    double Bmax[]
  )
  { /* The triangularization conditions should still hold: */
    check_triangularize(m, n, M, p, X_ref, total, Bmax);
    /* Check the {gsel_diagonalize} post_conditions: */
    int32_t i, j, k;
    for (i = 0; i < m; i++)
      { if (total)
          { /* Find leading nonzero column {j} in row {i}: */
            j = 0; 
            while ((j < n) && (M[i*n + j] == 0.0)) { j++; }
            if (j >= n) { j = INT32_MAX; }
          }
        else
          { /* Aut {M[i,i]} aut nihil: */
            j = ((i < n) && (M[i*n + i] != 0.0) ? i : INT32_MAX);
          }
        if (j < INT32_MAX) 
          { /* Element {M[i,j]} must be the only non-zero elem in column {j}: */
            for (k = 0; k < m; k++) 
              { if ((k != i) && (M[k*n + j] != 0.0))
                  { fprintf(stderr, "M[%d,%d] = %24.16e\n", k, j, M[k*n + j]);
                    demand(FALSE, "** matrix is not diagonal");
                  }
              }
          }
      }
  }

void check_determinant(int32_t m, double Amax, double detA, double detR)
  {
    double tol = max_determinant_roundoff(m, Amax);
    if (fabs(detA - detR) > tol) 
      { fprintf(stderr, "determinant = %24.16e  should be = %24.16e  tol = %24.16e\n", detA, detR, tol);
        demand(FALSE, "** determinant has changed");
      }
  }

void check_normalize
  ( int32_t m, int32_t n, double M[], 
    int32_t p, double X_ref[], 
    bool_t total, 
    double Bmax[]
  )
  { /* The diagonalization conditions should still hold: */
    check_diagonalize(m, n, M, p, X_ref, total, Bmax);
    /* Check the {gsel_normalize} post-conditions: */
    int32_t i, j;
    for (i = 0; i < m; i++)
      { if (total)
          { /* Find leading nonzero column {j} in row {i}: */
            j = 0; 
            while ((j < n) && (M[i*n + j] == 0.0)) { j++; }
            if (j >= n) { j = INT32_MAX; }
          }
        else
          { /* Normalize with {total=FALSE} leaves diagonal at 1.0 or 0.0: */
            j = ((i >= n) || (M[i*n + i] == 0.0) ? INT32_MAX : i);
          }
        if (j < INT32_MAX) 
          { /* Element {M[i,j]} must be 1.0: */
            if (M[i*n + j] != 1.0)
              { fprintf(stderr, "M[%d,%d] = %24.16e\n", i, j, M[i*n + j]);
                demand(FALSE, "** matrix is not normalized");
              }
          }
      }
  }

void check_extract_solution
  ( int32_t m, int32_t n, double M[], 
    int32_t p, double X[], 
    int32_t r,double X_ref[],
    bool_t total, 
    double Bmax[]
  )
  { /* The normalization conditions should still hold: */
    check_normalize(m, n, M, p, X_ref, total, Bmax);
    /* Check the {gsel_extract_solution} post-conditions: */
    double tol[p]; /* Max roundoff error for each column of {X} and {B}. */
    int32_t k; 
    for (k = 0; k < p; k++) { tol[k] = max_residual_roundoff(m, n, Bmax[k]); }

    int32_t q = n - p;
    int32_t i, jr;
    int32_t neq = 0; /* Number of non-zero equations. */
    for (i = 0, jr = 0; i < m; i++)
      { double Mij; /* Pivot element, or {0.0} if coefs are all zero. */
        while (jr < q)
          { Mij = M[i*n + jr];
            if (Mij != 0.0) { break; }
            /* Equation ignores {X[jr,0..p-1]}, so they should have been set to zero: */
            for (k = 0; k < p; k++) 
              { double Xjk = X[jr*p + k];
                if ((Xjk != 0.0))
                  { fprintf(stderr, "X[%d,%d] = %24.16e\n", jr, k, Xjk);
                    demand(FALSE, "** this element should be zero");
                  }
              }
            if (! total) { break; }
            jr++;
          }
        if (jr >= q) { jr = INT32_MAX; Mij = 0.0; }
         
        if (Mij != 0.0)
          { /* We have one more effective equation: */
            neq++;
            /* Residuals for row {i} should be zero: */
            for (k = 0; k < p; k++)
              { /* Compute the residual {Yik = (A X - B)[i,k]}: */
                double Yik = 0.0;
                int32_t j;
                for (j = 0; j < q; j++) { Yik += M[i*n + j]*X[j*p + k]; }
                Yik -= M[i*n + q + k];

                /* The residual must be zero or nearly so: */
                if (fabs(Yik) > tol[k])
                  { fprintf(stderr, "(A X - B)[%d,%d] = %24.16e  tol = %24.16e\n", i, k, Yik, tol[k]);
                    demand(FALSE, "** this residual should be zero");
                  }
              }
          }
        else
          { /* Equation was ignored, nothing to check: */ }

        /* In any case, the next row's pivot must be at least one col to the right: */
        if (jr < INT32_MAX) { jr++; }
      }
    
    /* Check number of equations used. */
    if (neq != r)
      { fprintf(stderr, " number of useful equations = %d  used = %d\n", neq, r);
        demand(FALSE, "** counts don't match");
      }
  }

void check_solution_with_reference
  ( int32_t m, int32_t n, int32_t p, 
    double X[], double X_ref[], 
    double Bmax[]
  )
  { int32_t k; 
    for (k = 0; k < p; k++)
      { double tol = max_residual_roundoff(m, n, Bmax[k]);
        int32_t j;
        for (j = 0; j < n; j++)
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


void check_solve(int32_t m, int32_t n, double A[], int32_t p, double B[], double X[], int32_t r, double Bmax[])
  { if (r < m)
      { fprintf(stderr, "  excess equations, solution not checked.\n"); }
    else
      { int32_t i, j, k;
        for (k = 0; k < p; k++) 
          { double tol = max_residual_roundoff(m, n, Bmax[k]);
            for (i = 0; i < m; i++)
              { /* Compute the residual {Yik = (A X - B)[i,k]}: */
                double Yik = 0.0;
                for (j = 0; j < n; j++) { Yik += A[i*n + j]*X[j*p + k]; }
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
  ( int32_t m, int32_t n, double A[], 
    int32_t p, double B[], 
    double X[],
    double R[],
    double Bmax[]
  )
  { int32_t i, j, k;
    for (k = 0; k < p; k++) 
      { /* The residual must be more accurate than the solution, right? */
        double tol = 0.1 * max_residual_roundoff(m, n, Bmax[k]);
        for (i = 0; i < m; i++)
          { /* Compute the residual {Yik = (A X - B)[i,k]}, carefully: */
            double sum = 0.0, corr = 0.0;
            for (j = 0; j <= n; j++) 
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

void check_quadratic_min(int32_t n, double A[], double b[], double x[], double bmax)
  { 
    double tol = max_residual_roundoff(n, n, bmax);
    int32_t i, j;
    for (i = 0; i < n; i++)
      { /* Compute the derivative {yi = (A x - b)[i]} of {Q} w.r.t {x[i]}: */
        double yi = 0.0;
        for (j = 0; j < n; j++) { yi += A[i*n + j]*x[j]; }
        yi -= b[i];
        /* Check conditions for {x[i]} and {(A x - b)[i]}: */
        if (x[i] < 0.0)
          { fprintf(stderr, "x[%d] = %24.16e\n", j, x[i]);
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

void throw_system(int32_t m, int32_t n, double A[], int32_t p, double B[], double X[], bool_t verbose)
  {
    /* Generate power-of-ten scale factors: */
    double Ascale = pow(10.0, rand()/(RAND_MAX/3));
    double Xscale = pow(10.0, rand()/(RAND_MAX/3));
    if (verbose) { fprintf(stderr, "  scales: A = %8.1e X = %8.1e\n", Ascale, Xscale); }
    
    int32_t i, j, k;
    /* Generate a random solution matrix {X_ref}: */
    for (j = 0; j < n; j++)
      { for (k = 0; k < p; k++) 
          { X[j*p + k] = Xscale * (double)((rand() % (n*p)) - (n*p)/2); }
      }
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  original solution:", n, p, X, ""); }

    /* Generate a random coefficient matrix {A}: */
    for (i = 0; i < m; i++) 
      { for (j = 0; j < n; j++)
          { A[i*n + j] = Ascale * (double)((rand() % (n*p)) - (n*p)/2); }
      }
    /* Compute the right-hand side {B = A X}: */
    for (i = 0; i < m; i++) 
      { for (k = 0; k < p; k++) 
          { double s = 0;
            for (j = 0; j < n; j++) { s += A[i*n + j]*X[j*p + k]; }
            B[i*p + k] = s;
          }
      }
    if (verbose) { gsel_print_system(stderr, "%12.6f", "  original system:", m, n, A, p, B, ""); }
  }

#define MAX_PHIS MAX_COLS
  /* Max number of sample points in least squares basis. */

void throw_quadratic_fn(int32_t n, double A[], double b[], bool_t verbose)
  {
    /* Generate power-of-ten scale factors: */
    double Fscale = pow(10.0, int32_abrandom(0, 2));
    double xscale = pow(10.0, int32_abrandom(0, 3));
    if (verbose) { fprintf(stderr, "  scales: F = %8.1e x = %8.1e\n", Fscale, xscale); }
    
    int32_t i, j, k;
    /* Generate a random goal point {x}: */
    double x[n]; 
    for (j = 0; j < n; j++) { x[j] = xscale * (double)(int32_abrandom(0,n-1) - n/2); }
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  original solution:", n, 1, x, ""); }

    /* Generate a random {q × n} basis matrix {F}: */
    int32_t q = (n >= MAX_PHIS ? n : n + int32_abrandom(0, MAX_PHIS - n));
    double F[q*n]; /* Least squares basis. */
    for (k = 0; k < q; k++) 
      { for (i = 0; i < n; i++)
          { F[k*n + i] = Fscale * (double)(int32_abrandom(0,n-1) - n/2.718281828); }
      }
      
    /* Compute {A = F' F: */
    for (i = 0; i < n; i++) 
      { for (j = 0; j < n; j++)
          { double s = 0;
            for (k = 0; k < q; k++) { s += F[k*n + i]*F[k*n + j]; }
            A[i*n + j] = s;
          }
      }

    /* Compute the right-hand side {b = A x}: */
    for (i = 0; i < n; i++) 
      { double s = 0;
        for (j = 0; j < n; j++) { s += A[i*n + j]*x[j]; }
        b[i] = s;
      }
    if (verbose) { gsel_print_system(stderr, "%12.6f", "  original system:", n, n, A, 1, b, ""); }
  }

double determinant(int32_t m, int32_t n, double M[], int32_t q)
  { bool_t debug = FALSE;
  
    if ((q > m) || (q > n)) { return 0.0; }

    int32_t perm[q];
    int32_t i;
    for (i = 0; i < q; i++) { perm[i] = i; }
    
    auto void add_terms(int32_t k, double factor);
      /* Generates all permutations of {perm[k..q-1]}, and adds to
        {det} the {factor} times {Q[r,perm[r]} for {r} in {k..q-1}.
        Upon return, {perm[k..q-1]} are back to their initial
        state. */
        
    double det = 0.0;
    if (debug) { fprintf(stderr, "\n"); fprintf(stderr, "    det =   %22.14e\n", det); }
    add_terms(0, 1.0);
    if (debug) { fprintf(stderr, "  det =   %22.14e\n", det); }
    return det;
    
    void add_terms(int32_t k, double factor)
      { if (k >= q)
          { if (debug) 
              { fprintf(stderr, "    det(");
                int32_t j;
                char *sep = "";
                for (j = 0; j < k; j++) { fprintf(stderr, "%s%d", sep, perm[j]);  sep = ","; }
                fprintf(stderr, ")");
                fprintf(stderr, " + %22.14e\n", factor);
              }
            det += factor;
          }
        else
          { int32_t j;
            for (j = k; j < q; j++)
              { if (j != k) { int32_t t = perm[k]; perm[k] = perm[j]; perm[j] = t; }
                int32_t row = k, col = perm[k];
                double subfactor = factor * (k == j ? 1.0 : -1.0) * M[row*n + col];
                add_terms(k+1, subfactor);
                if (j != k) { int32_t t = perm[k]; perm[k] = perm[j]; perm[j] = t; }
              }
          }
      }
  }

double max_abs_elem(int32_t m, int32_t n, double M[])
  { int32_t i, j;
    double emax = 0.0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { double Mij = fabs(M[i*n + j]);
            if (Mij > emax) { emax = Mij; }
          }
      }
    return emax;
  }

double max_abs_col_elem(int32_t m, int32_t n, double M[], int32_t j)
  { int32_t i;
    double emax = 0.0;
    for (i = 0; i < m; i++)
      { double Mij = fabs(M[i*n + j]);
        if (Mij > emax) { emax = Mij; }
      }
    return emax;
  }
  
double max_residual_roundoff(int32_t m, int32_t n, double bmax)
  { /* Assumes that the magnitude of the products {A[i,j]*x[j]}
      is comparable to the magnitude {bmax} of the right-hand-side
      vector. Also assumes that roundoff errors are independent: */
    return 1.0e-14 * bmax * sqrt(n);
  }

double max_determinant_roundoff(int32_t m, double emax)
  { double nf = m; /* Number of factors in each term */
    double nt = 1.0; /* Number of terms. */
    int32_t i; for (i = 1; i <= nf; i++) { nt *= i; }
    /* Assumes that the roundoff errors are independent: */
    return 1.0e-14 * pow(emax, nf) * sqrt(nt);
  }
