/* test_gauss_elim.c --- test program for gauss_elim.h  */
/* Last edited on 2024-11-08 11:23:23 by stolfi */

#define _GNU_SOURCE
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
#include <bool.h>
#include <affirm.h>

#include <qmin_simplex.h>
#include <gauss_elim.h>

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

void test_triangularize (int32_t trial, bool_t verbose);
  /* Tests {gsel_triangularize} and {gsel_triangular_det}. */

void test_gauss_elim(int32_t trial, bool_t verbose);
  /* Tests the low-level routines (except {gsel_solve} and {gsel_quadratic_min}). */

void test_solve(int32_t trial, bool_t verbose);
  /* Tests {gsel_solve}. */

void test_quadratic_min(int32_t trial, bool_t verbose);
  /* Tests {gsel_quadratic_min}. */

void choose_system
  ( int32_t trial,
    int32_t *m_P,
    int32_t *n_P,
    double **A_P,
    int32_t *p_P,
    double **B_P, 
    double **X_P,
    bool_t verbose
  );
  /* Choses the problem size {m,n,p} and the arrays {A,X}, then
    computes {B = A*X}.  If {verbose}, also prints the system
    to {stderr}. */
    
void throw_matrix(int32_t m, int32_t n, double A[], bool_t verbose);
  /* Generates a random {m × n} matrix {A}. */

void throw_system(int32_t m, int32_t n, double A[], int32_t p, double X[], bool_t verbose);
  /* Generates a random {m × n} coefficient matrix {A} and a random {m × p}
    solution matrix {X}. */

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
  ( int32_t m, 
    int32_t n,
    double M[], 
    int32_t p,
    double X_ref[], 
    double Bmax[]
  );
  /* Checks the outcome of {gsel_diagonalize}. The tests include
    {check_satisfaction(m,n,M,p,X_ref,Bmax)}. */

void check_normalize
  ( int32_t m, 
    int32_t n, 
    double M[], 
    int32_t p, 
    double X_ref[], 
    double Bmax[]
  );
  /* Checks the outcome of {gsel_normalize}. The tests include
    {check_satisfaction(m,n,M,p,X_ref,Bmax)}. */

void check_extract_solution
  ( int32_t m, 
    int32_t n, 
    double M[], 
    int32_t p, 
    double X[], 
    int32_t rank_ext,
    double X_ref[], 
    double Bmax[]
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
  { 
    for (int32_t i = 0; i < MAX_RUNS; i++) { test_triangularize(i, i < 10); }
    for (int32_t i = 0; i < MAX_RUNS; i++) { test_gauss_elim(i, i < 10); }
    for (int32_t i = 0; i < MAX_RUNS; i++) { test_solve(i, i < 5); }
    for (int32_t i = 0; i < MAX_RUNS; i++) { test_quadratic_min(i, i < 30); }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_triangularize (int32_t trial, bool_t verbose)
  { 
    srand(20230225 + 2*trial);
    srandom(20230225 + 2*trial);
    int32_t m = rand()/(RAND_MAX/MAX_ROWS) + 1; /* Rows (equations). */
    int32_t n = rand()/(RAND_MAX/MAX_COLS) + 1; /* Main columns (unknowns). */
    bool_t total = (rand()/8) % 2 == 0; /* Do half the tests with {total=TRUE}. */

    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);
    fprintf(stderr, "testing with m = %d  n = %d total = %d ...\n", m, n, total);
    
    double A[m*n]; /* Main systems matrix. */

    if (verbose) { fprintf(stderr, "  generating matrix...\n\n"); }
    throw_matrix(m, n, A, verbose);

    /* Compute the determinant of the first {max(m,n)} columns of {A}, recursively: */
    int32_t q = (m > n ? m : n);
    if (verbose) { fprintf(stderr, "  computing 'true' determinant...\n"); }
    double detA_ref = determinant(m, n, A, q);
    if (verbose) { fprintf(stderr, "  determinant = %22.14e\n\n", detA_ref); }
    
    /* Some procedures below need this: */
    double Amax = max_abs_elem(m, n, A);

    /* Call procedures and check results: */
    if (verbose) { fprintf(stderr, "  running tests...\n\n"); }

    gsel_triangularize(m, n, A, total, 0.0);
    if (verbose) 
      { gsel_print_array(stderr, 4, "%12.6f", "gsel_triangularize:", m, n,"A",A, ""); }
    check_triangularize(m, n, A, 0, NULL, total, NULL);

    /* Up to this point, the determinant of the first {n} cols should not have changed: */
    double detA = gsel_triangular_det(m, n, A, q);
    if (verbose) { fprintf(stderr, "  gsel_triangular_det: result = %24.16e\n\n", detA); }
    check_determinant(q, Amax, detA, detA_ref);

    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void test_gauss_elim (int32_t trial, bool_t verbose)
  { 
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);
    
    int32_t m, n, p;
    double *A, *B, *X_ref;
    choose_system(trial, &m, &n, &A, &p, &B, &X_ref, verbose);

    fprintf(stderr, "testing with m = %d  n = %d  p = %d ...\n", m, n, p);
    
    /* Compute the determinant of the first {max(m,n)} columns of {A}, recursively: */
    int32_t q = (m > n ? m : n);
    if (verbose) { fprintf(stderr, "  computing 'true' determinant...\n"); }
    double detA_ref = determinant(m, n, A, q);
    if (verbose) { fprintf(stderr, "  determinant = %22.14e\n\n", detA_ref); }
    
    /* Some procedures below need this: */
    double Amax = max_abs_elem(m, n, A);
    double Bmax[p];
    for (int32_t k = 0; k < p; k++) { Bmax[k] = max_abs_col_elem(n, p, B, k); }

    /* Pack arrays into a single array: */
    int32_t np = n+p;       /* Total columns of {A} and {B}. */
    double AB[m*np];    /* Combined {A} and {B} matrices, side by side. */
    for (int32_t i = 0; i < m; i++) 
      { for (int32_t j = 0; j < n; j++) { AB[i*np + j] = A[i*n + j]; }
        for (int32_t k = 0; k < p; k++) { AB[i*np + n + k] = B[i*p + k]; }
      }
    if (verbose) { gsel_print_array(stderr, 4, "%12.6f", "original packed array:", m, np,"AB",AB, ""); }
    
    /* Call procedures and check results: */
    if (verbose) { fprintf(stderr, "  running tests...\n\n"); }
    double X[n*p];      /* Computed solution. */

    gsel_triangularize(m, np, AB, TRUE, 0.0);
    if (verbose) 
      { gsel_print_array(stderr, 4, "%12.6f", "gsel_triangularize:", m, np,"AB",AB, ""); }
    check_triangularize(m, np, AB, p, X_ref, TRUE, Bmax);
    /* Find if there are all-zeros equations: */
    int32_t j = 0;
    double *ABij = AB;
    for (int32_t i = 0; i < m; i++)
      { while ((j < np) && ((*ABij) == 0)) { j++; ABij++; }
        ABij += np;
      }
    if (verbose)
      { fprintf(stderr, "lead(AB,%d) = %d%s\n", m-1, j, (j >= n ? " (some eqs are zero)" : "")); }
    /* The equations should be equivalent, so {X_ref} should still be a solution: */
    check_satisfaction(m, np, AB, p, X_ref, Bmax);

    gsel_diagonalize(m, np, AB);
    if (verbose) { gsel_print_array(stderr, 4, "%12.6f", "gsel_diagonalize:", m, np,"AB",AB, ""); }
    /* The triangularization conditions should still hold: */
    check_triangularize(m, np, AB, p, X_ref, TRUE, Bmax);
    /* Check the diagonalization conditions: */
    check_diagonalize(m, np, AB, p, X_ref, Bmax);

    /* Up to this point, the determinant of the first {n} cols should not have changed: */
    double detA = gsel_triangular_det(m, np, AB, q);
    if (verbose) { fprintf(stderr, "  gsel_triangular_det: result = %24.16e\n\n", detA); }
    check_determinant(q, Amax, detA, detA_ref);

    /* We should be able to normalize and solve: */
    gsel_normalize(m, np, AB);
    if (verbose) { gsel_print_array(stderr, 4, "%12.6f", "gsel_normalize:", m, np,"AB",AB, ""); }
    /* The triangularization conditions should still hold: */
    check_triangularize(m, np, AB, p, X_ref, TRUE, Bmax);
    /* The diagonalization conditions should still hold: */
    check_diagonalize(m, np, AB, p, X_ref, Bmax);
    /* Check the normalization condiitons: */
    check_normalize(m, np, AB, p, X_ref, Bmax);

    int32_t rank_ext = gsel_extract_solution(m, np, AB, p, X);
    if (verbose) { gsel_print_array(stderr, 4, "%12.6f", "gsel_extract_solution:", n, p,"X",X, ""); }
    if (verbose) { fprintf(stderr, "  used %d out of %d equations\n", rank_ext, m); }
    check_extract_solution(m, np, AB, p, X, rank_ext, X_ref, Bmax);

    if ((rank_ext == m) && (rank_ext >= n)) { check_solution_with_reference(m, n, p, X, X_ref, Bmax); }

    free(A); free(B); free(X_ref);
    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void test_solve (int32_t trial, bool_t verbose)
  { 
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);

    int32_t m, n, p;
    double *A, *B, *X_ref;
    choose_system(trial, &m, &n, &A, &p, &B, &X_ref, verbose);

    fprintf(stderr, "testing with m = %d  n = %d  p = %d ...\n", m, n, p);

    /* Compute the determinant of the first {max(m,n)} columns of {A}, recursively: */
    int32_t q = (m > n ? m : n);
    if (verbose) { fprintf(stderr, "  computing 'true' determinant...\n"); }
    double detA_ref = determinant(m, n, A, q);
    if (verbose) { fprintf(stderr, "  determinant = %22.14e\n\n", detA_ref); }
    
    /* Some procedures below need this: */
    double Amax = max_abs_elem(m, n, A);
    double Bmax[p];
    for (int32_t k = 0; k < p; k++) { Bmax[k] = max_abs_col_elem(n, p, B, k); }

    /* Call procedures and check results: */
    if (verbose) { fprintf(stderr, "  running tests...\n\n"); }

    /* Check the determinant: */
    double detA = gsel_determinant(m, n, A, q);
    if (verbose) { fprintf(stderr, "  gsel_determinant: result = %24.16e\n\n", detA); }
    check_determinant(q, Amax, detA, detA_ref);

    double X[n*p];      /* Computed solution. */
    int32_t r = gsel_solve(m, n, A, p, B, X, 0.0);
    if (verbose) { gsel_print_array(stderr, 4, "%12.6f", "solution {X}:", n, p,"X",X, ""); }
    if (verbose) { fprintf(stderr, "  used %d out of %d equations\n", r, m); }
    check_solve(m, n, A, p, B, X, r, Bmax);

    if ((r == m) && (r >= n)) 
      { check_solution_with_reference(m, n, p, X, X_ref, Bmax); }

    double R[m*p];      /* Residual {A X - B}. */
    gsel_residual(m, n, A, p, B, X, R);
    if (verbose) { gsel_print_array(stderr, 4, "%12.6f", "gsel_residual:", m, p,"R",R, ""); }
    check_residual(m, n, A, p, B, X, R, Bmax);

    free(A); free(B); free(X_ref);
    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }
    
void choose_system
  ( int32_t trial,
    int32_t *m_P,
    int32_t *n_P,
    double **A_P,
    int32_t *p_P,
    double **B_P, 
    double **X_P,
    bool_t verbose
  ) 
  { 
    srand(1665 + 12*trial);
    srandom(1665 + 12*trial);
    
    int32_t special_trials = 4;  /* Num trial with special systems */
    
    int32_t m, n, p;
    if (trial < special_trials) 
      { m = trial+1; n = trial+1; p = 2; }
    else
      { m = rand()/(RAND_MAX/MAX_ROWS) + 1; /* Rows (equations). */
        n = rand()/(RAND_MAX/MAX_COLS) + 1; /* Main columns (unknowns). */
        p = rand()/(RAND_MAX/MAX_PRBS) + 1; /* RHS columns (problems). */
      }

    double *A = rmxn_alloc(m,n); /* Main systems matrix. */
    double *X = rmxn_alloc(n,p);  /* True (mostly) solution. */

    if (verbose) { fprintf(stderr, "  generating system...\n\n"); }
    if (trial < special_trials)
      { for (int32_t j = 0; j < n; j++)
          { for (int32_t i = 0; i < m; i++) 
              { A[i*n + j] = i + n*sin(i*j + 1); }
            for (int32_t k = 0; k < p; k++) 
              { X[j*p + k] = (k+1)*(j + cos(k*j + 1))/n; }
          }
      }
    else
      { throw_system(m, n, A, p, X, verbose); }

    /* Compute the right-hand side {B = A X}: */
    double *B = rmxn_alloc(m,p);  /* True (mostly) solution. */
    for (int32_t i = 0; i < m; i++) 
      { for (int32_t k = 0; k < p; k++) 
          { double s = 0;
            for (int32_t j = 0; j < n; j++) { s += A[i*n + j]*X[j*p + k]; }
            B[i*p + k] = s;
          }
      }
    if (verbose) { gsel_print_system(stderr, 4, "%12.6f", "original system:", m, n,"A",A, p,"B",B, -1,NULL,NULL, ""); }

    (*m_P) = m;
    (*n_P) = n;
    (*A_P) = A;
    (*p_P) = p;
    (*B_P) = B;
    (*X_P) = X;
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
    if (verbose) { gsel_print_array(stderr, 4, "%12.6f", "gsel_quadratic_min:", n, 1,"x",x, ""); }
    check_quadratic_min(n, A, b, x, bmax);

    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void check_satisfaction(int32_t m, int32_t n, double M[], int32_t p, double X[], double Bmax[])
  { int32_t q = n - p;
    
    for (int32_t k = 0; k < p; k++) 
      { double tol = max_residual_roundoff(m, n, Bmax[k]);
        for (int32_t i = 0; i < m; i++)
          { double s = 0.0;
            for (int32_t j = 0; j < q; j++) { s += M[i*n + j]*X[j*p + k]; }
            s -= M[i*n + q + k];
            if (fabs(s) > tol)
              { fprintf(stderr, "(A X - B)[%d,%d] = %24.16e  tol = %24.16e\n", i, k, s, tol);
                demand(FALSE, "** system is no longer satisfied");
              }
          }
      }
  }

void check_triangularize
  ( int32_t m, 
    int32_t n, 
    double M[], 
    int32_t p, 
    double X_ref[], 
    bool_t total, 
    double Bmax[]
  )
  { /* Check the {gsel_triangularize} post-conditions: */
    
    int32_t jr = -1; /* Leading non-zero column in previous row */
    for (int32_t i = 0; i < m; i++)
      { /* Find the leading nonzero column {j} in row {i}: */
        int32_t j = 0; 
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
  ( int32_t m, 
    int32_t n, 
    double M[], 
    int32_t p, 
    double X_ref[], 
    double Bmax[]
  )
  { /* Check the {gsel_diagonalize} post_conditions: */
    
    for (int32_t i = 0; i < m; i++)
      { int32_t jr;
        /* Find leading nonzero column {j} in row {i}: */
        jr = 0; 
        while ((jr < n) && (M[i*n + jr] == 0.0)) { jr++; }
        if (jr >= n) { jr = INT32_MAX; }
        if (jr < INT32_MAX) 
          { /* Element {M[i,jr]} must be the only non-zero elem in column {j}: */
            for (int32_t k = 0; k < m; k++) 
              { if ((k != i) && (M[k*n + jr] != 0.0))
                  { fprintf(stderr, "M[%d,%d] = %24.16e\n", k, jr, M[k*n + jr]);
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
  ( int32_t m, 
    int32_t n, 
    double M[], 
    int32_t p, 
    double X_ref[], 
    double Bmax[]
  )
  { /* Check the {gsel_normalize} post-conditions: */
    for (int32_t i = 0; i < m; i++)
      { int32_t jr;
        /* Find leading nonzero column {jr} in row {i}: */
        jr = 0; 
        while ((jr < n) && (M[i*n + jr] == 0.0)) { jr++; }
        if (jr >= n) { jr = INT32_MAX; }
        if (jr < INT32_MAX) 
          { /* Element {M[i,jr]} must be 1.0: */
            if (M[i*n + jr] != 1.0)
              { fprintf(stderr, "M[%d,%d] = %24.16e\n", i, jr, M[i*n + jr]);
                demand(FALSE, "** matrix is not normalized");
              }
          }
      }
  }

void check_extract_solution
  ( int32_t m, 
    int32_t np, 
    double AB[], 
    int32_t p, 
    double X[], 
    int32_t rank_ext,
    double X_ref[],
    double Bmax[]
  )
  { /* The normalization conditions should still hold: */
    check_normalize(m, np, AB, p, X_ref, Bmax);
    
    /* Check the {gsel_extract_solution} post-conditions: */
    double tol[p]; /* Max roundoff error for each column of {X} and {B}. */
     
    for (int32_t k = 0; k < p; k++) { tol[k] = max_residual_roundoff(m, np, Bmax[k]); }

    int32_t n = np - p;
    int32_t  jr = 0; /* Index of leading nonzero elem in row {i}. */
    int32_t neq = 0; /* Number of non-zero equations. */
    for (int32_t i = 0; i < m; i++)
      { double ABij; /* Pivot element, or {0.0} if coefs are all zero. */
        while (jr < n)
          { ABij = AB[i*np + jr];
            assert(isfinite(ABij));
            if (ABij != 0.0) { break; }
            /* Equation ignores {X[jr,0..p-1]}, so they should have been set to zero: */
            for (int32_t k = 0; k < p; k++) 
              { double Xjk = X[jr*p + k];
                if ((Xjk != 0.0))
                  { fprintf(stderr, "X[%d,%d] = %24.16e\n", jr, k, Xjk);
                    demand(FALSE, "** this element should be zero");
                  }
              }
            jr++;
          }
        if (jr >= n) { jr = INT32_MAX; ABij = 0.0; }
        assert((jr >= n) || (ABij != 0));
         
        if (ABij != 0.0)
          { /* We have one more effective equation: */
            neq++;
            /* Residuals for row {i} should be zero: */
            for (int32_t k = 0; k < p; k++)
              { /* Compute the residual {Yik = (A X - B)[i,k]}: */
                double Yik = 0.0;
                
                for (int32_t j = 0; j < n; j++) { Yik += AB[i*np + j]*X[j*p + k]; }
                Yik -= AB[i*np + n + k];

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
    if (neq != rank_ext)
      { fprintf(stderr, " number of useful equations = %d  returned by {extract_solution} = %d\n", neq, rank_ext);
        demand(FALSE, "** counts don't match");
      }
  }

void check_solution_with_reference
  ( int32_t m, int32_t n, int32_t p, 
    double X[], double X_ref[], 
    double Bmax[]
  )
  {  
    for (int32_t k = 0; k < p; k++)
      { double tol = max_residual_roundoff(m, n, Bmax[k]);
        
        for (int32_t j = 0; j < n; j++)
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
      { for (int32_t k = 0; k < p; k++) 
          { double tol = max_residual_roundoff(m, n, Bmax[k]);
            for (int32_t i = 0; i < m; i++)
              { /* Compute the residual {Yik = (A X - B)[i,k]}: */
                double Yik = 0.0;
                for (int32_t j = 0; j < n; j++) { Yik += A[i*n + j]*X[j*p + k]; }
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
  { 
    for (int32_t k = 0; k < p; k++) 
      { /* The residual must be more accurate than the solution, right? */
        double tol = 0.1 * max_residual_roundoff(m, n, Bmax[k]);
        for (int32_t i = 0; i < m; i++)
          { /* Compute the residual {Yik = (A X - B)[i,k]}, carefully: */
            double sum = 0.0, corr = 0.0;
            for (int32_t j = 0; j <= n; j++) 
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
    
    for (int32_t i = 0; i < n; i++)
      { /* Compute the derivative {yi = (A x - b)[i]} of {Q} w.r.t {x[i]}: */
        double yi = 0.0;
        for (int32_t j = 0; j < n; j++) { yi += A[i*n + j]*x[j]; }
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

void throw_matrix(int32_t m, int32_t n, double A[], bool_t verbose)
  {
    /* Generate power-of-ten scale factors: */
    double Ascale = pow(10.0, rand()/(RAND_MAX/3));
    if (verbose) { fprintf(stderr, "  scale for A = %8.1e\n", Ascale); }
    
    /* Generate a random coefficient matrix {A}: */
    for (int32_t i = 0; i < m; i++) 
      { for (int32_t j = 0; j < n; j++)
         { A[i*n + j] = Ascale * (double)((rand() % (m*n)) - (m*n)/2); }
      }

   if (verbose) { gsel_print_array(stderr, 4, "%12.6f", "  original matrix:", m, n,"A",A, ""); }
  }

void throw_system(int32_t m, int32_t n, double A[], int32_t p, double X[], bool_t verbose)
  {
    /* Generate power-of-ten scale factors: */
    double Ascale = pow(10.0, rand()/(RAND_MAX/3));
    double Xscale = pow(10.0, rand()/(RAND_MAX/3));
    if (verbose) { fprintf(stderr, "  scales: A = %8.1e X = %8.1e\n", Ascale, Xscale); }
    
    
    /* Generate a random solution matrix {X_ref}: */
    for (int32_t j = 0; j < n; j++)
      { for (int32_t k = 0; k < p; k++) 
          { X[j*p + k] = Xscale * (double)((rand() % (n*p)) - (n*p)/2); }
      }
    if (verbose) { gsel_print_array(stderr, 4, "%12.6f", "original solution:", n, p,"X",X, ""); }

    /* Generate a random coefficient matrix {A}: */
    for (int32_t i = 0; i < m; i++) 
      { for (int32_t j = 0; j < n; j++)
          { A[i*n + j] = Ascale * (double)((rand() % (n*p)) - (n*p)/2); }
      }
  }

#define MAX_PHIS MAX_COLS
  /* Max number of sample points in least squares basis. */

void throw_quadratic_fn(int32_t n, double A[], double b[], bool_t verbose)
  {
    /* Generate power-of-ten scale factors: */
    double Fscale = pow(10.0, int32_abrandom(0, 2));
    double xscale = pow(10.0, int32_abrandom(0, 3));
    if (verbose) { fprintf(stderr, "  scales: F = %8.1e x = %8.1e\n", Fscale, xscale); }
    
    
    /* Generate a random goal point {x}: */
    double x[n]; 
    for (int32_t j = 0; j < n; j++) { x[j] = xscale * (double)(int32_abrandom(0,n-1) - n/2); }
    if (verbose) { gsel_print_array(stderr, 4, "%12.6f", "original solution:", n, 1,"x",x, ""); }

    /* Generate a random {q × n} basis matrix {F}: */
    int32_t q = (n >= MAX_PHIS ? n : n + int32_abrandom(0, MAX_PHIS - n));
    double F[q*n]; /* Least squares basis. */
    for (int32_t k = 0; k < q; k++) 
      { for (int32_t i = 0; i < n; i++)
          { F[k*n + i] = Fscale * (double)(int32_abrandom(0,n-1) - n/2.718281828); }
      }
      
    /* Compute {A = F' F: */
    for (int32_t i = 0; i < n; i++) 
      { for (int32_t j = 0; j < n; j++)
          { double s = 0;
            for (int32_t k = 0; k < q; k++) { s += F[k*n + i]*F[k*n + j]; }
            A[i*n + j] = s;
          }
      }

    /* Compute the right-hand side {b = A x}: */
    for (int32_t i = 0; i < n; i++) 
      { double s = 0;
        for (int32_t j = 0; j < n; j++) { s += A[i*n + j]*x[j]; }
        b[i] = s;
      }
    if (verbose) { gsel_print_system(stderr, 4, "%12.6f", "original system:", n, n,"A",A, 1,"b",b, -1,NULL,NULL, ""); }
  }

double determinant(int32_t m, int32_t n, double M[], int32_t q)
  { if ((q > m) || (q > n)) { return 0.0; }

    int32_t perm[q];
    
    for (int32_t i = 0; i < q; i++) { perm[i] = i; }
    
    auto void add_terms(int32_t k, double factor);
      /* Generates all permutations of {perm[k..q-1]}, and adds to
        {det} the {factor} times {Q[r,perm[r]} for {r} in {k..q-1}.
        Upon return, {perm[k..q-1]} are back to their initial
        state. */
        
    double det = 0.0;
    add_terms(0, 1.0);
    return det;
    
    void add_terms(int32_t k, double factor)
      { if (k >= q)
          { det += factor; }
        else
          { 
            for (int32_t j = k; j < q; j++)
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
  { 
    double emax = 0.0;
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { double Mij = fabs(M[i*n + j]);
            if (Mij > emax) { emax = Mij; }
          }
      }
    return emax;
  }

double max_abs_col_elem(int32_t m, int32_t n, double M[], int32_t j)
  { 
    double emax = 0.0;
    for (int32_t i = 0; i < m; i++)
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
    for (int32_t i = 1; i <= nf; i++) { nt *= i; }
    /* Assumes that the roundoff errors are independent: */
    return 1.0e-14 * pow(emax, nf) * sqrt(nt);
  }
