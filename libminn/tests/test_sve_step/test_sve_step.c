/* test_sve_step --- test program for {sve_minn_step} in sve_minn.h */
/* Last edited on 2017-03-13 21:17:59 by stolfilocal */

#define _GNU_SOURCE

#include <values.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include <vec.h>
#include <rn.h>
#include <rmxn_extra.h>
#include <rmxn.h>
#include <affirm.h>
#include <jsrandom.h>
#include <bool.h>

#include <sve_minn.h>

/* GENERAL PARAMETERS */

#define MAX_RUNS 100
  /* Max number of trials per test. */

#define MAX_DIM 10
  /* Max number of corners in simplex. */

#define MAX_VARS 10
  /* Max number of variables in goal function. */

/* INTERNAL PROTOTYPES */

int main (int argc, char **argv);

void test_sve_minn_step(int m, int n, int trial, bool_t nice, bool_t verbose);
  /* Tests the routine {sve_minn_step} with some function of {n}
    variables on some {m}-dimensional affine subspace of {R^n}. If
    {nice} is TRUE, uses a nice test function and simplex, otherwise
    uses random ones. */

void pick_nice_problem(int n, double A[], double Xref[], double *CP, bool_t verbose);
void throw_problem(int n, double A[], double Xref[], double *CP, bool_t verbose);
  /* These procedures pick the parameters {Xref,A,C} of a quadratic
     function with the formula {F(x)=(x-Xref)'*A*(x-Xref)+C}.
     
     Procedure {pick_nice_problem} yields relatively nice quadratic
     function, while {throw_problem} picks one with random
     coefficients.
     
     The stationary point {XRef} will be in the cube {U^n} where
     {U==[-1_+1]}. The entries of the coefficient matrix {A} and the
     constant term {C} will be in {U}. The parameter {C} will be
     stored into {*CP}. */

void throw_simplex(int m, int n, double v[], double x[], bool_t nice, bool_t verbose);
  /* Picks a regular {m}-dimensional simplex {v} of unit
    radius in {R^n}, such that the point {x[0..n-1]} lies in the affine span of
    {v} and inside its canonical ellipsoid. They store into {v[i*n +
    j]} the coordinate {j} of vertex {i}, for {i} in {0..m} and {j} in
    {0..n-1}.  If {nice} is TRUE, chooses a rather nice simplex. */

void check_function_values(int m, int n, double A[], double Xref[], double Xcmp[], double C, double v[]);
  /* Compares the values of {F(x)=(x-xc)'*A*(x-xc) + C} where {xc==Xref} (original
    function) and {xc==Xcmp} (function translated from {Xref} to {Xcmp}).
    The comparison is made at the corners and edge midpoints of the 
    {m}-dimensional simplex in {R^n} whose corners are stord in {v}. */

void check_optimum_position(int m, double Xref[], double Xcmp[]);
  /* Compares the computed optimum {Xcmp[0..m-1]} against the 
    `true' optimum { Xref[0..m-1]} */

void sample_function
  ( int m,
    int n,
    double A[],
    double Xref[], 
    double C, 
    double v[], 
    double fv[], 
    bool_t verbose
  );
  /* Evaluates {F(v[i,j])} where {v[i,j]} is the midpoint of 
    corners {i} and {j} of the simplex {v}. Stores the values
    into {fv[0..M-1]}, as expected by {sve_minn_step}. */

double eval_function(int n, double A[], double Xref[], double C, double x[]);
  /* Evaluates the function {F(x)=(x-Xref)'*A*(x-Xref) + C} at the point {x[0..n-1]}. */

void convert_to_cartesian(int m, int n, double cm[], double v[], double Xcmp[], bool_t verbose);
  /* Converts the computed stationary point from barycentric
    coordinates {cm[0..m]}, relative to the simplex {v}, to Cartesian
    coordinates {Xcmp[0..n-1]}. */

double max_abs_elem(int n, int m, double N[]);
  /* Returns the maximum absolute value of the elements in the {n × m} matrix {N}. */

double max_abs_col_elem(int n, int m, double N[], int j);
  /* Returns the maximum absolute value of the elements in column {j} of
     the {n × m} matrix {N}. */

/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  { int i;
    for (i = 0; i < MAX_RUNS; i++) 
      { /* Choose simplex dimension {m}, num of variables {n}, niceness, and verbosity: */
        int m = (i < 10 ? i : rand()/(RAND_MAX/MAX_DIM)) + 1;
        int n = m + (i < 10 ? i % 3 : rand()/(RAND_MAX/MAX_VARS) + 1);
        if (n > MAX_VARS) { n = MAX_VARS; }
        if (n < m) { n = m; }
        bool_t nice = (i < 5);
        bool_t verbose = (i < 6);
        /* Run one test: */
        test_sve_minn_step(m, n, i, nice, verbose);
      }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_sve_minn_step(int m, int n, int trial, bool_t nice, bool_t verbose)
  { srand(1665 + 2*trial);
    srandom(1665 + 2*trial);

    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);
    fprintf(stderr, "testing with %d-dimensional simplex and %d variables\n", m, n);

    if (verbose) { fprintf(stderr, "  generating problem...\n"); }
    double A[n*n];  /* Main coefficient matrix of {F}. */
    double Xref[n]; /* True stationary point. */
    double C;       /* Value of function {F} at {Xref}. */
    if (nice)
      { pick_nice_problem(n, A, Xref, &C, verbose); }
    else
      { throw_problem(n, A, Xref, &C, verbose); }

    if (verbose) { fprintf(stderr, "  choosing the probe simplex...\n"); }
    double v[(m+1)*n]; /* Elem {v[i*n + j]} is coord {j} of simplex vertex {i}. */
    throw_simplex(m, n, v, Xref, nice, verbose);

    if (verbose) { fprintf(stderr, "  gathering the function samples...\n"); }
    int nf = (m+1)*(m+2)/2;
    double fv[nf]; /* Function values at simplex corners and edge midpoints. */
    sample_function(m, n, A, Xref, C, v, fv, verbose);
    
    if (verbose) { fprintf(stderr, "  computing the stationary point...\n"); }
    double cm[m+1];
    sve_minn_step(m, fv, cm);
    if (verbose) 
      { rn_gen_print(stderr, m+1, cm, "%12.8f", "  barycentric coords:\n  [ ", "\n    ", " ]\n"); }
    
    if (verbose) { fprintf(stderr, "  converting barycentric to Cartesian...\n"); }
    double Xcmp[n]; 
    convert_to_cartesian(m, n, cm, v, Xcmp, verbose);
    
    if (verbose) { fprintf(stderr, "  checking function implied by computed optimum...\n"); }
    check_function_values(m, n, A, Xref, Xcmp, C, v);

    if (verbose) { fprintf(stderr, "  checking position of computed optimum...\n"); }
    check_optimum_position(n, Xref, Xcmp);
    
    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void throw_problem(int n, double A[], double Xref[], double *CP, bool_t verbose)
  {
    int i, j;
    /* Generate a random stationary point {Xref} in the signed unit cube: */
    for (j = 0; j < n; j++) { Xref[j] = 2*drandom() - 1.0; }

    if (verbose) 
      { rn_gen_print(stderr, n, Xref, "%12.8f", "  true optimum:\n  [ ", "\n    ", " ]\n"); }

    /* Generate a random coefficient matrix {A}: */
    for (i = 0; i < n; i++) 
      { for (j = 0; j < n; j++)
          { A[i*n + j] = 2*drandom() - 1.0; }
      }
    if (verbose) 
      { rmxn_gen_print
          ( stderr, n, n, A, "%12.8f", 
            "  coefficient matrix:\n  [ ", "\n    ", " ]\n",
            "[ ", " ", " ]"
            );
      }
    /* Choose the function value at {Xref}: */
    (*CP) = 2*drandom() - 1.0;
    if (verbose) { fprintf(stderr, "  true stationary value = %12.8f\n", *CP); }
  }

void pick_nice_problem(int n, double A[], double Xref[], double *CP, bool_t verbose)
  {
    int i, j;
    /* Pick the stationary point {Xref}: */
    for (j = 0; j < n; j++) { Xref[j] = ((double)j+1)/((double)n+1); }

    if (verbose) 
      { rn_gen_print(stderr, n, Xref, "%12.8f", "  true optimum:\n  [ ", "\n    ", " ]\n"); }

    /* Pick a nice coefficient matrix {A}: */
    for (i = 0; i < n; i++) 
      { for (j = 0; j < n; j++)
          { A[i*n + j] = 1.0/(1 + (i-j)*(i-j)); }
      }
    if (verbose) 
      { rmxn_gen_print
          ( stderr, n, n, A, "%12.8f", 
            "  coefficient matrix:\n  [ ", "\n    ", " ]\n",
            "[ ", " ", " ]"
            );
      }
    /* Choose the function value at {Xref}: */
    (*CP) = 0.5;
    if (verbose) { fprintf(stderr, "  true stationary value = %12.8f\n", *CP); }
  }

void throw_simplex(int m, int n, double v[], double x[], bool_t nice, bool_t verbose)
  { /* Generate a regular {m}-simplex in {R^n}, with unit radius and center at origin: */
    assert(n >= m); 
    
    if (n == m)
      { /* Generate a regular {m}-simplex, centered at the origin: */
        rmxn_regular_simplex(n, v);
      }
    else if (n > m)
      { 
        /* Generate the canonical {m}-simplex in {R^n}: */
        rmxn_canonical_simplex(m, n, v);
      }

    if (verbose) 
      { rmxn_gen_print
          ( stderr, m+1, n, v, "%12.8f", 
            "  simplex vertices (raw):\n  [ ", "\n    ", " ]\n",
            "[ ", " ", " ]"
            );
      }

    double y[n]; /* Point of simplex to be matched to {x}: */
    if (nice) 
      { /* Let the reference point {y} be the simplex's center: */
        int i, j;
        for (j = 0; j < n; j++)
          { double sum = 0;
            for (i = 0; i <= m; i++) { sum += v[i*n + j]; }
            y[j] = sum/(m+1);
          }
      }
    else
      { /* Rotate the simplex in a random orientation in {R^n}: */
        rmxn_spin_rows(m, n, v, v);
        /* Pick random point {b} in enclosing ball of the canonical {m}-simplex: */
        double b[m+1];
        rmxn_throw_canonical_simplex_ball(m, b);
        /* Compute the corresponding point of {R^n} relative to {v}: */
        rmxn_map_row(m+1, n, b, v, y);
      }
    
    /* Shift {v} so that {y} coincides with {x}: */
    double dx[n];
    rn_sub(n, x, y, dx);
    rmxn_shift_rows(m+1, n, v, dx, v);

    if (verbose) 
      { rmxn_gen_print
          ( stderr, m+1, n, v, "%12.8f", 
            "  simplex vertices (shifted):\n  [ ", "\n    ", " ]\n",
            "[ ", " ", " ]"
            );
      }
  }

void sample_function
  ( int m,
    int n,
    double A[],
    double Xref[], 
    double C, 
    double v[], 
    double fv[], 
    bool_t verbose
  )
  { double x[n];
    int i0, i1;
    for (i0 = 0; i0 <= m; i0++)
      { for (i1 = 0; i1 <= i0; i1++)
          { /* Set {x[0..n-1]} to the midpoint of simplex corners {i0,i1}: */
            int j;
            for (j = 0; j < n; j++) { x[j] = (v[i0*n + j] + v[i1*n + j])/2; }
            /* Get the function's value {F(x)} at {x}: */
            double Fx = eval_function(n, A, Xref, C, x);
            if (verbose) 
              { fprintf(stderr, "  x[%d,%d] = ", i0, i1); 
                rn_gen_print(stderr, n, x, "%12.8f", "[ ", " ", " ]");
                fprintf(stderr, "  F = %14.8f\n", Fx);
              }
            /* Store {F(x)} into {fv}: */
            fv[i0*(i0+1)/2 + i1] = Fx;
          }
      }
  }

void check_function_values(int m, int n, double A[], double Xref[], double Xcmp[], double C, double v[])
  { double tol = 1.0e-6;
    double x[n];
    int i0, i1;
    for (i0 = 0; i0 <= m; i0++)
      { for (i1 = 0; i1 <= i0; i1++)
          { /* Set {x[0..n-1]} to the midpoint of cartesian simplex corners {i,j}: */
            int j;
            for (j = 0; j < n; j++) { x[j] = (v[i0*n + j] + v[i1*n + j])/2; }
            /* Get the original and new function values {F(x)} at {x}: */
            double FXref = eval_function(n, A, Xref, C, x);
            double FXcmp = eval_function(n, A, Xcmp, C, x);
            double diff = FXref - FXcmp;
            if (fabs(diff) > tol)
              { fprintf(stderr, "**\n");
                rn_gen_print(stderr, n, x, "%12.8f", "  argument vector x:\n  [ ", "\n    ", " ]\n");
                fprintf(stderr, "F(x) ref =   %22.16e\n", FXref);
                fprintf(stderr, "F(x) cmp =   %22.16e\n", FXcmp);
                fprintf(stderr, "difference = %22.16e\n", diff);
                fprintf(stderr, "tolerance =  %22.16e\n", tol);
                demand(FALSE, "** computed optimum does not match reference sol");
              }
          }
      }
  }

void check_optimum_position(int n, double Xref[], double Xcmp[])
  { double tol = 1.0e-6 * sqrt((rn_norm_sqr(n, Xref) + rn_norm_sqr(n, Xcmp))/2);
    int j; 
    for (j = 0; j < n; j++)
      { double Xcmpj = Xcmp[j];
        double Xrefj = Xref[j];
        double diffj = Xcmpj - Xrefj;
        if (fabs(diffj) > tol)
          { fprintf(stderr, "Xcmp[%d] =   %22.16e\n", j, Xcmp[j]);
            fprintf(stderr, "Xref[%d] =   %22.16e\n", j, Xref[j]);
            fprintf(stderr, "difference = %22.16e\n", diffj);
            fprintf(stderr, "tolerance =  %22.16e\n", tol);
            demand(FALSE, "** computed optimum does not match reference sol");
          }
      }
  }

double eval_function(int n, double A[], double Xref[], double C, double x[])
  { int i, j;
    double F = C;
    for (i = 0; i < n; i++)
      { for (j = 0; j < n; j++)
          { F += (x[i] - Xref[i])*A[i*n + j]*(x[j] - Xref[j]); }
      }
    return F;
  }

void convert_to_cartesian(int m, int n, double cm[], double v[], double Xcmp[], bool_t verbose)
  { int i, j;
    for (j = 0; j < n; j++)
      { double Xsum = 0;
        for (i = 0; i <= m; i++) { Xsum += cm[i] * v[i*n + j]; }
        Xcmp[j] = Xsum;
      }
    if (verbose) 
      { rn_gen_print(stderr, n, Xcmp, "%12.8f", "  computed optimum:\n  [ ", "\n    ", " ]\n"); }
  }

double max_abs_elem(int n, int m, double N[])
  { int i, j;
    double emax = 0.0;
    for (i = 0; i < n; i++)
      { for (j = 0; j < m; j++)
          { double Mij = fabs(N[i*m + j]);
            if (Mij > emax) { emax = Mij; }
          }
      }
    return emax;
  }

double max_abs_col_elem(int n, int m, double N[], int j)
  { int i;
    double emax = 0.0;
    for (i = 0; i < n; i++)
      { double Mij = fabs(N[i*m + j]);
        if (Mij > emax) { emax = Mij; }
      }
    return emax;
  }
