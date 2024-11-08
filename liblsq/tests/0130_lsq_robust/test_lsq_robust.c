/* test program for {lsq_robust} */
/* Last edited on 2024-11-07 00:46:57 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsmath.h>
#include <rn.h>
#include <rmxn.h>

#include <lsq.h>
#include <lsq_robust.h>

/* GENERAL PARAMETERS */

#define MAX_DEG 4
  /* Max degree of approximation. */

#define N_POINTS 100
  /* Number of cases to generate. */

#define MAX_RUNS 100
  /* Max number of trials per test. */

/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);

void test_lsq_robust_fit(int32_t g, int32_t ixpoly, double dev_gud, double pri_bad, double dev_bad, bool_t verbose);
  /* Tests the least squares fit. Chooses a random polynomial {R} of
    degree {g} with {test_lsq_robust_throw_poly} with arguments
    {g,ixpoly} and generates data values using {test_lsq_robust_throw_data}
    with arguments {dev_gud,pri_bad,dev_bad}. */

void test_lsq_robust_throw_poly(int32_t g, int32_t ixpoly, double R[]);
  /* If {ixpoly} is negative of greater than {g}, generates a random polynomial coeff vector {R[0..g]}.
    If {ixpoly} is in {0..g}, returns in {R[0..g]} the coeffs of the Bernstein
    polynomial of degree {g} and index {i} for the argument range {[-1 _ +1]}. */

void test_lsq_robust_gen_poly_bezier(int32_t g, int32_t ix, double a, double b, double R[]);
  /* Sets {R[0..g]} to the coefficients of the Bernstein-Bezier polynomial
    of degree {g} and index {ix} for the domain interval {[a_b]}. */

void test_lsq_robust_gen_args(int32_t nt, double x[]);
  /* Fills {x[0..nt-1]} with equally spaced values in {[-1 _ +1]}. */ 

double test_lsq_robust_poly_eval(int32_t g, double R[], double x);
  /* Evaluates the polynomial with coeffs {R[0..g]} at the given {x} value. */
  
void test_lsq_robust_poly_eval_multi(int32_t g, double R[], int32_t nt, double x[], double yt[]);
  /* Assumes that {x} and {yt} have {nt} elements. Sets {yt[k]} to the
    value of the polynomial with coeffs {R[0..g]} at the argument
    {x[k]}, for {k} in {0..nt-1]}. */

void test_lsq_robust_gen_weights(int32_t nt, double x[], double w[]);
  /* Fills {w[0..nt-1]} with importance weight suitable for arguments {x[0..nt-1]}. */ 

void test_lsq_robust_throw_data
  ( int32_t nt, 
    double yt[],
    double dev_gud,
    double pri_bad,
    double dev_bad,
    double yn[]
  );
  /* Fills {yn[0..nt-1]} with perturbed versions of {yt[0..nt-1]}. Each data
    value {yn[k]} is an inlier with probability {1-pri_bad} and an
    outlier with probability {pri_bad}. In the first case {yn[k]} is
    {yt[k]} plus a Gaussian error with deviation {dev_gud}, In the second
    case {yn[k]} is a Gaussian error with mean 0 and deviation
    {dev_bad}. */

void test_lsq_robust_print_data
  ( char *name,
    int32_t g,       /* Degree of target polynomial. */
    int32_t ixpoly,  /* Index of target polynomial. */
    int32_t iter,    /* Iteration count, or {-1} for the target and final states. */
    int32_t nt,      /* Number of sampling points. */
    double x[],      /* Argument values ({nt} elements). */
    double Pc[],     /* Assumed inlier probs for current iteration ({nt} elements). */ 
    double yc[],     /* Data values adjusted for current iteration ({nt} elements). */
    double yf[],     /* Values of current approximation at arguments {x[0..nt-1]}. */
    char *fmt,       /* Format for data values. */
    bool_t verbose
  );
  /* Writes fitting data to a file called
    "out/{name}_g{GG}_ix{TTT}_i{III}.txt", where {GG} is the two-digit degree {g}, 
    {TTT} is the three-digit polynomial index {ixpoly}, and {III} the three-digit 
    iteration index. 
    
    When {iter} is {-1}, expects {yf[0..nt-1]} to be the true data
    values {ty[0..nt-1]} from the target polynomial, without noise;
    {yc[0..nt-1]} to be the given data values {yn[0..nt-1]}, with noise
    and outliers; and {Pc} to be null.  In this case the file name
    will be just "out/{name}_g{GG}_ix{TTT}.txt" without the "_i{III}" part.
    
    If {iter} is non-negative, expects {yf[0..nt-1]} to be the values of
    the current approximating polynomial, without noise; {Pc[0..nt-1]}
    to be the assumed probabilities of each data value being inlier; and
    {yc[0..nt-1]} to be the adjusted data values computed from the true
    values, {yf}, and {Pc}.
    
    Each line of the file will have five values {k,x[k],yf[k],yc[k],Pc[k]} for {k} in
    {0..nt-1}, in a {gnuplot}-compatible format. The values {x[k]},
    {yf[k]}, and {yc[k]} are printed with the format {fmt}. If {yc} is
    null, prints "nan" on the corresponding column of every row; ditto for {Pc}. */

void test_lsq_robust_fill_basis_matrix(int32_t nt, double x[], int32_t g, double X[]);
  /* Assumes that {X} is an array with {nt} rows and {nv=g+1} columns, stored by rows.
    Fills row {k} of {X} with the powers of {x[k]}, namely 
    {X[k*nv+i] = x[k]^i}, for {i} in {0..g}. */

void test_lsq_robust_check_fit(int32_t nt, double yt[], double yf[], double dev_gud, double pri_bad);
  /* Checks whether the fitted polynomial values {yf[0..nt-1]} match the true 
    values {yt[0..nt-1]}, assuming that the inliers had noise with
    deviation {dev_gud} and the outlier probability was {pri_bad}. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { 
    int32_t iarg = 1;
    int32_t g = atoi(argv[iarg]); iarg++; /* Degree of polynomial. */
    int32_t ixpoly = atoi(argv[iarg]); iarg++; /* Polynomial index within degree {g}, from 0. */
    
    srand(1665 + 2*ixpoly);
    srandom(1665 + 2*ixpoly);
    
    /* For ixpoly between 0 and {nbez} we have only inliers, then a mix of inliers/outliers: */
    double pri_bad = (ixpoly <= 0 ? 0.0 : 0.2); /* Probability of outliers. */
    double dev_gud = 0.1;       /* Deviation of inlier noise. */
    double dev_bad = 6*dev_gud; /* Deviation of outlier values. */
    
    bool_t verbose = (ixpoly <= 3);
    
    test_lsq_robust_fit(g, ixpoly, dev_gud, pri_bad, dev_bad, verbose); 
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_lsq_robust_fit(int32_t g, int32_t ixpoly, double dev_gud, double pri_bad, double dev_bad, bool_t verbose)
  { 
    int32_t nt = N_POINTS;  /* Num of data points. */
    int32_t nv = g+1;       /* Num of coefficients to fit. */

    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, ixpoly);
    fprintf(stderr, "testing with g = %d  ixpoly = %d  nt = %d  nv = %d\n", g, ixpoly, nt, nv);
    fprintf(stderr, "  dgud = %12.5f  pbad = %9.7f  dbad = %12.5f\n", dev_gud, pri_bad, dev_bad);
    
    /* Choose the true polynomial: */
    if (verbose) { fprintf(stderr, "  generating polynomial...\n"); }
    double R[nv]; /* The true polynomial. */
    test_lsq_robust_throw_poly(g, ixpoly, R);
    if (verbose) { rn_gen_print(stderr, nv, R, "%+12.6f", "  R = ", " ", "\n"); }

    /* Generate the argument values: */
    if (verbose) { fprintf(stderr, "  choosing argument values...\n"); }
    double x[nt];   /* Argument values */
    test_lsq_robust_gen_args(nt, x);
    
    /* Generate the true data values: */
    if (verbose) { fprintf(stderr, "  generating the true data values...\n"); }
    double yt[nt];   /* True data values */
    test_lsq_robust_poly_eval_multi(g, R, nt, x, yt);
    
    /* Generate the perturbed data values: */
    if (verbose) { fprintf(stderr, "  generating the perturbed data values...\n"); }
    double yn[nt];   /* Noisy data values */
    test_lsq_robust_throw_data(nt, yt, dev_gud, pri_bad, dev_bad, yn);

    test_lsq_robust_print_data("init", g, ixpoly, -1, nt, x, NULL, yn, yt, "%+12.6f", verbose);
    
    /* Choose the importance weights: */
    if (verbose) { fprintf(stderr, "  choosing the importance weights...\n"); }
    double w[nt];
    test_lsq_robust_gen_weights(nt, x, w);
    
    /* Create the basis matrix {X}: */
    double *X = rmxn_alloc(nt,nv);
    test_lsq_robust_fill_basis_matrix(nt, x, g, X);
    
    /* Fit the data: */
    if (verbose) { fprintf(stderr, "  computing the best polynomial fit...\n"); }
    double Q[nv];  /* Fitted polynomial. */
    double P[nt];  /* Inlier probabilities. */
    { auto void report(int32_t iter, int32_t ntc, int32_t nvc, double Pc[], double yc[], double Qc[]);
        /* Called by {lsq_robust_fit} at each iteration.  Assumes that {Pc[0..nt-1]} are 
          the assumed inlier probs, {yc[0..nt-1]} are the adjusted function values,
          and {Qc[0..nv-1]} are the coefficients of the approximation fitted to {Pc,yc}. */

      int32_t maxiter = 12;
      lsq_robust_fit(nt, nv, X, yn, w, maxiter, Q, P, &report, verbose);

      void report(int32_t iter, int32_t ntc, int32_t nvc, double Pc[], double yc[], double Qc[])
        { assert(ntc == nt);
          assert(nvc == nv);
          double yf[nt];  /* Current approximation evaluated at sampling arguments. */
          test_lsq_robust_poly_eval_multi(g, Qc, nt, x, yf);
          test_lsq_robust_print_data("iter", g, ixpoly, iter, nt, x, Pc, yc, yf, "%+12.6f", verbose);
        }
    }
    if (verbose) { rn_gen_print(stderr, nv, Q, "%+12.6f", "  Q = ", " ", "\n"); }
    
    /* Evaluate the polynomial: */
    double yf[nt];  /* Fitted polynomial values. */
    test_lsq_robust_poly_eval_multi(g, Q, nt, x, yf);
    test_lsq_robust_print_data("term", g, ixpoly, -1, nt, x, NULL, NULL, yf, "%+12.6f", verbose);
    
    test_lsq_robust_check_fit(nt, yt, yf, dev_gud, pri_bad);

    free(X);
    if (verbose) { fprintf(stderr, "done.\n"); }
    fprintf(stderr, "======================================================================\n");
  }

void test_lsq_robust_throw_poly(int32_t g, int32_t ixpoly, double R[])
  { 
    double a = -1.0;
    double b = +1.0;
    if ((ixpoly >= 0) && (ixpoly <= g))
      { test_lsq_robust_gen_poly_bezier(g, ixpoly, a, b, R); }
    else
      { double B[g+1];
        for (int32_t j = 0; j <= g; j++) { R[j] = 0; }
        for (int32_t i = 0; i <= g; i++)
          { /* Choose a Bézier coeff {C} in {[-1_+1]}: */
            double C = 2*drandom() - 1;
            /* Add to {R} the Bernstein polynomial with index {i}, times {C}: */
            test_lsq_robust_gen_poly_bezier(g, i, a, b, B);
            for (int32_t j = 0; j <= g; j++) { R[j] = R[j] + C*B[j]; }
          }
      }
  }

void test_lsq_robust_gen_poly_bezier(int32_t g, int32_t ix, double a, double b, double R[])
  {
    demand((0 <= ix) && (ix <= g), "invalid Bezier exponent or index");
    double h = b - a;
    double ah = a/h, bh = b/h;
    R[0] = (double)comb(g,ix); /* Hopefully there is no overflow or rounding. */
    /* Multiply {R[0..0]} by {((x-a)/h)^ix} yielding {R[0..ix]}: */
    for (int32_t k = 0; k < ix; k++)
      { /* Multiply {R[0..k]} by {(x-a)/h} yielding {R[0..k+1]}: */
        R[k+1] = R[k]/h;
        for (int32_t r = k; r > 0; r--) { R[r] = R[r-1]/h - R[r]*ah; }
        R[0] = -R[0]*ah;
      }
    /* Multiply {R[0..ix]} by {((b-x)/h)^(g-ix)} yielding {R[0..g]}: */
    for (int32_t k = ix; k < g; k++)
      { /* Multiply {R[0..k]} by {(b-x)/h} yielding {R[0..k+1]}: */
        R[k+1] = -R[k]/h;
        for (int32_t r = k; r > 0; r--) { R[r] = -R[r-1]/h + R[r]*bh; }
        R[0] = R[0]*bh;
      }
  }

double test_lsq_robust_poly_eval(int32_t g, double R[], double x)
  {
    demand(g >= 0, "invalid power");
    double y = R[g];
    for (int32_t i = g-1; i >= 0; i--) { y = R[i] + x*y; }
    return y;
  }

void test_lsq_robust_poly_eval_multi(int32_t g, double R[], int32_t nt, double x[], double yt[])
  {
    demand(g >= 0, "invalid power");
    for (int32_t k = 0; k < nt; k++) 
      { yt[k] = test_lsq_robust_poly_eval(g, R, x[k]); }
  }

void test_lsq_robust_throw_data
  ( int32_t nt, 
    double yt[],
    double dev_gud,
    double pri_bad,
    double dev_bad,
    double yn[]
  )
  {
    for (int32_t k = 0; k < nt; k++)
      { double toss = drandom();
        if (toss >= pri_bad)
          { /* Inlier: */
            yn[k] = yt[k] + dev_gud*dgaussrand();
          }
        else
          { /* Outlier: */
            yn[k] = dev_bad*dgaussrand();
          }
      }
  }

void test_lsq_robust_check_fit(int32_t nt, double yt[], double yf[], double dev_gud, double pri_bad)
  { double tol = dev_gud/sqrt((1 - pri_bad)*nt);
    double erms = rn_dist(nt, yf, yt); 
    fprintf(stderr, "  rms error of fitted polynomial = %12.6f expected = %12.6f\n", erms, tol);
    if (erms > 3*tol) { fprintf(stderr, "** fit is not very good\n"); }
  }
  
void test_lsq_robust_gen_args(int32_t nt, double x[])
  { int32_t k;
    for(k = 0; k < nt; k++) { x[k] = 2*(k + 0.5)/nt - 1.0; }
  }
  
void test_lsq_robust_gen_weights(int32_t nt, double x[], double w[])
  { for(int32_t k = 0; k < nt; k++)
      { double xk = x[k];
        assert((xk > -1.0) && (xk < +1.0));
        w[k] = 1.0 - 0.9*xk*xk;
      }
  }
  
void test_lsq_robust_fill_basis_matrix(int32_t nt, double x[], int32_t g, double X[])
  {
    demand(g >= 0, "invalid power");
    int32_t nv = g + 1;
    for (int32_t  k = 0; k < nt; k++)
      { double *Xk = &(X[k*nv]);
        double p = 1;
        Xk[0] = 1;
        double xk = x[k];
        for (int32_t i = 1; i <= g; i++) { p *= xk; Xk[i] = p; }
      }
  }

void test_lsq_robust_print_data
  ( char *name,
    int32_t g,
    int32_t ixpoly, 
    int32_t iter, 
    int32_t nt, 
    double x[], 
    double Pc[],
    double yc[],
    double yf[], 
    char *fmt,
    bool_t verbose
  )
  { char *fname = NULL;
    if (iter >=0)
      { asprintf(&fname, "out/%s_g%02d_ix%03d_i%03d.txt", name, g, ixpoly, iter); }
    else
      { asprintf(&fname, "out/%s_g%02d_ix%03d.txt", name, g, ixpoly); }
    FILE *wr = open_write(fname, TRUE);
    for (int32_t iwr = 0; iwr <= 1; iwr++)
      { FILE *wri = (iwr == 0 ? wr : stderr);
        for(int32_t k = 0; k < nt; k++)
          { fprintf(wri, "%5d ", k);
            fprintf(wri, fmt, x[k]);
            fprintf(wri, " ");
            fprintf(wri, fmt, (yf == NULL ? NAN : yf[k]));
            fprintf(wri, " ");
            fprintf(wri, fmt, (yc == NULL ? NAN : yc[k]));
            fprintf(wri, " ");
            fprintf(wri, "%9.7f", (Pc == NULL ? NAN : Pc[k]));
            fprintf(wri, "\n");
          }
      }
    fclose(wr);
    free(fname);
  }
