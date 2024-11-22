/* test program for {neuromat_poly_fit} */
/* Last edited on 2024-11-20 05:57:36 by stolfi */

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
#include <rn.h>
#include <rmxn.h>

#include <neuromat_poly.h>

/* GENERAL PARAMETERS */

#define MAX_DEG 4
  /* Max degree of approximation. */

#define N_POINTS 100
  /* Number of cases to generate. */

#define MAX_RUNS 100
  /* Max number of trials per test. */

/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);

void test_poly_fit(int32_t trial, int32_t g, int32_t bezix, double dev_gud, double pri_bad, double dev_bad, bool_t verbose);
  /* Tests the least squares fit. Chooses a random polynomial {P} of degree {g} 
    with {throw_poly} with arguments {g,bezix} and generates data values using {throw_data} with
    arguments {dev_gud,pri_bad,dev_bad}. */

void throw_poly(int32_t g, int32_t bezix, double P[]);
  /* If {bezix} is negative, generates a random polynomial coeff vector {P[0..g]}.
    If {bezix} is non-negative, returns in {P[0..g]} the coeffs of the Bernstein
    polynomial of degree {g} and index {i} for the argument range {[-1 _ +1]}. */

void gen_args(int32_t n, double x[]);
  /* Fills {x[0..n-1]} with equally spaced values in {[-1 _ +1]}. */ 

void gen_weights(int32_t n, double x[], double w[]);
  /* Fills {w[0..n-1]} with importance weight suitable for arguments {x[0..n-1]}. */ 

void throw_data
  ( int32_t n, 
    double t[],
    double dev_gud,
    double pri_bad,
    double dev_bad,
    double y[]
  );
  /* Fills {y[0..n-1]} with perturbed versions of {t[0..n-1]}. Each data
    value {y[k]} is an inlier with probability {1-pri_bad} and an
    outlier with probability {pri_bad}. In the first case {y[k]} is
    {t[k]} plus a Gaussian error with deviation {dev_gud}, In the second
    case {y[k]} is a Gaussian error with mean 0 and deviation
    {dev_bad}. */

void print_data(char *name, int32_t trial, int32_t iter, int32_t n, double x[], double s[], double y[], char *fmt);
  /* Writes fitting data to a file called
    "out/{name}_t{TTT}_i{III}.txt", wherer {TTT} is the three-digit
    {trial} and {III} the three-difit iteration index. If {iter} is
    {-1}, omits the "_i{III}" part.
    
    Each line will have a quadruplet {k,x[k],s[k],y[k]} for {k} in
    {0..n-1}, in a {gnuplot}-complatible format. The values {x[k]},
    {s[k]}, and {y[k]} are printed with the format {fmt}. If {y} is
    null, prints "nan" on every row. */

void check_fit(int32_t n, double t[], double s[], double dev_gud, double pri_bad);
  /* Checks whether the fitted polynomial values {s[0..n-1]} match the true 
    values {t[0..n-1]}, assuming that the inliers had noise with
    deviation {dev_gud} and the outlier probability was {pri_bad}. */

int32_t randeg(void);
  /* Returns a random degree in {0..MAX_DEG}. */
  
int32_t bez_count(int32_t g);
  /* Returns the number of Bernstein polynomials with degree {g} or less. */
  
int32_t bez_deg(int32_t ibez);
  /* Assumes {ibez} is a global index of a Bernstein polynomial, that
    is, an index in the list of all Bernstein polynomials sorted by
    degree then by Bernstein index within those with same degree.
    Returns the degree of that polynomial. */
  
int32_t bez_index(int32_t g, int32_t ibez);
  /* Assumes {ibez} is a global index of a Bernstein polynomial.
    Returns the Bernstein index of that polynomial within those
    of degree {g}. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { 
    int32_t iarg = 1;
    int32_t trial = atoi(argv[iarg]); iarg++; /* Trial number, from 0. */
    
    srand(1665 + 2*trial);
    srandom(1665 + 2*trial);
    
    int32_t nbez = bez_count(MAX_DEG); /* Number of Bernstein polys with degree at most {MAX_DEG}. */

    /* For {trial} between 0 and {nbez}, we use the Bernstein polys with degree up to {MAX_DEG}: */
    int32_t g = (trial < nbez ? bez_deg(trial) : randeg()); /* Poly degree. */
    int32_t bezix = (trial < nbez ? bez_index(g, trial) : -1); /* Bernstein index. */
    
    /* For trial between 0 and {nbez} we have only inliers, then a mix of inliers/outliers: */
    double pri_bad = (trial < nbez ? 0.0 : 0.2); /* Probability of outliers. */
    double dev_gud = 0.1;       /* Deviation of inlier noise. */
    double dev_bad = 6*dev_gud; /* Deviation of outlier values. */
    
    bool_t verbose = (trial <= nbez+3);
    
    test_poly_fit(trial, g, bezix, dev_gud, pri_bad, dev_bad, verbose); 
    fclose(stderr);
    fclose(stdout);
    return (0);
  }
  
int32_t randeg(void)
  { return int32_abrandom(0, MAX_DEG); }
  
int32_t bez_count(int32_t g)
  { return (g+1)*(g+2)/2; }
  
int32_t bez_deg(int32_t ibez)
  { int32_t g = 0;
    while (ibez >= bez_count(g)) { g++; }
    return g;
  }
  
int32_t bez_index(int32_t g, int32_t ibez)
  { int32_t iblo = (g == 0 ? 0 : bez_count(g-1));
    int32_t ibhi = bez_count(g) - 1;
    demand((iblo <= ibez) && (ibez <= ibhi), "invalid degree or global index");
    return ibez - iblo;
  }

void test_poly_fit(int32_t trial, int32_t g, int32_t bezix, double dev_gud, double pri_bad, double dev_bad, bool_t verbose)
  { 
    int32_t n = N_POINTS;                      /* Num of data points. */

    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);
    fprintf(stderr, "testing with g = %d  bezix = %d  n = %d\n", g, bezix, n);
    fprintf(stderr, "  dgud = %12.5f  pbad = %9.7f  dbad = %12.5f...\n", dev_gud, pri_bad, dev_bad);
    
    /* Choose the true polynomial: */
    if (verbose) { fprintf(stderr, "  generating polynomial...\n"); }
    int32_t g1 = g+1;
    double P[g1]; /* The true polynomial. */
    throw_poly(g, bezix, P);
    if (verbose) { rn_gen_print(stderr, g1, P, "%+12.6f", "  P = ", " ", "\n"); }

    /* Generate the argument values: */
    if (verbose) { fprintf(stderr, "  choosing argument values...\n"); }
    double x[n];   /* Argument values */
    gen_args(n, x);
    
    /* Generate the true data values: */
    if (verbose) { fprintf(stderr, "  generating the true data values...\n"); }
    double t[n];   /* True data values */
    neuromat_poly_eval_multi(g, P, n, x, t);
    
    /* Generate the perturbed data values: */
    if (verbose) { fprintf(stderr, "  generating the perturbed data values...\n"); }
    double y[n];   /* Noisy data values */
    throw_data(n, t, dev_gud, pri_bad, dev_bad, y);

    print_data("init", trial, -1, n, x, t, y, "%+12.6f");
    
    /* Choose the importance weights: */
    if (verbose) { fprintf(stderr, "  choosing the importance weights...\n"); }
    double w[n];
    gen_weights(n, x, w);
    
    /* Fit the data: */
    if (verbose) { fprintf(stderr, "  computing the best polynomial fit...\n"); }
    double Q[g1];  /* Fitted polynomial. */
    if (pri_bad == 0)
      { neuromat_poly_fit(n, x, y, w, g, Q); }
    else
      { auto void report(int32_t iter, double R[], double z[]);
      
        int32_t maxiter = 12;
        neuromat_poly_fit_robust(n, x, y, w, maxiter, g, Q, &report);
        
        void report(int32_t iter, double R[], double z[])
          { double r[n];
            neuromat_poly_eval_multi(g, R, n, x, r);
            print_data("iter", trial, iter, n, x, r, z, "%+12.6f");
          }
      }
    if (verbose) { rn_gen_print(stderr, g1, Q, "%+12.6f", "  Q = ", " ", "\n"); }
    
    /* Evaluate the polynomial: */
    double s[n];  /* Fitted polynomial values. */
    neuromat_poly_eval_multi(g, Q, n, x, s);
    print_data("term", trial, -1, n, x, s, NULL, "%+12.6f");
    
    check_fit(n, s, t, dev_gud, pri_bad);

    if (verbose) { fprintf(stderr, "done.\n"); }
    fprintf(stderr, "======================================================================\n");
  }

void throw_poly(int32_t g, int32_t bezix, double P[])
  { 
    double a = -1.0;
    double b = +1.0;
    if (bezix >= 0)
      { neuromat_poly_bezier(g, bezix, a, b, P); }
    else
      { double B[g+1];
        for (int32_t i = 0; i <= g; i++)
          { /* Choose a Bézier coeff {C} in {[-1_+1]}: */
            double C = 2*drandom() - 1;
            /* Add to {P} the Bernstein polynomial with index {i}, times {C}: */
            neuromat_poly_bezier(g, i, a, b, B);
            for (int32_t j = 0; j <= g; j++) { P[j] = P[j] + C*B[j]; }
          }
      }
  }

void throw_data
  ( int32_t n, 
    double t[],
    double dev_gud,
    double pri_bad,
    double dev_bad,
    double y[]
  )
  {
    for (int32_t k = 0; k < n; k++)
      { double toss = drandom();
        if (toss >= pri_bad)
          { /* Inlier: */
            y[k] = t[k] + dev_gud*dgaussrand();
          }
        else
          { /* Outlier: */
            y[k] = dev_bad*dgaussrand();
          }
      }
  }

void check_fit(int32_t n, double t[], double s[], double dev_gud, double pri_bad)
  { double tol = dev_gud/sqrt((1 - pri_bad)*n);
    double erms = rn_dist(n, s, t); 
    fprintf(stderr, "  rms error of fitted polynomial = %12.6f expected = %12.6f\n", erms, tol);
    if (erms > 3*tol) { fprintf(stderr, "** fit is not very good\n"); }
  }
  
void gen_args(int32_t n, double x[])
  { 
    for(int32_t k = 0; k < n; k++) { x[k] = 2*(k + 0.5)/n - 1.0; }
  }
  
void gen_weights(int32_t n, double x[], double w[])
  { 
    for(int32_t k = 0; k < n; k++)
      { double xk = x[k];
        assert((xk > -1.0) && (xk < +1.0));
        w[k] = 1.0 - 0.9*xk*xk;
      }
  }
  
void print_data(char *name, int32_t trial, int32_t iter, int32_t n, double x[], double s[], double y[], char *fmt)
  { char *fname;
    if (iter >=0)
      { fname = jsprintf("out/%s_t%03d_i%03d.txt", name, trial, iter); }
    else
      { fname = jsprintf("out/%s_t%03d.txt", name, trial); }
    FILE *wr = open_write(fname, TRUE);
    for(int32_t k = 0; k < n; k++)
      { fprintf(wr, "%5d ", k);
        fprintf(wr, fmt, x[k]);
        fprintf(wr, " ");
        fprintf(wr, fmt, s[k]);
        fprintf(wr, " ");
        fprintf(wr, fmt, (y == NULL ? NAN : y[k]));
        fprintf(wr, "\n");
      }
    fclose(wr);
    free(fname);
  }
