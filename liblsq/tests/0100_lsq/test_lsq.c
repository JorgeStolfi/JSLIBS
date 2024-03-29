/* test_lsq --- test program for {lsq.h}  */
/* Last edited on 2022-10-20 06:31:13 by stolfi */

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
#include <jsmath.h>
#include <rn.h>
#include <rmxn.h>
#include <gauss_elim.h>

#include <lsq.h>
#include <lsq_array.h>

/* GENERAL PARAMETERS */

#define N_CASES 100
  /* Number of cases to generate. */

#define MAX_RUNS 100
  /* Max number of trials per test. */

#define MAX_VARS 10
  /* Max number of independent variables (argument coordinates per data point). */

#define MAX_FUNS 10
  /* Max number of dependent variables (function samples per data point). */

/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);

void test_lsq_fit(int32_t trial, double eps, bool_t verbose);
  /* Tests the least-squares fitter with noise magnitude {eps}. */

void test_lsq_throw_linear_fn(int32_t nx, int32_t nf, double M[]);
  /* Generates a random {nx � nf} linear transformation matrix {M}. */

void test_lsq_throw_data_point(int32_t nx, int32_t nf, double M[], double eps, double v[], double f[]);
  /* Generates a case number from some sample set,
    with {nx} independent variables (which are returned in {v[0..nx-1]})
    and {nf} dependent variables (which are returned in {f[0..nx-1]}).
    
    The independent variables will be a random point in the unit ball 
    of {R^nx}. 
    
    The dependent variables are generated by computing the
    vector-matrix product {f = v M}, and adding to each {f[i]} a
    random perturbation uniformly distributed in {[-eps _ +eps]}. */

/* CHECKING LSQ FITTING */

void test_lsq_check_fit(int32_t nx, int32_t nf, double M[], double eps, double U[]);
  /* Checks whether the solution {U[]} matches {M}, assuming that the 
    data vectors were perturbed by {�eps}. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { int32_t i;
    for (i = 0; i < MAX_RUNS; i++) 
      { double eps = pow(0.1, 2*(i % 5) + 3);
        test_lsq_fit(i, eps, i < 5);
      }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_lsq_fit(int32_t trial, double eps, bool_t verbose)
  { 
    srand(1665 + 2*trial);
    srandom(1665 + 2*trial);
    int32_t nx = rand()/(RAND_MAX/MAX_VARS) + 1; /* Number of independent variables (argument coords per data point). */
    int32_t nf = rand()/(RAND_MAX/MAX_FUNS) + 1; /* Number of dependent variables (function samples per data point). */
    int32_t nt = (int32_t)imax(10000, 2*nx + 10);    /* Number of data points. */
    
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);
    fprintf(stderr, "testing with nt = %d  nx = %d  nf = %d  eps = %g ...\n", nt, nx, nf, eps);
    
    double *M = rmxn_alloc(nx, nf); /* Correct linear map matrix. */
    
    if (verbose) { fprintf(stderr, "  generating true map...\n\n"); }
    test_lsq_throw_linear_fn(nx, nf, M);
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  true map matrix:", nx, nf, M, "\n"); }
    
    /* Data arrays for {lsq_array_fit}: */
    double *X = rmxn_alloc(nt, nx);
    double *F = rmxn_alloc(nt, nf);
    double *W = rn_alloc(nt);
    int32_t i, k;
    for (k = 0; k < nt; k++)
      { double *Xk = &(X[k*nx]);
        double *Fk = &(F[k*nf]);
        test_lsq_throw_data_point(nx, nf, M, eps, Xk, Fk);
        W[k] = drandom();
      }
    
    /* Case generator for {lsq_fit}: */
    auto void gen_data_point(int32_t k, int32_t nxg, double Xkg[], int32_t nfg, double Fkg[], double *WkgP);
    void gen_data_point(int32_t k, int32_t nxg, double Xkg[], int32_t nfg, double Fkg[], double *WkgP)
      { 
        assert(nxg == nx);
        assert(nfg == nf);
        assert((k >= 0) && (k < nt));
        double *Xk = &(X[k*nx]);
        double *Fk = &(F[k*nf]);
        for (i = 0; i < nx; i++) { Xkg[i] = Xk[i]; }
        for (i = 0; i < nf; i++) { Fkg[i] = Fk[i]; }
        (*WkgP) = W[k];
      }
    
    /* Call procedures in {lsq} and check results: */
    double *U = rmxn_alloc(nx, nf); /* Fitted linear map matrix. */
    int32_t rank; /* Rank of least squares system. */
    if (verbose) { fprintf(stderr, "  calling {lsq_fit}...\n\n"); }
    rank = lsq_fit(nt, nx, nf, gen_data_point, U, verbose);
    demand(rank == nx, "could not solve the least squares system");
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  fitted map matrix:", nx, nf, U, "\n"); }
    test_lsq_check_fit(nx, nf, M, eps, U);

    /* Call procedures in {lsq_array} and check results: */
    if (verbose) { fprintf(stderr, "  calling {lsq_array_fit}...\n\n"); }
    rank = lsq_array_fit(nt, nx, nf, X, F, W, U, verbose);
    demand(rank == nx, "could not solve the least squares system");
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  fitted map matrix:", nx, nf, U, "\n"); }
    test_lsq_check_fit(nx, nf, M, eps, U);

    /* Cleanup: */
    free(M);
    free(U);
    free(X);
    free(F);
    free(W);

    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void test_lsq_check_fit(int32_t nx, int32_t nf, double M[], double eps, double U[])
  { int32_t iv, jf;
    double tol = 10*eps;
    for (iv = 0; iv < nx; iv++)
      { for (jf = 0; jf < nf; jf++) 
          { double Mij = M[iv*nf + jf];
            double Uij = U[iv*nf + jf];
            double s = Uij - Mij;
            if (fabs(s) > tol)
              { fprintf
                  ( stderr,
                    "(U - M)[%d,%d] = %24.16e  eps = %24.16e tol = %24.16e\n", 
                    iv, jf, s, eps, tol
                  );
                demand(FALSE, "** fit is not very good");
              }
          }
      }
    fprintf(stderr, "fit is good with tolerance %23.16e\n", tol);
  }

void test_lsq_throw_linear_fn(int32_t nx, int32_t nf, double M[])
  {
    /* Generate power-of-ten scale factor: */
    double Mscale = pow(10.0, rand()/(RAND_MAX/2));
    
    /* Generate a random linear map matrix {M}: */
    int32_t iv, jf;
    for (iv = 0; iv < nx; iv++)
      { for (jf = 0; jf < nf; jf++) 
          { M[iv*nf + jf] = Mscale * (2*drandom() - 1); }
      }
  }

void test_lsq_throw_data_point
  ( int32_t nx, 
    int32_t nf, 
    double M[],
    double eps,
    double v[],
    double f[]
  )
  {
    rn_throw_ball(nx, v);
    rmxn_map_row(nx, nf, v, M, f);
    int32_t jf;
    for (jf = 0; jf < nf; jf++) 
      { f[jf] += eps * (2*drandom() - 1); }
  }
  
