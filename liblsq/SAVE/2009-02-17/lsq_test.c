/* lsq_test --- test program for gauss_elim.h  */
/* Last edited on 2009-02-17 18:40:15 by stolfi */

#include <lsq.h>

#include <affirm.h>
#include <bool.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <rn.h>
#include <rmxn.h>
#include <gauss_elim.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

/* GENERAL PARAMETERS */

#define N_CASES 100
  /* Number of cases to generate. */

#define MAX_RUNS 100
  /* Max number of trials per test. */

#define MAX_VARS 10
  /* Max number of independent variables. */

#define MAX_FUNS 10
  /* Max number of dependent variables. */

/* INTERNAL PROTOTYPES */

int main (int argc, char **argv);

void test_lsq_fit(int trial, bool_t verbose);
  /* Tests the least-squares fitter. */

void throw_linear_fn(int nv, int nf, double M[]);
  /* Generates a random {nv � nf} linear transformation matrix {M}. */

void throw_case
  ( int nv, 
    int nf, 
    double M[],
    double eps,
    double v[],
    double f[]
  );
  /* Generates a case number from some sample set,
    with {nv} independent variables (which are returned in {v[0..nv-1]})
    and {nf} dependent variables (which are returned in {f[0..nv-1]}).
    
    The independent variables will be a random point in the unit ball 
    of {R^nv}. 
    
    The dependent variables are generated by computing the
    vector-matrix product {f = v M}, and adding to each {f[i]} a
    random perturbation uniformly distributed in {[-eps _ +eps]}. */

/* CHECKING LSQ FITTING */

void check_lsq_fit(int nv, int nf, double M[], double eps, double U[]);
  /* Checks whether the solution {U[]} matches {M}, assuming that the 
    data vectors were perturbed by {�eps}. */

/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  { int i;
    for (i = 0; i < MAX_RUNS; i++)
      { test_lsq_fit(i, i < 5); }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_lsq_fit(int trial, bool_t verbose)
  { 
    srand(1665 + 2*trial);
    srandom(1665 + 2*trial);
    int nv = rand()/(RAND_MAX/MAX_VARS) + 1; /* Rows (independent variables). */
    int nf = rand()/(RAND_MAX/MAX_FUNS) + 1; /* Columns (functions). */
    int nt = imax(10000, 2*nv + 10); /* Sample cases. */

    double eps = 0.001;
    
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);
    fprintf(stderr, "testing with nv = %d  nf = %d  eps = %g ...\n", nv, nf, eps);
    
    double M[nv*nf]; /* Correct linear map matrix. */
    
    if (verbose) { fprintf(stderr, "  generating true map...\n\n"); }
    throw_linear_fn(nv, nf, M);
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  true map matrix:", nv, nf, M, "\n"); }

    auto void gen_case(int i, int nv, double v[], int nf, double f[]);
      /* Case generator. */
      
    void gen_case(int i, int nv, double v[], int nf, double f[])
      { throw_case(nv, nf, M, eps, v, f);  }
    
    /* Call procedures and check results: */
    if (verbose) { fprintf(stderr, "  generating best fit...\n\n"); }
    double U[nv*nf]; /* Fitted linear map matrix. */
    int rank = lsq_fit(nv, nf, nt, gen_case, U, verbose);
    demand(rank == nv, "could not solve the least squares system");
    if (verbose) { gsel_print_array(stderr, "%12.6f", "  fitted map matrix:", nv, nf, U, "\n"); }

    /* Up to this point, the determinant of the first {nf} cols should not have changed: */
    check_lsq_fit(nv, nf, M, eps, U);

    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void check_lsq_fit(int nv, int nf, double M[], double eps, double U[])
  { int iv, jf;
    double tol = 10*eps;
    for (iv = 0; iv < nv; iv++)
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
  }

void throw_linear_fn(int nv, int nf, double M[])
  {
    /* Generate power-of-ten scale factor: */
    double Mscale = pow(10.0, rand()/(RAND_MAX/2));
    
    /* Generate a random linear map matrix {M}: */
    int iv, jf;
    for (iv = 0; iv < nv; iv++)
      { for (jf = 0; jf < nf; jf++) 
          { M[iv*nf + jf] = Mscale * (2*drandom() - 1); }
      }
  }

void throw_case
  ( int nv, 
    int nf, 
    double M[],
    double eps,
    double v[],
    double f[]
  )
  {
    rn_throw_ball(nv, v);
    rmxn_map_row(nv, nf, v, M, f);
    int jf;
    for (jf = 0; jf < nf; jf++) 
      { f[jf] += eps * (2*drandom() - 1); }
  }
  
