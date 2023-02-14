/* test_lsq --- test program for constrained {lsq_solve_system}  */
/* Last edited on 2019-12-18 22:28:20 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <rn.h>
#include <rmxn.h>
#include <rmxn_extra.h>
#include <gauss_elim.h>

#include <lsq.h>
#include <lsq_array.h>

/* GENERAL PARAMETERS */

#define N_CASES 100
  /* Number of cases to generate. */

#define MAX_RUNS 20
  /* Max number of trials per test. */

#define MAX_VARS 10
  /* Max number of independent variables (argument coordinates per data point). */

#define MAX_FUNS 3
  /* Max number of dependent variables (function samples per data point). */

/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);

void test_lsq_solve_system(int32_t trial, double tol, bool_t verbose);
  /* Tests the least-squares fitter with tolerance {tol}. */

void test_lsq_throw_linear_system
  ( int32_t nx, 
    int32_t nf,
    int32_t nc,
    bool_t scramble,
    double **AP,
    double **BP, 
    double **RP,
    double **SP,
    double **UexpP
  );
  /* Generates a linear system for {lsq_solve_system}
    with {nx} unknowns (coefficients to fit) for {nf} functions to
    fit,  {nc} constraints.
    
    Namely, allocates and returns in {*AP,*BP,*RP,*SP,*UP} the following 
    matrices:  
    
      (A}, the {nx×nx} moment matrix;
      {B}, the {nx×nf} right-hand-side matrix;
      {R}, the {nc×nx} constraint coefficient matrix;
      {S}, the {nc×nf} right-hand-side of the constraints;
      {Uexp}, the {nx×nf} expected constrained optimum.
      
    If {scramble} is false, returns a canonical problem.  If {scramble}
    is true, modifies that canonical problem by random linear mappings.
      
  */

/* CHECKING LSQ FITTING */

void test_lsq_check_solution
  ( int32_t nx, 
    int32_t nf, 
    double tol, 
    double Uexp[], 
    double Ucmp[]
  );
  /* Checks whether the computed solution matrix {Ucmp[]} matches the expected 
    solution matrix {Uexp}, within the tolerance {tol}. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { for (int32_t i = 0; i < MAX_RUNS; i++) 
      { double tol = pow(0.1, 2*(i % 5) + 2);
        test_lsq_solve_system(i, tol, i < 5);
      }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_lsq_solve_system(int32_t trial, double tol, bool_t verbose)
  { 
    srand(1665 + 2*trial);
    srandom(1665 + 2*trial);
    int32_t nx = rand()/(RAND_MAX/MAX_VARS) + 1; /* Number of indep variables. */
    int32_t nf = rand()/(RAND_MAX/MAX_FUNS) + 1; /* Number of dep variables. */
    int32_t nc = rand()/(RAND_MAX/nx);    /* Number of constraints. */
    
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);
    fprintf(stderr, "testing with nx = %d  nf = %d  nc = %d   tol = %23.16e ...\n", nx, nf, nc, tol);
    
    /* Generate a system with known solution {Uexp}: */
    double *A, *B, *R, *S, *Uexp;
    bool_t scramble = (trial > 0);
    test_lsq_throw_linear_system(nx, nf, nc, scramble, &A, &B, &R, &S, &Uexp);
    
    /* Solve the system with {lsq_solve_system}: */
    double *Ucmp = rmxn_alloc(nx, nf); /* Computed solution. */;
    double *Lcmp = (nc > 0 ? rmxn_alloc(nc, nf) : NULL); /* Lagrange multipliers. */
    
    int32_t rank; /* Rank of least squares system. */
    if (verbose) { fprintf(stderr, "  calling {lsq_solve_system}...\n\n"); }
    rank = lsq_solve_system(nx, nf, A, B, nc, R, S, Ucmp, Lcmp, verbose);
    demand(rank == nx + nc, "could not solve the least squares system");
    if (verbose) 
      { gsel_print_array(stderr, "%12.6f", "  fitted map matrix:", nx, nf, Ucmp, "\n");
        if (nc > 0) { gsel_print_array(stderr, "%12.6f", "  Lagrange multipliers:", nc, nf, Lcmp, "\n"); }
      }
    test_lsq_check_solution(nx, nf, tol, Uexp, Ucmp);

    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void test_lsq_throw_linear_system
  ( int32_t nx, 
    int32_t nf,
    int32_t nc,
    bool_t scramble,
    double **AP,
    double **BP, 
    double **RP,
    double **SP,
    double **UexpP
  )
  {
    /* Create a canonical minimization problem {A * U = B} with 
      constraints {R * U = S}  and solution {U}: */
    double *A = rmxn_alloc(nx, nx);
    double *B = rmxn_alloc(nx, nf);
     double *U = rmxn_alloc(nx, nf);
    
    /* Fill {A} with the identity, and cols {jf} of {B} with all {jf+1}: */
    for (int32_t ix = 0; ix < nx; ix++)
      { for (int32_t jx = 0; jx < nx; jx++)
          { A[ix*nx + jx] = (ix == jx ? 1.0 : 0.0); }
        for (int32_t jf = 0; jf < nf; jf++) 
          { B[ix*nf + jf] = (double)(jf + 1); }
      }
      
    double *R = NULL;
    double *S = NULL;
    if (nc > 0)
      { R = rmxn_alloc(nc, nx);
        S = rmxn_alloc(nc, nf);
        /* Each constraint requires one variable {U[ic,jf]} to be {-jf-1}: */
        for (int32_t ic = 0; ic < nc; ic++)
          { for (int32_t jx = 0; jx < nx; jx++)
              { R[ic*nx + jx] = (ic == jx ? 1.0 : 0.0); }
            for (int32_t jf = 0; jf < nf; jf++)
              { S[ic*nf + jf] = - (double)(jf + 1); }
          }
      }
      
    /* The constrained optimum: */
    for (int32_t ix = 0; ix < nx; ix++)
      { for (int32_t jf = 0; jf < nf; jf++)
          { U[ix*nf + jf] = (ix < nc ? -1.0 : +1.0)*((double)(jf + 1)); }
      }
      
    /* Should check...: */ 
    
    if (scramble)
      { /* Apply linear transformation to the system: */
        /* Allocate temporary matrices: */
        double *Aw = rmxn_alloc(nx, nx);
        double *Bw = rmxn_alloc(nx, nf);
        double *Uw = rmxn_alloc(nx, nf);

        /* Now apply linear transformation to the canonical problem: */
        double *M = rmxn_alloc(nx, nx); rmxn_throw_ortho(nx, M); rmxn_perturb_unif(nx, nx, 0.1, M);
        double *Minv = rmxn_alloc(nx, nx);
        rmxn_inv(nx, M, Minv);

        double *Q = rmxn_alloc(nf, nf); rmxn_throw_ortho(nf, Q); rmxn_perturb_unif(nf, nf, 0.1, Q);

        rmxn_mul(nx, nx, nx, M, A, Aw); rmxn_mul(nx, nx, nx, Aw, M, A);
        rmxn_mul(nx, nx, nf, M, B, Bw); rmxn_mul(nx, nf, nf, Bw, Q, B);

        if (nc > 0)
          { double *N = rmxn_alloc(nc, nc); rmxn_throw_ortho(nc, N); rmxn_perturb_unif(nc, nc, 0.1, N);
            double *Rw = rmxn_alloc(nc, nx);
            double *Sw = rmxn_alloc(nc, nf);

            rmxn_mul(nc, nc, nx, N, R, Rw); rmxn_mul(nc, nx, nx, Rw, M, R);
            rmxn_mul(nc, nc, nf, N, S, Sw); rmxn_mul(nc, nf, nf, Sw, Q, S);
            free(Rw); free(Sw);  free(N);
          }

        rmxn_mul(nx, nx, nf, Minv, U, Uw); rmxn_mul(nx, nf, nf, Uw, Q, U);

        /* Free work matrices: */
        free(Aw); free(Bw); free(Uw); free(Q);
        free(M); free(Minv);

      }

    /* Return results: */
    (*AP) = A; (*BP) = B;  (*RP) = R; (*SP) = S;  (*UexpP) = U;
  }

void test_lsq_check_solution(int32_t nx, int32_t nf, double tol, double Uexp[], double Ucmp[])
  { for (int32_t ix = 0; ix < nx; ix++)
      { for (int32_t jf = 0; jf < nf; jf++) 
          { double Uexpij = Uexp[ix*nf + jf];
            double Ucmpij = Ucmp[ix*nf + jf];
            double s = Ucmpij - Uexpij;
            if (fabs(s) > tol)
              { fprintf
                  ( stderr,
                    "(Ucmp - Uexp)[%d,%d] = %24.16e  tol = %24.16e\n", 
                    ix, jf, s, tol
                  );
                demand(FALSE, "** fit is not very good");
              }
          }
      }
    fprintf(stderr, "fit is good with tolerance %23.16e\n", tol);
  }
  
