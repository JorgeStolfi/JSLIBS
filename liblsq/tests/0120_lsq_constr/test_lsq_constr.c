/* test_lsq --- test program for constrained {lsq_solve_system}  */
/* Last edited on 2024-11-08 16:42:48 by stolfi */

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
   
#define MAT_ABS_PERT 1.0e-12
#define MAT_REL_PERT 1.0e-6
  /* Absolute and relative perturbation on orthogonal matrix elems
    for scrambling variables and equations. */

/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);

void tlsq_build_and_solve_system
  ( int32_t trial,
    int32_t nx,
    int32_t nf,
    int32_t nc,
    double tol,
    bool_t scrambleVars,
    bool_t scrrambleEqs,
    bool_t verbose
  );
  /* Tests the least-squares fitter with tolerance {tol}, {nx}
     variables, {nf} functions, {nc} constraints. The number {nt} of
     data points is chosen internally.
     
     First builds a trivial problem with a simple solution.
     Then modifies it as requested by {scrambleVars} and/or {scrambleEqs}. */

void tlsq_make_trivial_problem
  ( int32_t nx, 
    int32_t nf,
    int32_t nc,
    double **AP,
    double **BP, 
    double **RP,
    double **SP,
    double **UP,
    double **LP,
    bool_t verbose
  );
  /* Creates a system of equations {A*U + R'*L = B}, {R*U = S}, suitable for {lsq_solve_system},
    for a trivial constrained least squares problem. The matrices are allocated by the
    procedure and returned in {*AP}, {*BP}, etc.
    
    The system will have {nx} unknowns (coefficients to fit), {nf} functions to
    fit, and {nc} constraints.  Namely,
    
      (A} {nx×nx} is the moment matrix;
      {B} {nx×nf} is the right-hand-side matrix;
      {R} {nc×nx} is the constraint coefficient matrix;
      {S} {nc×nf} is the right-hand-side of the constraints;
      {U} {nx×nf} is the expected solution. 
      {L} {nc×nf} is the Lagrangian matrix of the constraints.
      
    The combined matrix {[[A,R'],[R,0]]} will be a permutation of the identity,
    and the solution {U} will be the numbers from 1 to {nx*nf}. */

void tlsq_scramble_problem_variables
  ( int32_t nx, 
    int32_t nf,
    int32_t nc,
    double A[],
    double B[], 
    double R[],
    double S[],
    double U[],
    double L[],
    bool_t verbose
  );
  /* Modifies the matrices {A,B,R,S,U,L} of a constrained least squares system problem 
    equivalent to applying a random non-singular linear mapping to the variables. */

void tlsq_scramble_problem_equations
  ( int32_t nx, 
    int32_t nf,
    int32_t nc,
    double A[],
    double B[], 
    double R[],
    double S[],
    double U[],
    double L[],
    bool_t verbose
  );
  /* Modifies the matrices {A,B,R,S,U,L} of a constrained least squares system problem 
    by applying random non-singular linear mappings to the equations. */
  
/* CHECKING LSQ FITTING */

void tlsq_compare_with_expected_soln(int32_t nx, int32_t nf, double tol, double Uexp[], double Ucmp[]);
  /* Checks whether the computed solution matrix {Ucmp[]} matches the expected 
    solution matrix {Uexp}, within the tolerance {tol}. */

void tlsq_check_main_system(int32_t nx, int32_t nc, int32_t nf, double tol, double A[], double U[], double R[], double L[], double B[], bool_t verbose);
  /* Computes the residual {A*U + R'*L - B} of the main system and checks if the max abs elem does not exceed {tol}. */
  
void tlsq_check_constraints(int32_t nx, int32_t nc, int32_t nf, double tol, double R[], double S[], double U[], bool_t verbose);
  /* Computes the residual {R*U - S} of the constraint equations,
    prints max elem. */

void tlsq_check_system
  ( int32_t nx,
    int32_t nc,
    int32_t nf,
    double tol,
    double A[],
    double B[],
    double R[],
    double S[],
    double U[],
    double L[],
    bool_t verbose
  );
  /* Calls {tlsq_check_main_system} and {tlsq_check_constraints} with the given arguments. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { int32_t trial = 0;
    tlsq_build_and_solve_system(trial, 1, 1, 0, 1.0e-9, FALSE, FALSE,  TRUE); trial++; /* Try 1x1 trivial, no constraints. */
    tlsq_build_and_solve_system(trial, 2, 3, 0, 1.0e-9, FALSE, FALSE,  TRUE); trial++; /* Try 2x3 trivial, no constraints. */
    tlsq_build_and_solve_system(trial, 2, 1, 1, 1.0e-8, FALSE, FALSE,  TRUE); trial++; /* Try 2x1 trivial, 1 constraint. */ 
    tlsq_build_and_solve_system(trial, 2, 1, 1, 1.0e-8, TRUE,  FALSE,  TRUE); trial++; /* Try 2x1 scrambled vars, 1 constraint. */
    tlsq_build_and_solve_system(trial, 2, 1, 1, 1.0e-8, FALSE, TRUE,   TRUE); trial++; /* Try 2x1 scrambled vars, 1 constraint. */
    while (trial < MAX_RUNS) 
      { srand(1665 + 2*trial);
        srandom(1665 + 2*trial);
        int32_t max_nx = 2*trial;
        int32_t max_nf = trial+2;
        int32_t nx = int32_abrandom(max_nx/2, max_nx); /* Number of indep variables. */
        int32_t nf = int32_abrandom(max_nf/2, max_nf); /* Number of dep variables. */
        int32_t nc = int32_abrandom(nx/2, nx-1);       /* Number of constraints. */
        double tol = 1.0e-8*nx;
        bool_t scrambleVars = TRUE;
        bool_t scrambleEqs = TRUE;
        bool_t verbose = (nx+nc+nf <= 12);
        if ((nx  == 14) && (nf == 9) && (nc == 8)) { verbose = TRUE; }
        tlsq_build_and_solve_system(trial+2, nx, nf, nc, tol, scrambleVars, scrambleEqs, verbose);
        trial++;
      }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void tlsq_build_and_solve_system
  ( int32_t trial,
    int32_t nx,
    int32_t nf,
    int32_t nc,
    double tol,
    bool_t scrambleVars,
    bool_t scrambleEqs,
    bool_t verbose
  )
  { 
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s (%d)\n", __FUNCTION__, trial);
    fprintf(stderr, "testing with nx = %d  nf = %d  nc = %d tol = %23.16e", nx, nf, nc, tol);
    fprintf(stderr, " scrambleVars = %c", "FT"[scrambleVars]);
    fprintf(stderr, " scrambleEqs = %c", "FT"[scrambleEqs]);
    fprintf(stderr, " ...\n");
    
    /* Generate a system with known solution {Uexp}: */
    double *A, *B, *R, *S, *Uexp, *Lexp;
    tlsq_make_trivial_problem(nx, nf, nc, &A, &B, &R, &S, &Uexp, &Lexp, verbose);
    assert(A != NULL); assert(B != NULL); assert(Uexp != NULL);
    if (nc > 0) { assert(R != NULL); assert(S != NULL); assert(Lexp != NULL); }
    if (scrambleVars) { tlsq_scramble_problem_variables(nx, nf, nc, A, B, R, S, Uexp, Lexp, verbose); }
    if (scrambleEqs) { tlsq_scramble_problem_equations(nx, nf, nc, A, B, R, S, Uexp, Lexp, verbose); }
    
    /* Solve the system with {lsq_solve_system}: */
    double *Ucmp = rmxn_alloc(nx, nf); /* Computed solution. */;
    double *Lcmp = (nc > 0 ? rmxn_alloc(nc, nf) : NULL); /* Lagrange multipliers. */
    
    int32_t rank; /* Rank of least squares system. */
    if (verbose) { fprintf(stderr, "  calling {lsq_solve_system}...\n\n"); }
    rank = lsq_solve_system(nx, nf, A, B, nc, R, S, Ucmp, Lcmp, verbose);
    demand(rank == nx + nc, "could not solve the least squares system");
    tlsq_check_system(nx, nc, nf, tol, A, B, R, S, Ucmp, Lcmp, verbose);
    
    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void tlsq_make_trivial_problem
  ( int32_t nx, 
    int32_t nf,
    int32_t nc,
    double **AP,
    double **BP, 
    double **RP,
    double **SP,
    double **UP,
    double **LP,
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "> --- %s -----------------------------------------------------\n", __FUNCTION__); }
    
    demand(nx >= 1, "invalid num of variables");
    demand(nc <= nx, "too many constraints");
    demand(nf >= 1, "invalid num of fitting problems");
    
    /* Create a canonical minimization problem {A * U = B} with 
      constraints {R * U = S}, solution {U}, Lagrangian {L}: */
    double *A = rmxn_alloc(nx, nx);
    double *B = rmxn_alloc(nx, nf);
    double *U = rmxn_alloc(nx, nf);

    double *R = NULL;
    double *S = NULL;
    double *L = NULL;
      
    /* Define {A,B,U} with a trivial system and distinctive solution: */
    for (uint32_t ix = 0;  ix < nx; ix++)
      { for (uint32_t jx = 0;  jx < nx; jx++)
          { /* Fill {A} with identity except first {nc} rows and columns: */
            A[ix*nx + jx] = ((ix >= nc ) && (ix == jx) ? 1.0 : 0.0);
          }
        for (uint32_t jf = 0;  jf < nf; jf++)
          { /* Each element of {U} is a distinct nonzero integer: */
            U[ix*nf + jf] = ((double)(jf*nx + ix + 1));
            /* The first {nc} rows of {B} are more distinct integers, the rest are the same as {U}: */
            B[ix*nf + jf] = (ix < nc ? nx*nf + jf*nc + ix : U[ix*nf + jf]);
          }
      }
      
    /* Define {R,L,S} : */
    if (nc > 0)
      { R = rmxn_alloc(nc, nx);
        L = rmxn_alloc(nc, nf);
        S = rmxn_alloc(nc, nf);
        for (uint32_t ic = 0;  ic < nc; ic++)
          { for (uint32_t jx = 0;  jx < nx; jx++)
              { /* The first {nc} columns of {R} are the identity, rest zeros: */
                R[ic*nx + jx] = (ic == jx ? 1.0 : 0.0);
              }
            for (uint32_t jf = 0;  jf < nf; jf++)
              { /* The constraint RHS {L} is the same as the first {nc} rows of {U}: */
                S[ic*nf + jf] = U[ic*nf + jf];
                /* The Lagrangian {L} is the first {nc} rows of {B}: */
                L[ic*nf + jf] = B[ic*nf + jf];
              }
          }
      }
      
    if (verbose) { lsq_print_problem(stderr, 4, "%12.6f", "trivial problem:\n", nx, nc, nf, A, B, R, S, "Uexp", U, "Lexp", L); }
    tlsq_check_system(nx, nc, nf, 1.0e-16, A, B, R, S, U, L, verbose);
    if (verbose) { fprintf(stderr, "\n"); }

    /* Return results: */
    (*AP) = A; (*BP) = B;  (*RP) = R; (*SP) = S;  (*UP) = U; (*LP) = L;

    if (verbose) { fprintf(stderr, "< --- %s -----------------------------------------------------\n", __FUNCTION__); }
  }

void tlsq_scramble_problem_variables
  ( int32_t nx, 
    int32_t nf,
    int32_t nc,
    double A[],
    double B[], 
    double R[],
    double S[],
    double U[],
    double L[],
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "> --- %s -----------------------------------------------------\n", __FUNCTION__); }
    demand(nx >= 1, "invalid num of variables");
    demand(nc <= nx, "too many constraints");
    demand(nf >= 1, "invalid num of fitting problems");

    { /* Apply linear transformations to the variables and equations: */

      /* Apply linear transformation {Q} to the variables: */
      double *Q = rmxn_alloc(nx, nx); 
      rmxn_throw_ortho(nx, Q); 
      rmxn_perturb_unif(nx, nx, MAT_ABS_PERT, MAT_REL_PERT, Q);
      double *Qinv = rmxn_alloc(nx, nx);
      rmxn_inv(nx, Q, Qinv);

      /* Temporary matrices: */
      double *Aw = rmxn_alloc(nx, nx);
      double *Uw = rmxn_alloc(nx, nf);

      rmxn_mul(nx, nx, nx, A, Qinv, Aw);
      rmxn_copy(nx, nx, Aw, A);

      rmxn_mul(nx, nx, nf, Q, U, Uw);
      rmxn_copy(nx, nf, Uw, U);
      
      if (nc > 0)
        { double *Rw = rmxn_alloc(nc, nx);
          double *Rd = rmxn_alloc(nc, nx);
          rmxn_mul(nc, nx, nx, R, Qinv, Rw);
          rmxn_sub(nc, nx, Rw, R, Rd);
          rmxn_copy(nc, nx, Rw, R);
          
          /* Recompute {B} to preserve the Lagrangian: */
          double *Bd = rmxn_alloc(nx, nf);
          rmxn_tr_mul(nc, nx, nf, Rd, L, Bd);
          rmxn_add(nx, nf, Bd, B, B);
          
          free(Rw);
          free(Rd);
          free(Bd);
        }

      free(Q);
      free(Qinv);
      free(Aw);
      free(Uw);
    }

    if (verbose) { lsq_print_problem(stderr, 4, "%12.6f", "scrambled variables problem:\n", nx, nc, nf, A, B, R, S, "Uexp", U, "Lexp", L); }
    tlsq_check_system(nx, nc, nf, nx*1.0e-12, A, B, R, S, U, L, verbose);
    if (verbose) { fprintf(stderr, "\n"); }

    if (verbose) { fprintf(stderr, "< --- %s -----------------------------------------------------\n", __FUNCTION__); }
  }

void tlsq_scramble_problem_equations
  ( int32_t nx, 
    int32_t nf,
    int32_t nc,
    double A[],
    double B[], 
    double R[],
    double S[],
    double U[],
    double L[],
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "> --- %s -----------------------------------------------------\n", __FUNCTION__); }
    demand(nx >= 1, "invalid num of variables");
    demand(nc <= nx, "too many constraints");
    demand(nf >= 1, "invalid num of fitting problems");

    { /* Apply linear transformations to the variables and equations: */

      /* Temporary matrices: */
      double *Aw = rmxn_alloc(nx, nx);
      double *Bw = rmxn_alloc(nx, nf);
      double *Uw = rmxn_alloc(nx, nf);

      /* Apply linear transformation {M} to the main equations: */
      double *M = rmxn_alloc(nx, nx); 
      rmxn_throw_ortho(nx, M); 
      rmxn_perturb_unif(nx, nx, MAT_ABS_PERT, MAT_REL_PERT, M);
      double *Minv = rmxn_alloc(nx, nx);
      rmxn_inv(nx, M, Minv);

      rmxn_mul(nx, nx, nx, M, A, Aw);
      rmxn_copy(nx, nx, Aw, A);

      rmxn_mul(nx, nx, nf, M, B, Bw);
      rmxn_copy(nx, nf, Bw, B);

      if (nc > 0)
        { double *Rw = rmxn_alloc(nc, nx);
          rmxn_mul_tr(nc, nx, nx, R, M, Rw);
          rmxn_copy(nc, nx, Rw, R);
          free(Rw);
        }

      free(M);
      free(Minv);
      free(Aw);
      free(Bw);
      free(Uw);
    }

  if (nc > 0)
    { /* Apply linear transformation {N} to the constraint equations: */
      double *N = rmxn_alloc(nc, nc); 
      rmxn_throw_ortho(nc, N); 
      rmxn_perturb_unif(nc, nc, MAT_ABS_PERT, MAT_REL_PERT, N);
      double *Ninv = rmxn_alloc(nc, nc);
      rmxn_inv(nc, N, Ninv);
      
      double *Rw = rmxn_alloc(nc, nx);
      double *Rd = rmxn_alloc(nc, nx);
      
      rmxn_mul(nc, nc, nx, N, R, Rw);
      rmxn_sub(nc, nx, Rw, R, Rd);
      rmxn_copy(nc, nx, Rw, R);
      
      double *Bd = rmxn_alloc(nx, nf);
      rmxn_tr_mul(nc, nx, nf, Rd, L, Bd);
      rmxn_add(nx, nf, Bd, B, B);
    
      /* Recompute the RHS {S} of the constraints: */
      rmxn_mul(nc, nx, nf, R, U, S);
    }

    if (verbose) { lsq_print_problem(stderr, 4, "%12.6f", "scrambled equations problem:\n", nx, nc, nf, A, B, R, S, "Uexp", U, "Lexp", L); }
    tlsq_check_system(nx, nc, nf, nx*1.0e-12, A, B, R, S, U, L, verbose);
    if (verbose) { fprintf(stderr, "\n"); }

    if (verbose) { fprintf(stderr, "< --- %s -----------------------------------------------------\n", __FUNCTION__); }
  }

void tlsq_check_system
  ( int32_t nx,
    int32_t nc,
    int32_t nf,
    double tol,
    double A[],
    double B[],
    double R[],
    double S[],
    double U[],
    double L[],
    bool_t verbose
  )
  {
    tlsq_check_main_system(nx, nc, nf, tol, A, U, R, L, B, verbose);
    if (nc > 0) { tlsq_check_constraints(nx, nc, nf, tol, R, S, U, verbose); }
  }

void tlsq_compare_with_expected_soln(int32_t nx, int32_t nf, double tol, double Uexp[], double Ucmp[])
  { for (uint32_t ix = 0;  ix < nx; ix++)
      { for (uint32_t jf = 0;  jf < nf; jf++) 
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
    fprintf(stderr, "  matches expected solution good with tolerance %23.16e\n", tol);
  }

void tlsq_check_main_system(int32_t nx, int32_t nc, int32_t nf, double tol, double A[], double U[], double R[], double L[], double B[], bool_t verbose)
  { 
    double *Y = rmxn_alloc(nx, nf);
    rmxn_mul(nx, nx, nf, A, U, Y);
    double *Z = rmxn_alloc(nx, nf);
    rmxn_tr_mul(nc, nx, nf, R, L, Z);
    double *BE = rmxn_alloc(nx, nf);
    rmxn_add(nx, nf, Y, Z, BE);
    rmxn_sub(nx, nf, B, BE, BE);
    if (verbose)
      { gauss_elim_print_array(stderr, 4, "%12.6f", "residual of main system:", nx, nf, "BE", BE, ""); }
    double maxE = rmxn_max_abs_elem(nx, nf, BE);
    fprintf(stderr, "  max main system residual %23.16e\n", maxE);
    demand(maxE <= tol, "main system residual is too large");
  }

void tlsq_check_constraints(int32_t nx, int32_t nc, int32_t nf, double tol, double R[], double S[], double U[], bool_t verbose)
  { double *Y = rmxn_alloc(nc, nf);
    rmxn_mul(nc, nx, nf, R, U, Y);
    double *SE = rmxn_alloc(nc, nf);
    rmxn_sub(nc, nf, Y, S, SE);
    if (verbose)
      { gauss_elim_print_array(stderr, 4, "%12.6f", "residual of constraints:", nc, nf, "SE", SE, ""); }
    double maxE = rmxn_max_abs_elem(nc, nf, SE);
    fprintf(stderr, "  max constraint residual %23.16e\n", maxE);
    demand(maxE <= tol, "constraint residual is too large");
  }
