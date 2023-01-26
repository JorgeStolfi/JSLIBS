/* See {lsq.h} */
/* Last edited on 2023-01-24 10:12:24 by stolfi */

#define lsq_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <rmxn.h>
#include <jsmath.h>
#include <gauss_elim.h>

#include <lsq.h>

/* INTERNAL PROTOTYPES */

void lsq_assemble_constrained_system
  ( int32_t nx, 
    int32_t nf, 
    double A[], 
    double B[], 
    int32_t nc, 
    double R[], 
    double S[], 
    double AR[], 
    double BS[]
  );
  /* Builds the linear system of constrained quadratic minimization, 
    given the matrices {A} ({nx} by {nx}) and  {B} ({nx} by {nf})
    of the unconstrained minimization system {A * U = B}, and the matrices
    {R} ({nc} by {nx}) and {S} ({nc} by {nf}) of the constraint
    equations {R * U = S}.  
    
    Stores the augmented matries into {AR} ({nxc} by {nxc}, where {nxc = nx+nc})
    and {BS} ({nxc} by {nf}). */

void lsq_split_constrained_solution
  ( int32_t nx, 
    int32_t nf, 
    int32_t nc, 
    double UL[], 
    double U[], 
    double L[]
  );
  /* Splits the augmented solution matrix {UL} ({nxc} by {nf}, where {nxc = nx+nc})
    into the solution matrix proper {U} ({nx} by {nf}) and the Lagrange
    multiplier matrix {L} ({nc} by {nf}).  
    
    The parameter {L} may be {NULL}, in which case the Lagrange 
    multipliers are ignored. */

/* IMPLEMENTATIONS */

int32_t lsq_fit
  ( int32_t nt,     /* Number of data points to generate. */
    int32_t nx,     /* Number of independent variables (argument coordinates per data point). */
    int32_t nf,     /* Number of dependent variables (functions samples per data point). */
    lsq_gen_data_point_t *gen_data_point,
    double U[], /* Fitted linear transformation matrix. */
    bool_t verbose
  )
  { 
    /* Lsq fitting system {A U = B} for {f[0..nf-1]} in terms of {x[0..nx-1]}: */
    double *A = rmxn_alloc(nx,nx);
    double *B = rmxn_alloc(nx,nf);
    lsq_compute_matrix_and_rhs(nt, nx, nf, gen_data_point, A, B, verbose);
      
    int32_t rank = lsq_solve_system(nx, nf, A, B, 0,NULL,NULL, U,NULL, verbose);
    free(B);
    free(A);
    return rank;
  }
      
void lsq_compute_matrix_and_rhs
  ( int32_t nt, 
    int32_t nx, int32_t nf, 
    lsq_gen_data_point_t *gen_data_point, 
    double A[], 
    double B[], 
    bool_t verbose
  )
  {
    /* Allocate storage for data points: */
    double Xk[nx];  /* Independent variables (argument coordinates). */
    double Fk[nf];  /* Dependent variables (function samples). */

    /* Generate all test data points, accumulate statistics: */
    rmxn_zero(nx, nx, A);
    rmxn_zero(nx, nf, B);
    for (int32_t kt = 0; kt < nt; kt++)
      { 
        bool_t verbacc = verbose & (kt < 20); /* Debug the stats accumulator? */
        
        /* Obtain data point number {kt} in {Xk,Fk,Wk}: */
        double Wk = NAN;
        gen_data_point(kt, nx, Xk, nf, Fk, &Wk);
        if (verbacc) 
          { fprintf(stderr, "  k = %5d", kt);
            fprintf(stderr, "  Xk =");
            lsq_debug_double_vec(nx, Xk, "%9.6f");
            fprintf(stderr, "  Fk =");
            lsq_debug_double_vec(nf, Fk, "%9.6f");
            fprintf(stderr, "  Wk = %10.6f", Wk);
            fprintf(stderr, "\n");
          }
        demand(Wk >= 0, "reliability weight must be non-negative");
        demand(Wk < +INF, "infinte reliability weights not implemented yet");
        
        /* Accumulate scalar products on matrix: */
        for (int32_t ix = 0; ix < nx; ix++)
          { for (int32_t jx = 0; jx < nx; jx++)
              { A[ix*nx+jx] += Xk[ix]*Wk*Xk[jx]; }
            for (int32_t jf = 0; jf < nf; jf++)
              { B[ix*nf+jf] += Xk[ix]*Wk*Fk[jf]; }
          }
          
        if (verbacc) { fprintf(stderr, "\n"); }
      }
  }
  
int32_t lsq_solve_system
  ( int32_t nx, 
    int32_t nf, 
    double A[], 
    double B[], 
    int32_t nc, 
    double R[], 
    double S[], 
    double U[], 
    double L[], 
    bool_t verbose
  )
  {
    assert(nx >= 0);
    assert(nf >= 0);
    assert(nc >= 0);
    assert((A != NULL) && (B != NULL) && (U != NULL));
    if (nc > 0) { assert((R != NULL) && (S != NULL)); }
      
    if (verbose)
      { /* Print the least squares system: */
        fprintf(stderr, "  main least squares system:\n");
        for (int32_t ix = 0; ix < nx; ix++)
          { fprintf(stderr, "  %-4s", (ix == nx/2 ? "A = " : ""));
            lsq_debug_double_vec(nx, &(A[ix*nx]), "%12.5f");
            fprintf(stderr, "  %-4s", (ix == nx/2 ? "B = " : ""));
            lsq_debug_double_vec(nf, &(B[ix*nf]), "%12.5f");
            fprintf(stderr, "\n");
          }
        fprintf(stderr, "\n");
        if (nc > 0)
          { fprintf(stderr, "  constraints:\n");
            for (int32_t ix = 0; ix < nc; ix++)
              { fprintf(stderr, "  %-4s", (ix == nc/2 ? "R = " : ""));
                lsq_debug_double_vec(nx, &(R[ix*nx]), "%12.5f");
                fprintf(stderr, "  %-4s", (ix == nx/2 ? "S = " : ""));
                lsq_debug_double_vec(nf, &(S[ix*nf]), "%12.5f");
                fprintf(stderr, "\n");
              }
            fprintf(stderr, "\n");
          }
      }
          
    /* Solve the least squares system: */
    int32_t rank;
    if (nc == 0)
      { /* Just solve {A * U = B}: */
        rank = gsel_solve(nx, nx, A, nf, B, U, 0.0);
      }
    else
      { /* Build the augmented matrices {AR}, {BS}, {UL}: */
        int32_t nxc = nx + nc;
        double *AR = notnull(malloc(nxc*nxc*sizeof(double)), "no mem");
        double *BS = notnull(malloc(nxc*nf*sizeof(double)), "no mem");
        double *UL = notnull(malloc(nxc*nf*sizeof(double)), "no mem");
        lsq_assemble_constrained_system(nx, nf, A, B, nc, R, S, AR, BS);
        rank = gsel_solve(nxc, nxc, AR, nf, BS, UL, 0.0);
        lsq_split_constrained_solution(nx, nf, nc, UL, U, L);
        free(UL); free(BS); free(AR);
      } 
    if (verbose)
      { fprintf(stderr, "  rank = %d", rank);
        if (rank < nx) { fprintf(stderr, " (%s)", "indeterminate"); }
        fprintf(stderr, "  system solution:\n");
        for (int32_t ix = 0; ix < nx; ix++)
          { fprintf(stderr, "  %-4s", (ix == nx/2 ? "U = " : ""));
            lsq_debug_double_vec(nf, &(U[ix*nf]), "%12.5f");
            fprintf(stderr, "\n");
          }
        fprintf(stderr, "\n");
      }
    return rank;
  }

void lsq_assemble_constrained_system
  ( int32_t nx, 
    int32_t nf, 
    double A[], 
    double B[], 
    int32_t nc, 
    double R[], 
    double S[], 
    double AR[], 
    double BS[]
  )
  { int32_t nxc = nx + nc;
    for (int32_t ixc = 0; ixc < nxc; ixc++)
      { for (int32_t jxc = 0; jxc < nxc; jxc++)
          { double *ARij = &(AR[ixc*nxc +jxc]);
            if (ixc < nx)
              { if (jxc < nx) 
                  { /* Copy {A}: */  (*ARij) = A[ixc*nx + jxc]; }
                else
                  { /* Copy {R} transp.: */ (*ARij) = R[(jxc-nx)*nx + ixc]; }
              }
            else
              { if (jxc < nx) 
                  { /* Copy {R}: */  (*ARij) = R[(ixc-nx)*nx + jxc]; }
                else
                  { /* Zero: */ (*ARij) = 0.0; }
              }
          }
        for (int32_t jf = 0; jf < nf; jf++)
          { double *BSij = &(BS[ixc*nf + jf]);
            if (ixc < nx)
              { /* Copy {B}: */ (*BSij) = B[ixc*nf + jf]; }
            else
              { /* Copy {S}: */ (*BSij) = S[(ixc-nx)*nf + jf]; }
          }
      }
  }

void lsq_split_constrained_solution
  ( int32_t nx, 
    int32_t nf, 
    int32_t nc, 
    double UL[], 
    double U[], 
    double L[]
  )
  { int32_t nxc = nx + nc;
    for (int32_t ixc = 0; ixc < nxc; ixc++)
      { for (int32_t jf = 0; jf < nf; jf++)
          { double *ULij = &(UL[ixc*nf + jf]);
            if (ixc < nx)
              { /* Copy {U}: */ U[ixc*nf + jf] = (*ULij); }
            else if (L  != NULL)
              { /* Copy {L}: */  U[(ixc-nx)*nf + jf] = (*ULij); }
          }
      }
  }

void lsq_debug_double_vec(int32_t nx, double x[], char *fmt)
  { fprintf(stderr, "[");
    int32_t j;
    for (j = 0; j < nx; j++)
      { fprintf(stderr, " "); fprintf(stderr, fmt, x[j]); }
    fprintf(stderr, " ]");
  }
    
void lsq_debug_int32_vec(int32_t nx, int32_t x[], char *fmt)
  { fprintf(stderr, "[");
    int32_t j;
    for (j = 0; j < nx; j++)
      { fprintf(stderr, " "); fprintf(stderr, fmt, x[j]); }
    fprintf(stderr, " ]");
  }
