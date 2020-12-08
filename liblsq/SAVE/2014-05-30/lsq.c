/* See {lsq.h} */
/* Last edited on 2014-05-25 21:08:48 by stolfilocal */

#define lsq_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <rmxn.h>
#include <jsmath.h>
#include <gauss_elim.h>

#include <lsq.h>

int lsq_fit
  ( int nt,     /* Number of cases to generate. */
    int nx,     /* Number of independent variables. */
    int nf,     /* Number of dependent variables (functions to fit). */
    lsq_gen_case_t *gen_case,
    double U[], /* Fitted linear transformation matrix. */
    bool_t verbose
  )
  { 
    /* Lsq fitting system {A U = B} for {f[0..nf-1]} in terms of {x[0..nx-1]}: */
    double *A = rmxn_alloc(nx,nx);
    double *B = rmxn_alloc(nx,nf);
    lsq_compute_matrix_and_rhs(nt, nx, nf, gen_case, A, B, verbose);
      
    int rank = lsq_solve_system(nx, nf, A, B, U, verbose);
    free(B);
    free(A);
    return rank;
  }
      
void lsq_compute_matrix_and_rhs(int nt, int nx, int nf, lsq_gen_case_t *gen_case, double A[], double B[], bool_t verbose)
  {
    /* Allocate storage for data records cases: */
    double xi[nx];  /* Independent variables. */
    double fi[nf];  /* Dependent variables. */

    /* Generate all test cases, accumulate statistics: */
    rmxn_zero(nx, nx, A);
    rmxn_zero(nx, nf, B);
    int it;
    for (it = 0; it < nt; it++)
      { 
        bool_t verbacc = verbose & (it < 20); /* Debug the stats accumulator? */
        
        /* Obtain datum number {it} in {xi,fi,wi}: */
        double wi = NAN;
        gen_case(it, nx, xi, nf, fi, &wi);
        if (verbacc) 
          { fprintf(stderr, "  i = %5d", it);
            fprintf(stderr, "  xi =");
            lsq_debug_double_vec(nx, xi, "%6.3f");
            fprintf(stderr, "  fi =");
            lsq_debug_double_vec(nf, fi, "%6.3f");
            fprintf(stderr, "\n");
          }
        demand(wi >= 0, "case weight must be non-negative");
        demand(wi < +INF, "infinte case weights not implemented yet");
        
        /* Accumulate scalar products on matrix: */
        int iv;
        for (iv = 0; iv < nx; iv++)
          { int jv, jf;
            for (jv = 0; jv < nx; jv++)
              { A[iv*nx+jv] += xi[iv]*wi*xi[jv]; }
            for (jf = 0; jf < nf; jf++)
              { B[iv*nf+jf] += xi[iv]*wi*fi[jf]; }
          }
          
        if (verbacc) { fprintf(stderr, "\n"); }
      }
  }
  
int lsq_solve_system(int nx, int nf, double A[], double B[], double U[], bool_t verbose) 
  {
    if (verbose)
      { /* Print the least squares system: */
        fprintf(stderr, "  least squares systems:\n");
        int iv;
        for (iv = 0; iv < nx; iv++)
          { fprintf(stderr, "  %-4s", (iv == nx/2 ? "A = " : ""));
            lsq_debug_double_vec(nx, &(A[iv*nx]), "%12.5f");
            fprintf(stderr, "  %-4s", (iv == nx/2 ? "B = " : ""));
            lsq_debug_double_vec(nf, &(B[iv*nf]), "%12.5f");
            fprintf(stderr, "\n");
          }
        fprintf(stderr, "\n");
      }
          
    /* Solve the least squares system: */
    int rank = gsel_solve(nx, nx, A, nf, B, U, 0.0);
    if (verbose)
      { fprintf(stderr, "  rank = %d", rank);
        if (rank < nx) { fprintf(stderr, " (%s)", "indeterminate"); }
        fprintf(stderr, "\n");
      }
    return rank;
  }

void lsq_debug_double_vec(int nx, double x[], char *fmt)
  { fprintf(stderr, "[");
    int j;
    for (j = 0; j < nx; j++)
      { fprintf(stderr, " "); fprintf(stderr, fmt, x[j]); }
    fprintf(stderr, " ]");
  }
    
void lsq_debug_int_vec(int nx, int x[], char *fmt)
  { fprintf(stderr, "[");
    int j;
    for (j = 0; j < nx; j++)
      { fprintf(stderr, " "); fprintf(stderr, fmt, x[j]); }
    fprintf(stderr, " ]");
  }
