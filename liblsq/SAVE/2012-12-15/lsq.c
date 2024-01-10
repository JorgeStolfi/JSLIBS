/* See {lsq.h} */
/* Last edited on 2012-12-15 09:51:21 by stolfilocal */

#define lsq_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <lsq.h>

#include <bool.h>
#include <affirm.h>
#include <gauss_elim.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void lsq_debug_double_vec(double *v, int nv, char *fmt);
void lsq_debug_int_vec(int *v, int nv, char *fmt);
  /* These procedures print {v[0..nv-1]} to {stderr}, each with format {fmt},
    separated by spaces and bracketed by '[' and ']'. */

int lsq_fit
  ( int nv,     /* Number of independent variables. */
    int nf,     /* Number of dependent variables (functions to fit). */
    int nt,     /* Number of cases to generate. */
    lsq_gen_case_t *gen_case,
    double U[], /* Fitted linear transformation matrix. */
    bool_t verbose
  )
  { 
    /* Allocate storage for sample cases: */
    double v[nv];  /* Independent variables. */
    double f[nf];  /* Dependent variables. */

    /* Lsq fitting system {A U = B} for {f[0..nf-1]} in terms of {v[0..nv-1]}: */
    double A[nv*nv];  /* Has {nv} rows and {nv} columns. */
    double B[nv*nf];  /* Has {nv} rows and {nf} columns. */
    int k;
    for (k = 0; k < nv*nv; k++) { A[k] = 0.0; }
    for (k = 0; k < nv*nf; k++) { B[k] = 0.0; }
    
    /* Generate all test cases, accumulate statistics: */
    int it;
    for (it = 0; it < nt; it++)
      { 
        bool_t verbacc = verbose & (it < 20); /* Debug the stats accumulator? */
        
        /* Compute {X,Y,Z} for a random seq pair {x,y}: */
        gen_case(it, nv, v, nf, f);
        if (verbacc) 
          { fprintf(stderr, "  v =");
            lsq_debug_double_vec(v, nv, "%6.3f");
            fprintf(stderr, "  f =");
            lsq_debug_double_vec(f, nf, "%6.3f");
            fprintf(stderr, "\n");
          }
        
        /* Accumulate scalar products on matrix: */
        int iv;
        for (iv = 0; iv < nv; iv++)
          { int jv, jf;
            for (jv = 0; jv < nv; jv++)
              { A[iv*nv+jv] += v[iv]*v[jv]; }
            for (jf = 0; jf < nf; jf++)
              { B[iv*nf+jf] += v[iv]*f[jf]; }
          }
          
        if (verbacc) { fprintf(stderr, "\n"); }
      }
      
    if (verbose)
      { /* Print the least squares system: */
        fprintf(stderr, "  least squares systems:\n");
        int iv;
        for (iv = 0; iv < nv; iv++)
          { fprintf(stderr, "  %-4s", (iv == nv/2 ? "A = " : ""));
            lsq_debug_double_vec(&(A[iv*nv]), nv, "%12.5f");
            fprintf(stderr, "  %-4s", (iv == nv/2 ? "B = " : ""));
            lsq_debug_double_vec(&(B[iv*nf]), nf, "%12.5f");
            fprintf(stderr, "\n");
          }
        fprintf(stderr, "\n");
      }
          
    /* Solve the least squares system: */
    int rank = gsel_solve(nv, nv, A, nf, B, U, 0.0);
    return rank;
  }

void lsq_debug_double_vec(double *v, int nv, char *fmt)
  { fprintf(stderr, "[");
    int j;
    for (j = 0; j < nv; j++)
      { fprintf(stderr, " "); fprintf(stderr, fmt, v[j]); }
    fprintf(stderr, " ]");
  }
    
void lsq_debug_int_vec(int *v, int nv, char *fmt)
  { fprintf(stderr, "[");
    int j;
    for (j = 0; j < nv; j++)
      { fprintf(stderr, " "); fprintf(stderr, fmt, v[j]); }
    fprintf(stderr, " ]");
  }
