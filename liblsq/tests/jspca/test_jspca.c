/* test_jspca --- test program for {jspca.h}  */
/* Last edited on 2023-02-02 22:41:43 by stolfi */

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
#include <rmxn_extra.h>

#include <jspca.h>

/* GENERAL PARAMETERS */

#define MAX_POINTS 10000
  /* Number of cases to generate. */

#define MAX_RUNS 100
  /* Max number of trials per test. */

#define MAX_VARS 10
  /* Max number of independent variables (argument coordinates per data point). */

/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);

void tpca_do_one_test(int32_t trial, bool_t verbose);

void tpca_throw_components(int32_t nv, double E[], double maxMag, double e[]);
  /* Generates a random {nv × nv} orthonormal matrix {E} and its presumed 
    square-rooted eigenvalues {e[0..nv-1]} in strictly decreasing magnitude order
    starting with {maxMag}. */

void tpca_throw_data_points(int32_t nd, int32_t nv, double E[], double e[], double d[], double D[]);
  /* Generates a set {D} of {nd} data points with {nv} components each, as a matrix {D} 
    with {nd} rows and {nv} columns.  Each row {D[id]} is the barycenter {d[0..nv-1]}
    plus {\SUM{ ie in 0..nv-1 : rnd()*e[ie]*E[ie] }} where {rnd()} is a normal random variable
    and {E[ie]} is row {ie} of {E}. */

/* CHECKING PCAS AND DECOMPOSITION */

void tpca_check_pcas
  ( int32_t nv, double minMag, 
    int32_t ne_exp, double E_exp[], double e_exp[], 
    int32_t ne_cmp, double E_cmp[], double e_cmp[]
  );
  /* Checks whether the computed principal components and magnitudes {E_cmp,e_cmp} 
    match the expected values {E_exp,e_exp}. */

void tpca_check_decomp(int32_t nd, int32_t nv, double D[], double d[], int32_t ne, double E[], double C[], double P[], double R[]);
  /* Checks whether the decomposition results {C,P,R} are consistent with the data {D,d,E}. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { int32_t i;
    for (i = 0; i < MAX_RUNS; i++) { tpca_do_one_test(i, i < 5); }
    return 0;
  }

void tpca_do_one_test(int32_t trial, bool_t verbose)
  { 
    srand(1665 + 2*trial);
    srandom(1665 + 2*trial);
    
    int32_t nv = int32_abrandom(1, MAX_VARS);            /* Number of coords of data points. */
    int32_t nd = int32_abrandom(2*nv + 10, MAX_POINTS);  /* Number of data points. */
    int32_t ne_exp = int32_abrandom(0, nv);              /* Number of significant pcas. */

    if (verbose) { fprintf(stderr, "choosing barycenter {d}...\n\n"); }
    double d[nv];
    rn_throw_ball(nv, d);
    rn_scale(nv, 9.0, d, d);

    if (verbose) { fprintf(stderr, "choosing the principal directions {E} and deviations {e}...\n\n"); }
    double *E_exp = rmxn_alloc(nv, nv); /* Chosen PCA matrix. */
    double e_exp[nv]; 
    double maxMag = 2.0;
    tpca_throw_components(nv, E_exp, maxMag, e_exp);
    assert(e_exp[nv-1] > 0);
    
    double minMag = (ne_exp > 0 ? e_exp[ne_exp-1] : maxMag + 1); 
    if (ne_exp < nv) { assert(e_exp[ne_exp] < minMag); }
    if (verbose) { fprintf(stderr, "minimum PC magnitude {minMag} = %12.8f...\n\n", minMag); }
    
    if (verbose) { fprintf(stderr, "generating the data array {D} and weight vector {w}...\n\n"); }
    double *D = rmxn_alloc(nd,nv);
    double *w = rn_alloc(nd);
    for (int32_t id = 0; id < nd; id++)
      { double *Di = &(D[id*nv]);
        for (int32_t jv = 0; jv < nv; jv++) { Di[jv] = d[jv]; }
        for (int32_t ke = 0; ke < nv; ke++)
          { double *Ek = &(E_exp[ke*nv]);
            double rk = dgaussrand()*e_exp[ke];
            for (int32_t jv = 0; jv < nv; jv++) { Di[jv] *= rk*Ek[jv]; }
          }
        w[id] = drandom();
      }
    
    /* Compute the PCAs from the data: */
    double *E_cmp = rmxn_alloc(nv, nv); /* Chosen PCA matrix. */
    double e_cmp[nv]; 
    int32_t ne_cmp = jspca_compute_components(nd, nv, D, w, d, minMag, E_cmp, e_cmp, verbose);
    
    /* Compare the results: */
    tpca_check_pcas(nv, minMag, ne_exp, E_exp, e_exp, ne_cmp, E_cmp, e_cmp);
    
    /* Try the PCA decomposition: */
    double *C = rmxn_alloc(nd, ne_cmp); /* PCA coeffs. */
    double *P = rmxn_alloc(nd, nv);     /* Projection of {D - u*d} on {E} row subspace. */
    double *R = rmxn_alloc(nv, nv);      /* Residual {D - u*d - P}. */
    jspca_decompose_data(nd, nv, D, d, ne_cmp, E_cmp, C, P, R, verbose); 
    
    /* Check the results: */
    tpca_check_decomp(nd, nv, D, d, ne_cmp, E_cmp, C, P, R);

    free(E_exp);
    free(D);
    free(w);
    free(E_cmp);
    free(C);
    free(P);
    free(R);

    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }

void tpca_throw_components(int32_t nv, double E[], double maxMag, double e[])
  { rmxn_throw_ortho(nv, E);
    double mag = maxMag;
    double att = (nv == 1 ? 1.0 : exp(log(0.1)/(nv-1)));
    for (int32_t ke = 0; ke < nv; ke++) 
      { e[ke] = mag;
        mag = att*mag;
      }
  }


void tpca_check_pcas
  ( int32_t nv, double minMag, 
    int32_t ne_exp, double E_exp[], double e_exp[], 
    int32_t ne_cmp, double E_cmp[], double e_cmp[]
  )
  { 
    fprintf(stderr, "!!! %s NOT IMPLEMENTED !!!\n", __FUNCTION__);
  }

void tpca_check_decomp(int32_t nd, int32_t nv, double D[], double d[], int32_t ne, double E[], double C[], double P[], double R[])
  { 
    fprintf(stderr, "!!! %s NOT IMPLEMENTED !!!\n", __FUNCTION__);
  }

