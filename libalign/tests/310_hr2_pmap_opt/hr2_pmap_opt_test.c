#define PROG_NAME "hr2_pmap_opt_test"
#define PROG_DESC "test of {hr2_pmap_opt.h} and related modules"
#define PROG_VERS "1.0"

/* Last edited on 2024-09-17 21:53:54 by stolfi */ 
/* Created on 2020-07-11 by J. Stolfi, UNICAMP */

#define hr2_pmap_opt_test_COPYRIGHT \
  "Copyright © 2020  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <ix.h>
#include <r2.h>
#include <hr2.h>
#include <rn.h>
#include <r3x3.h>
#include <bool.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>
#include <argparser.h>

#include <hr2_pmap_throw_type.h>
#include <hr2_pmap_from_many_pairs.h>
#include <hr2_pmap_encode.h>

#include <hr2_pmap_opt.h>

/* The program takes one command line option, the type of map to 
  test ("IDENTITY", "TRANSLATION", etc.). */

void test_hr2_pmap_opt_quadratic__hr2_pmap_opt_1D_plot(hr2_pmap_type_t type, bool_t verbose, char *outPrefix);
  /* Tests {hr2_pmap_opt_quadratic}, and, if {outPrefix} is not {NULL}, also
    {hr2_pmap_opt_1D_plot}, creating files with prefix {outPrefix}. */

int32_t main(int32_t argn, char **argv);

void hpot_do_test_opt_plot
  ( hr2_pmap_type_t type,
    sign_t sgn,
    bool_t tight,
    bool_t ident,
    bool_t verbose,
    char *outPrefix
  );
  /* Tests {hr2_pmap_opt_quadratic} and, if {outPrefix} is not {NULL},
    also {hr2_opt_1D_plot}.
  
    The goal function for a map {M} is the discrepancy squared between a certain number 
    {np} of random point pairs, mapped forward and backwards by {M}.
    The points are chosen so that the optimum map has the specified {sgn}.
    
    If {tight}, tests with the minimum number of points needed
    to determine the map of given {type}.
    
    If {ident}, the points are exactly related by the identity
    transformation (that is, each point is to be mapped to itself) if
    {sgn} is {+1} or by a Y coordinate negation if {sgn} is {-1}.
    
    Plots the goal function for the initial and final maps.  The plots
    will extend {urad} away from the map in question.*/

void hpot_do_test_1D_plot
  ( char *outPrefix,
    hr2_pmap_type_t type,
    sign_t sgn,
    bool_t tight,
    bool_t ident,
    hr2_pmap_t *M,
    hr2_pmap_opt_func_t *f2,
    double urad,
    char *stage,
    bool_t verbose
  );
  /* Plots the goal function {f2} in the neighborhood of the projective map {M},
    which must be of the given {type} and {sgn}.
    
    The plots will extend {urad} away from {M} in encoding parameter space.
    The file name will be "{outPrefix}-{type}-{xsgn}-tight{tight}-ident{ident}-{stage}.txt where
    {xsign} will be "m", "o", or "p" for {-1}, 0, and {+1} respectively.  */
  
void hpot_choose_r2_point_pairs
  ( hr2_pmap_type_t type,
    sign_t sgn,
    bool_t tight,
    bool_t ident,
    int32_t *np_P,
    r2_t **p1_P,
    r2_t **p2_P,
    double **w_P
  );
  /* Chooses a suitable number {np} of point pairs {p1[i],p2[i]}
    and weights {w[i]}, with {i} in {0..np-1}, for testing the quadratic optimization
    with given {type} and {sgn}.  
    
    If {tight}, {np} will be the minimum number {nr} of points needed to
    define such a map, and that map will be exact.
    
    If {ident}, points {p1[i]} and {p2[i]} will be the same point if {sgn} is {+1}.
    or related by a Y coordinate reversal if {sgn} is {-1}, with no noise. 
    
    Returns {np,p1,p2,w} in {*np_P,*p1_P,*p2_P,*w_P}.  Sometimes,
    the weight table {w} will be {NULL}. */

void hpot_check_opt_pmap
  ( hr2_pmap_t *M,
    hr2_pmap_type_t type,
    sign_t sgn,
    hr2_pmap_opt_func_t *f2,
    double f2Stop,
    double peps,
    double f2eps,
    bool_t verbose
  );
  /* If {f2(M)} does not exceed {f2Stop}, does nothing.
    Else evaluates {f2} at various perturbed versions {N} of {M},
    and fails if {M} does not look optimum at all -- specifically,
    if {f2(N) < f2(M) - f2eps}. The normof the perturbations, in 
    encoding space, will be {peps}.  */

void hpot_normalize_matrix(r3x3_t *A);
  /* Tries to scale all elements of {A} so that element {A[0][0]} becomes 1.
    If that element is too small, scales so that the sum of squares of elements is 1.
    if all elements are too small, leaves the matrix unchanged. */

void hpot_print_map(char *name, hr2_pmap_t *M, double f2M);
  /* Prints the map {M} (direct and inverse) to {stderr},
    prefixed by a line"  {name} = ".  if {f2M} is not {NAN}
    prints it too. */

r3x3_t hpot_throw_r3x3(double mag);
  /* Returns a random 3x3 matrix with entries in {[-mag _ +mag]}. */

void hpot_print_encoding(int32_t ny, double y[]);
  /* Prints {y[0..ny-1]} on one line, preceded by "encoding". */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srand(4615*417);
    
    char *xtype = argv[1];
    hr2_pmap_type_t type;
    demand(hr2_pmap_type_from_string(xtype, &type), "invalid map type option");

    for (int32_t i = 0; i < 10; i++)
      { bool_t verbose = (i < 3);
        char *outPrefix = NULL;
        if (i < 2) { asprintf(&outPrefix, "out/test-%03d", i); }
        test_hr2_pmap_opt_quadratic__hr2_pmap_opt_1D_plot(type, verbose, outPrefix);
        free(outPrefix);
      }
    return 0;
  }
  
void test_hr2_pmap_opt_quadratic__hr2_pmap_opt_1D_plot(hr2_pmap_type_t type, bool_t verbose, char *outPrefix)
  {
    char *xtype = hr2_pmap_type_to_string(type);
    fprintf(stderr, "> --- %s type = %s ---\n", __FUNCTION__, xtype);
    
    bool_t is_signed = ((type != hr2_pmap_type_IDENTITY) && (type != hr2_pmap_type_TRANSLATION));
                
    for (int32_t kident = TRUE; kident >= FALSE; kident--)
      { bool_t ident = (bool_t)kident;
        if (ident || (type != hr2_pmap_type_IDENTITY))
          { for (int32_t ktight = TRUE; ktight >= FALSE; ktight--)
              { bool_t tight = (bool_t)ktight;
                if (tight || (type != hr2_pmap_type_IDENTITY))
                  { for (int32_t sgn = +1; sgn >= -1; sgn -= 2)
                      { if (is_signed || (sgn == +1))
                          { hpot_do_test_opt_plot(type, sgn, tight, ident, verbose, outPrefix); }
                      }
                  }
              }
          }
      }

    fprintf(stderr, "< --- %s ---\n", __FUNCTION__);
  }
  
void hpot_do_test_opt_plot
  ( hr2_pmap_type_t type,
    sign_t sgn,
    bool_t tight,
    bool_t ident,
    bool_t verbose,
    char *outPrefix
  )
  {
    if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
    char *xtype = hr2_pmap_type_to_string(type);
    if (verbose) { fprintf(stderr, "  > --- %s type = %s sgn = %+d tight = %c ident = %c--- \n", __FUNCTION__, xtype, sgn, "FT"[tight], "FT"[ident]); }
    
    int32_t np;
    r2_t *p1 = NULL;
    r2_t *p2 = NULL;
    double *w = NULL;
    hpot_choose_r2_point_pairs(type, sgn, tight, ident, &np, &p1, &p2, &w);
    if (verbose) 
      { fprintf(stderr, "    testing with %d point pairs\n", np);
        for (int32_t kp = 0; kp < np; kp++)
          { fprintf(stderr, "      %03d", kp);
            r2_gen_print(stderr, &(p1[kp]), "%+14.6f", " [ ", " ", " ]");
            r2_gen_print(stderr, &(p2[kp]), "%+14.6f", " [ ", " ", " ]");
            if (w != NULL) { fprintf(stderr, " w = %20.16f", w[kp]); }
            fprintf(stderr, "\n");
          }
      }
    
    hr2_pmap_t M = hr2_pmap_throw_type(type, sgn); /* Initial guess matrix. */
      
    /* Compute max distance between paired points mapped  by the initial guess: */
    double maxDist = 0.0;
    for (int32_t kp = 0; kp < np; kp++)
      { r2_t *p1k = &(p1[kp]);
        r2_t *p2k = &(p2[kp]);
        r2_t q1 = hr2_pmap_r2_point(p1k, &M);
        double dk12 = r2_dist(&q1, p2k);
        if (dk12 > maxDist) { maxDist = dk12; }
        r2_t q2 = hr2_pmap_inv_r2_point(p2k, &M);
        double dk21 = r2_dist(&q2, p1k);
        if (dk21 > maxDist) { maxDist = dk21; }
      }
    if (verbose) { fprintf(stderr, "    maxDist = %20.14f\n", maxDist); }
    
    auto double goal_func(hr2_pmap_t *P);
    
    int32_t maxIter = 20;
    double maxMod = 5.0*fmax(1.0, maxDist);
    double f2Stop = 0.000001;
    if (verbose) { fprintf(stderr, "    maxIter = %d maxMod = %20.14f  f2Stop = %20.14f\n", maxIter, maxMod, f2Stop); }
    
    double f2M = goal_func(&M);
    if (verbose) { hpot_print_map("initial", &M, f2M); }
    if ((outPrefix != NULL) && (type != hr2_pmap_type_IDENTITY))
      { hpot_do_test_1D_plot(outPrefix, type, sgn, tight, ident, &M, goal_func, maxMod, "ini", verbose); }
    
    hr2_pmap_opt_quadratic(type, sgn, &goal_func, maxIter, maxMod, f2Stop, &M, &f2M, verbose);
    
    if (verbose) { hpot_print_map("final", &M, f2M); }
    if ((outPrefix != NULL) && (type != hr2_pmap_type_IDENTITY))
      { hpot_do_test_1D_plot(outPrefix, type, sgn, tight, ident, &M, goal_func, 0.01*maxMod, "opt", verbose); }
    
    if (tight) { demand(f2M <= 2*f2Stop, "{hr2_pmap_opt_quadratic} did not find exact solution"); }
    
    if (tight && ident)
      { /* Solution should be the identity matrix: */
        hr2_pmap_t I = hr2_pmap_identity();
        hr2_pmap_set_sign(&I, sgn);
        double d2 = hr2_pmap_diff_sqr(&M, &I);
        demand (d2 <= 1.0e-6, "solution should be the identity map");
      }

    double peps = 0.001;      /* Size of perturbations to encoded map. */
    double f2eps = 10*f2Stop; /* Tolerance for value of {f2} being minimum. */
    hpot_check_opt_pmap(&M, type, sgn, &goal_func, f2Stop, peps, f2eps, verbose);
    
    return;
    
    double goal_func(hr2_pmap_t *P)
      { return hr2_pmap_mismatch_sqr(P, np, p1, p2, w); }

    if (verbose) { fprintf(stderr, "  < --- %s ---\n", __FUNCTION__); }
    if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
  }

void hpot_choose_r2_point_pairs
  ( hr2_pmap_type_t type,
    sign_t sgn,
    bool_t tight,
    bool_t ident,
    int32_t *np_P,
    r2_t **p1_P,
    r2_t **p2_P,
    double **w_P
  )
  { 
    int32_t nr = hr2_pmap_from_many_pairs_required_rank(type);
    int32_t np = (tight? nr : nr + 10);
    r2_t *p1 = talloc(np, r2_t);
    r2_t *p2 = talloc(np, r2_t);
    double *w = (drandom() < 0.2 ? NULL : talloc(np, double));
    
    /* Get maps {M1,M2} with the right type and sign: */
    hr2_pmap_t M1 = hr2_pmap_throw_type((ident ? hr2_pmap_type_IDENTITY : type), +1);
    hr2_pmap_t M2 = hr2_pmap_throw_type((ident ? hr2_pmap_type_IDENTITY : type), sgn);
    
    /* Generate points {p1[0..np-1]} and weights {w[0..np-1]}: */
    for (int32_t kp = 0; kp < np; kp++)
      { /* Fill {p1} with {nr} well-spaced points and some random points: */
        r2_t q;
        if (kp < nr)
          { double x = ((kp&1) != 0 ? -1.0 : +1.0);
            double y = ((kp&2) != 0 ? -1.0 : +1.0);
            q = (r2_t){{ x, y }};
          }
        else
          { q = (r2_t){{ 10.0*(2*drandom()-1),  10.0*(2*drandom()-1) }}; }
        p1[kp] = q;
        p2[kp] = q;
        /* Set {p1,p2} to the images of {p1} by {M1,M2} plus noise: */
        p1[kp] = hr2_pmap_r2_point(&(p1[kp]), &M1);
        p2[kp] = hr2_pmap_r2_point(&(p2[kp]), &M2);
        if ((! ident) && ((! tight) || (kp > nr)))
          { /* Add random noise, small enough to preserve the handedness: */
            r2_t e1 = (r2_t){{ 0.1 * (drandom()-0.5),  0.1 * (drandom()-0.5) }};
            r2_add(&e1, &(p1[kp]), &(p1[kp]));
            r2_t e2 = (r2_t){{ 0.1 * (drandom()-0.5),  0.1 * (drandom()-0.5) }};
            r2_add(&e2, &(p2[kp]), &(p2[kp]));
            /* Set the weight to random, less than 1: */
            if (w != NULL) { w[kp] = fmax(0.0, 1.01*(drandom()*drandom()) - 0.01); }
          }
        else
          { /* Set the weight to random, large, less than 1: */
            if (w != NULL) { w[kp] = 0.7 + 0.3*drandom(); }
          }
      }
    (*np_P) = np;
    (*p1_P) = p1;
    (*p2_P) = p2;
    (*w_P) = w;
  }

void hpot_do_test_1D_plot
  ( char *outPrefix,
    hr2_pmap_type_t type,
    sign_t sgn,
    bool_t tight,
    bool_t ident,
    hr2_pmap_t *M,
    hr2_pmap_opt_func_t *f2,
    double urad,
    char *stage,
    bool_t verbose
  ) 
  { 
    if (verbose) { fprintf(stderr, "    > --- %s ---\n", __FUNCTION__); }
    char *xtype = hr2_pmap_type_to_string(type);
    char *xsgn = (sgn == 0 ? "o" : (sgn < 0 ? "m" : "p"));
    char *tag = NULL;
    asprintf(&tag, "%s-%s-tight%c-ident%c-%s", xtype, xsgn, "FT"[tight], "FT"[ident], stage);
    int32_t nu = 15;
    int32_t ns = 50;
    hr2_pmap_opt_1D_plot(outPrefix, tag, M, type, sgn, f2, nu, urad, ns);
    free(tag);

    if (verbose) { fprintf(stderr, "    < --- %s ---\n", __FUNCTION__); }
  }

void hpot_check_opt_pmap
  ( hr2_pmap_t *M,
    hr2_pmap_type_t type,
    sign_t sgn,
    hr2_pmap_opt_func_t *f2,
    double f2Stop,
    double peps,
    double f2eps,
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "    > --- %s ---\n", __FUNCTION__); }
    if (verbose) { fprintf(stderr, "    f2Stop = %.18f  peps = %.18f  f2eps = %.18f\n", f2Stop, peps, f2eps); }
    double f2M = f2(M);
    if (verbose) { fprintf(stderr, "    f2(M) = %.12f = %24.16e\n", f2M, f2M); }
    if (f2M <= f2Stop) 
      { /* Good enough: */ 
        if (verbose) { fprintf(stderr, "    M is good enough!\n"); }
      }
    else
      { if (verbose) { fprintf(stderr, "    M is not good enough; probng around it\n"); }
        int32_t ny = hr2_pmap_encode_num_parameters(type);
        if (ny > 0)
          { double y[ny];  /* Encoding of alleged optimum matrix. */
            hr2_pmap_encode(M, type, ny, y);

            double u[ny];   /* Perturbation direction. */
            double yt[ny];  /* Perturbed encoding vector. */
            hr2_pmap_t N;   /* Perturbed map */
            int32_t nd = 2*ny;  /* Number of directions for probing. */
            bool_t ok = TRUE;
            hr2_pmap_t N_best;
            double f2N_best = f2M;
            for (int32_t kd = 0; kd < nd; kd++)
              { rn_throw_dir(ny, u);
                for (int32_t ky = 0; ky < ny; ky++)
                  { yt[ky] = y[ky] + peps*u[ky]; }
                hr2_pmap_decode(ny, yt, type, sgn, &N);
                double f2N = f2(&N);
                if (f2N < f2M - f2eps)
                  { if (verbose) 
                      { fprintf(stderr, "    oops, found a better matrix:\n");
                        hpot_print_map("N", &N, f2N);
                      }
                    ok = FALSE;
                    if (f2N < f2N_best) { f2N_best = f2N; N_best = N; }
                  }
              }
            if (! ok) 
              { fprintf(stderr, "    M is NOT OPTIMAL\n");
                hpot_print_map("N_best", &N_best, f2N_best);
              }
            else
              { if (verbose) { fprintf(stderr, "    M seems to be optimal\n"); } }
          }
      }
    if (verbose) { fprintf(stderr, "    < --- %s ---\n", __FUNCTION__); }
  }
  
void hpot_normalize_matrix(r3x3_t *A)
  { double norm = r3x3_norm(A);
    double head = fabs(A->c[0][0]);
    double scale = (head > 0.01*norm ? head : norm); 
    if ((scale == 1.0) || (scale < 1.0e-200)) { return; }
    for(int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++) 
          { A->c[i][j] /= scale; }
      }
  }

r3x3_t hpot_throw_r3x3(double mag)
  { r3x3_t A;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { A.c[i][j] = mag*(2*drandom() - 1.0); }
      }
    return A;
  }

void hpot_print_map(char *name, hr2_pmap_t *M, double f2M)
  {
    fprintf(stderr, "  %s", name);
    if (! isnan(f2M)) { fprintf(stderr, " ( f2 = %18.10f )", f2M); }
    fprintf(stderr, " =\n");
    hr2_pmap_gen_print(stderr, M, "%12.7f", "    ", "[ ","  "," ]\n    ", "[ "," "," ]", "\n");
  }

void hpot_print_encoding(int32_t ny, double y[])
  { 
    fprintf(stderr, "  encoding factors:\n");
    fprintf(stderr, "    ");
    for (int32_t kv = 0; kv < ny; kv++)
      { fprintf(stderr, " %+12.8f", y[kv]); }
    fprintf(stderr, "\n");
  }
