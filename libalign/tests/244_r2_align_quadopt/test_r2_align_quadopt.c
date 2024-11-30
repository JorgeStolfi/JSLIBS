#define PROG_NAME "test_r2_align_quadopt"
#define PROG_DESC "test of {r2_align_quadopt.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-08 11:19:23 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test__r2_align_quadopt_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <ix.h>
#include <r2.h>
#include <rn.h>
#include <bool.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>
#include <r2_align.h>

#include <r2_align_quadopt.h>

int32_t main(int32_t argn, char **argv);

void talq_do_test(int32_t ni, bool_t bal, bool_t verbose);
  /* Performs a test with alignment vectors of {ni} points,
    required to be balanced or not depending on {bal}. */

void talq_choose_arad(int32_t ni, r2_t arad[], bool_t verbose);
  /* Stores into {arad[0..ni-1]} a random radius vector for the basic ellipsoid {\RE}. */

void talq_choose_ctr(int32_t ni, r2_t ctr[], bool_t verbose);
  /* Stores into {ctr[0..ni-1]} a random center for the basic ellipsoid {\RE}. */

void talq_test_align_quadopt
  ( int32_t ni,
    r2_t ctr[],
    r2_t arad[],
    bool_t bal,
    double tol,
    r2_t u0[],
    r2_t u1[],
    r2_t popt[],
    bool_t verbose
  );
  /* Tests the {r2_align_quadopt} function in the ellipsoid {\RF} defined
    by center {ctr[0..ni-1]}, radius vector {arad[0..ni-1]}, and balancing
    option {bal}. Searches on a grid with step {tol}. The goal function is 
    the squared distance from {popt[0..ni-1]}. 
    
    If {u0,u1} are not NULL, they should be unit vectors in the linear
    subspace {\RV} parallel to {\RF}. Then it also writes a file
    "out/f.dat" with 3D plot data of the goal function, with one line
    for each sample point {psmp} in the search ellipsoid, whose indep
    variables are the projections of {psmp-ctr} along those two
    vectors */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srandom(4615*417);

    bool_t verbose = TRUE;
    for (uint32_t kbal = 0;  kbal <= 1; kbal++)
      { bool_t bal = kbal > 0;
        talq_do_test(2, bal, verbose);
        talq_do_test(5, bal, verbose);
      }
    
    return 0;
  }
  
void talq_do_test(int32_t ni, bool_t bal, bool_t verbose)
  { 
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "--- testing {r2_align_quadopt} ni = %d ---\n", ni);

    if (verbose) { fprintf(stderr, "... choosing ellipsoid center {ctr} ...\n"); }
    r2_t ctr[ni];   /* Central alignment vector. */
    talq_choose_ctr(ni, ctr, verbose);

    if (verbose) { fprintf(stderr, "... choosing basic ellipsoid radius {arad} ...\n"); }
    r2_t arad[ni];  /* Radius vector if basic ellipsoid {\RE}. */
    talq_choose_arad(ni, arad, verbose);
    
    if (verbose) { fprintf(stderr, "... Computing the main axes and radii of the search ellipsoid {\\RF} ...\n"); }
    i2_t nv = r2_align_count_variable_coords (ni, arad);
    if (verbose) { fprintf(stderr, "  num of variable coords nv = (%d,%d)\n", nv.c[0], nv.c[1]); }
    int32_t nd = r2_align_count_degrees_of_freedom(nv, bal);
    if (verbose) { fprintf(stderr, "  search dimensions nd = %d\n", nd); }
    r2_t U[ni*nd];
    double urad[nd];
    r2_align_compute_search_ellipsoid (ni, arad, bal, nd, U, urad);

    if (verbose) { fprintf(stderr, "... Finding the largest dimension of {\\RF} ...\n"); }
    double ursup = 0.0;
    for (uint32_t k = 0;  k < nd; k++)  { if (urad[k] > ursup) { ursup = urad[k]; } }
    if (verbose) { fprintf(stderr, "largest radius of {\\RF} = %.8f\n", ursup); }

    if (verbose) { fprintf(stderr, "... Choosing the optimum point {popt} in {\\RF} ...\n"); }
    double tol = 0.30;
    r2_t popt[ni];
    if (nd > 0)
      { /* Choose a point away from center mut not too close to edge of {\RF}: */
        double demin = 2.0*tol;    /* Min Euclidean distance from {ctr} to {popt}. */
        demand(demin <= 0.5*ursup, "{tol} too big for ellipsoid {\\RF}");
        double drmax = 0.75;       /* Max distance from {ctr} rel to {\RF}. */
        int32_t ntry = 0;
        double de = NAN; /* Euclidean distance from {ctr} to {popt}. */
        double dr = NAN; /* Distance from {ctr} to {popt} relative to ellipsoid {\RE}. */
        while (TRUE)
          { assert(ntry < 1000); /* Should not happen... */
            double b[nd]; /* A random vector in the unit ball. */
            rn_throw_ball(nd, b);
            for (uint32_t i = 0;  i < ni; i++) 
              { for (uint32_t j = 0;  j < 2; j++) 
                  { double dij = 0.0;
                    for (uint32_t k = 0;  k < nd; k++) 
                      { r2_t *uk = &(U[k*ni]);
                        dij += b[k]*urad[k]*uk[i].c[j];
                      }
                    popt[i].c[j] = ctr[i].c[j] + dij;
                  }
              }
            ntry++;
            dr = sqrt(r2_align_rel_disp_sqr(ni, popt, ctr, arad));
            de = sqrt(r2_align_dist_sqr(ni, popt, ctr));
            if (verbose) { fprintf(stderr, "  dr = %.8f  de = %.8f\n", dr, de); }
            /* Check if well within the ellipsoid and not too close to {ctr}: */
            if ((dr <= drmax) && (de >= demin)) { break; }
          }
        if (verbose) { fprintf(stderr, "chose {popt} with %d attempts dr = %.8f de = %.8f\n", ntry, dr, de); }
        if (verbose) { r2_align_print_vector(stderr, ni, "popt", -1, popt); }
        demand(dr <= 1.0 + 1.0e-8, "{popt} outside the ellipsoid");
      }
    else
      { /* Use the center: */
        for (uint32_t i = 0;  i < ni; i++) { popt[i] = ctr[i]; }
      }

    /* If {nd} is 2, define the axis vectors of the indep plot variables: */
    r2_t *u0 = (nd == 2 ? &(U[0*ni]) : NULL);
    r2_t *u1 = (nd == 2 ? &(U[1*ni]) : NULL);
    talq_test_align_quadopt(ni, ctr, arad, bal, tol, u0, u1, popt, verbose);
    return;
  }
    
void talq_choose_arad(int32_t ni, r2_t arad[], bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "... choosing the basic domain radius {arad} ...\n"); }
    double rmax = 4.999;
    int32_t nvmin = 2;
    r2_align_throw_arad(ni, rmax, nvmin, arad, verbose);
    return;
  }  
  
void talq_choose_ctr(int32_t ni, r2_t ctr[], bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "... choosing the center {ctr} ...\n"); }
    /* r2_align_throw_ball_vector(ni, 0.0, 1.995, ctr);  */
    for (uint32_t i = 0;  i < ni; i++) { ctr[i] = (r2_t){{ 1.0, 2.0 }}; }
    if (verbose) { r2_align_print_vector(stderr, ni, "ctr", -1, ctr); }
    return;
  }
    
void talq_test_align_quadopt
  ( int32_t ni,
    r2_t ctr[],
    r2_t arad[],
    bool_t bal,
    double tol,
    r2_t u0[],
    r2_t u1[],
    r2_t popt[],
    bool_t verbose
  )
  {
    bool_t debug = FALSE;
    
    /* Count variable coords along each axis: */
    i2_t nv = r2_align_count_variable_coords (ni, arad); 
    
    /* Determine the dimension of the search space: */
    int32_t nd = r2_align_count_degrees_of_freedom(nv, bal);
    FILE *wr = NULL;
    
    /* Find index {ip[j]} of alignment coord to plot from each axis, or {-1} if none: */
    if ((u0 != NULL) && (u1 != NULL))  
      { assert(nd == 2);
        wr = open_write("out/f.dat", TRUE);
      }
      
    bool_t plot = FALSE;
    
    auto double FD2(int32_t ni, r2_t q[]);
      /* Quadratic distance from {q} to {popt}.  If {plot} is true and {wr} is not {NULL},
        also writes the data points to {wr}. */
      
    r2_t psol[ni];
    for (uint32_t i = 0;  i < ni; i++) { psol[i] = ctr[i]; }
    double F2ini = FD2(ni, psol);
    fprintf(stderr, "  F2 (ini) = %12.6f\n", F2ini);

    double F2opt = FD2(ni, popt);
    fprintf(stderr, "  F2 (opt) = %12.6f\n", F2opt);
    
    double F2sol = NAN;
    plot = TRUE;
    r2_align_quadopt(ni, &FD2, arad, bal, tol, psol, &F2sol);
    
    fprintf(stderr, "  F2 (sol) = %12.6f\n", F2sol);
    for (uint32_t i = 0;  i < ni; i++) 
      { fprintf(stderr, "  psol[%d] = (", i);
        for (uint32_t j = 0;  j < 2; j++) 
          { fprintf(stderr, " %12.6f", psol[i].c[j]); }
        fprintf(stderr, " ) popt = (");
        for (uint32_t j = 0;  j < 2; j++) 
          { fprintf(stderr, " %12.6f", popt[i].c[j]); }
        fprintf(stderr, " ) diff = (");
        for (uint32_t j = 0;  j < 2; j++) 
          { fprintf(stderr, " %12.6f", psol[i].c[j] - popt[i].c[j]); }
        fprintf(stderr, " )\n");
      }

    double dr = sqrt(r2_align_rel_disp_sqr(ni, psol, ctr, arad));
    demand(dr < 1.0 + 1.0e-8, "solution outside the ellipsoid");
    
    double de = sqrt(r2_align_dist_sqr(ni, psol, popt));
    fprintf(stderr, "error = %12.6f (%12.8f * tol)\n", de, de/tol);
    double deexp = tol*sqrt(nd)/2;  /* Assumes orthogonal grid with step {tol} in search space. */
    fprintf(stderr, "expected max %12.6f (%12.8f * tol)\n", deexp, deexp/tol);

    if (wr != NULL) { fclose(wr); }
    
    return;
    
    double FD2(int32_t ni, r2_t q[])
      { if (debug) { r2_align_print_vector(stderr, ni, "  psmp", -1, q); }
        double F2val = r2_align_dist_sqr(ni, q, popt);
        if (plot && (wr != NULL))
          { r2_t dsmp[ni];
            for (uint32_t i = 0;  i < ni; i++) 
              { for (uint32_t j = 0;  j < 2; j++) 
                  { dsmp[i].c[j] = q[i].c[j] - ctr[i].c[j]; }
              }
            double s0 = r2_align_dot(ni, dsmp, u0);
            double s1 = r2_align_dot(ni, dsmp, u1);
            fprintf(wr, "%12.6f %12.6f %12.6f\n", s0, s1, F2val);
          }
        return F2val;
      }
  }
