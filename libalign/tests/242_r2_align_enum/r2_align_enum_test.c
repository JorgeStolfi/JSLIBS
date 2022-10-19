#define PROG_NAME "r2_align_enum_test"
#define PROG_DESC "test of {r2_align_enum.h}"
#define PROG_VERS "1.0"

/* Last edited on 2022-02-28 13:12:30 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_align_COPYRIGHT \
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

#include <r2_align_enum.h>

int32_t main(int32_t argn, char **argv);

void ralent_do_test(int32_t ni);
  /* Performs a test with alignment vectors of {ni} points. */

void ralent_choose_arad(int32_t ni, r2_t arad[]);
  /* Stores into {arad[0..ni-1]} a random radius vector for the basic ellipsoid {\RE}. */

void ralent_choose_ctr(int32_t ni, r2_t ctr[]);
  /* Stores into {ctr[0..ni-1]} a random center for the basic ellipsoid {\RE}. */

void ralent_test_align_enum(int32_t ni, r2_t ctr[], r2_t arad[], double tol, r2_t u0[], r2_t u1[], r2_t popt[]);
  /* Tests the {r2_align_enum} function in the ellipsoid {\RF} defined
    by center {ctr[0..ni-1]} and radius vector {arad[0..ni-1]}. Searches
    on a grid with step {tol}. The goal function is the squared distance
    from {popt[0..ni-1]}.
    
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

    ralent_do_test(2);
    ralent_do_test(5);
    
    return 0;
  }
  
void ralent_do_test(int32_t ni)
  { 
    fprintf(stderr, "--- testing {r2_align_enum} ni = %d ---\n", ni);

    fprintf(stderr, "... choosing ellipsoid center {ctr} ...\n");
    r2_t ctr[ni];   /* Central alignment vector. */
    ralent_choose_ctr(ni, ctr);

    fprintf(stderr, "... choosing basic ellipsoid radius {arad} ...\n");
    r2_t arad[ni];  /* Radius vector if basic ellipsoid {\RE}. */
    ralent_choose_arad(ni, arad);
    
    fprintf(stderr, "... Computing the main axes and radii of the search ellipsoid {\\RF} ...\n");
    int32_t nd = r2_align_count_degrees_of_freedom(ni, arad);
    fprintf(stderr, "search dimensions nd = %d\n", nd);
    r2_t U[ni*nd];
    double urad[nd];
    r2_align_compute_search_ellipsoid (ni, arad, nd, U, urad);

    fprintf(stderr, "... Finding the largest dimension of {\\RF} ...\n");
    double ursup = 0.0;
    for (int32_t k = 0; k < nd; k++)  { if (urad[k] > ursup) { ursup = urad[k]; } }
    fprintf(stderr, "largest radius of {\\RF} = %.8f\n", ursup);

    fprintf(stderr, "... Choosing the optimum point {popt} in {\\RF} ...\n");
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
            for (int32_t i = 0; i < ni; i++) 
              { for (int32_t j = 0; j < 2; j++) 
                  { double dij = 0.0;
                    for (int32_t k = 0; k < nd; k++) 
                      { r2_t *uk = &(U[k*ni]);
                        dij += b[k]*urad[k]*uk[i].c[j];
                      }
                    popt[i].c[j] = ctr[i].c[j] + dij;
                  }
              }
            ntry++;
            dr = sqrt(r2_align_rel_dist_sqr(ni, popt, ctr, arad));
            de = sqrt(r2_align_dist_sqr(ni, popt, ctr));
            fprintf(stderr, "  dr = %.8f  de = %.8f\n", dr, de);
            /* Check if well within the ellipsoid and not too close to {ctr}: */
            if ((dr <= drmax) && (de >= demin)) { break; }
          }
        fprintf(stderr, "chose {popt} with %d attempts dr = %.8f de = %.8f\n", ntry, dr, de);
        r2_align_print_vector(stderr, ni, "popt", -1, popt, FALSE);
        demand(dr <= 1.0 + 1.0e-8, "{popt} outside the ellipsoid");
      }
    else
      { /* Use the center: */
        for (int32_t i = 0; i < ni; i++) { popt[i] = ctr[i]; }
      }

    /* If {nd} is 2, define the axis vectors of the indep plot variables: */
    r2_t *u0 = (nd == 2 ? &(U[0*ni]) : NULL);
    r2_t *u1 = (nd == 2 ? &(U[1*ni]) : NULL);
    ralent_test_align_enum(ni, ctr, arad, tol, u0, u1, popt);
    return;
  }
    
void ralent_choose_arad(int32_t ni, r2_t arad[])
  { 
    fprintf(stderr, "... choosing the basic domain radius {arad} ...\n");
    double rmax = 4.999;
    double rmin = 1.500;
    r2_t zfrac = (ni == 2 ? (r2_t){{ 0.00, 0.00 }} : (r2_t){{ 0.25, 0.75 }});
    r2_align_throw_arad(ni, zfrac, rmin, rmax, arad);
    r2_align_print_vector(stderr, ni, "arad", -1, arad, TRUE);
    return;
  }  
  
void ralent_choose_ctr(int32_t ni, r2_t ctr[])
  { 
    fprintf(stderr, "... choosing the center {ctr} ...\n");
    /* r2_align_throw_ball_vector(ni, 0.0, 1.995, ctr);  */
    for (int32_t i = 0; i < ni; i++) { ctr[i] = (r2_t){{ 1.0, 2.0 }}; }
    r2_align_print_vector(stderr, ni, "ctr", -1, ctr, FALSE);
    return;
  }
    
void ralent_test_align_enum(int32_t ni, r2_t ctr[], r2_t arad[], double tol, r2_t u0[], r2_t u1[], r2_t popt[])
  {
    bool_t debug = FALSE;
    
    int32_t nd = r2_align_count_degrees_of_freedom(ni, arad);
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
    for (int32_t i = 0; i < ni; i++) { psol[i] = ctr[i]; }
    double F2ini = FD2(ni, psol);
    fprintf(stderr, "  F2 (ini) = %12.6f\n", F2ini);

    double F2opt = FD2(ni, popt);
    fprintf(stderr, "  F2 (opt) = %12.6f\n", F2opt);
    
    double F2sol = NAN;
    plot = TRUE;
    r2_align_enum(ni, &FD2, arad, tol, psol, &F2sol);
    
    fprintf(stderr, "  F2 (sol) = %12.6f\n", F2sol);
    for (int32_t i = 0; i < ni; i++) 
      { fprintf(stderr, "  psol[%d] = (", i);
        for (int32_t j = 0; j < 2; j++) 
          { fprintf(stderr, " %12.6f", psol[i].c[j]); }
        fprintf(stderr, " ) popt = (");
        for (int32_t j = 0; j < 2; j++) 
          { fprintf(stderr, " %12.6f", popt[i].c[j]); }
        fprintf(stderr, " ) diff = (");
        for (int32_t j = 0; j < 2; j++) 
          { fprintf(stderr, " %12.6f", psol[i].c[j] - popt[i].c[j]); }
        fprintf(stderr, " )\n");
      }

    double dr = sqrt(r2_align_rel_dist_sqr(ni, psol, ctr, arad));
    demand(dr < 1.0 + 1.0e-8, "solution outside the ellipsoid");
    
    double de = sqrt(r2_align_dist_sqr(ni, psol, popt));
    fprintf(stderr, "error = %12.6f (%12.8f * tol)\n", de, de/tol);
    double deexp = tol*sqrt(nd)/2;  /* Assumes orthogonal grid with step {tol} in search space. */
    fprintf(stderr, "expected max %12.6f (%12.8f * tol)\n", deexp, deexp/tol);

    if (wr != NULL) { fclose(wr); }
    
    return;
    
    double FD2(int32_t ni, r2_t q[])
      { if (debug) { r2_align_print_vector(stderr, ni, "  psmp", -1, q, FALSE); }
        double F2val = r2_align_dist_sqr(ni, q, popt);
        if (plot && (wr != NULL))
          { r2_t dsmp[ni];
            for (int32_t i = 0; i < ni; i++) 
              { for (int32_t j = 0; j < 2; j++) 
                  { dsmp[i].c[j] = q[i].c[j] - ctr[i].c[j]; }
              }
            double s0 = r2_align_dot(ni, dsmp, u0);
            double s1 = r2_align_dot(ni, dsmp, u1);
            fprintf(wr, "%12.6f %12.6f %12.6f\n", s0, s1, F2val);
          }
        return F2val;
      }
  }
