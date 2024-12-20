#define PROG_NAME "test_r2_align_enum"
#define PROG_DESC "test of {r2_align_enum.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-05 10:20:18 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_align_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

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

void ralent_do_test(int32_t ni, bool_t bal);
  /* Performs a test with alignment vectors of {ni} points. 
    The alignments will be balanced if {bal} is true.. */

void ralent_test_align_enum(int32_t ni, r2_t ctr[], r2_t arad[], bool_t bal, double tol, r2_t u0[], r2_t u1[], r2_t popt[]);
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

void ralent_choose_arad(int32_t ni, r2_t arad[]);
void ralent_choose_ctr(int32_t ni, r2_t ctr[]);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srandom(4615*417);

    for (uint32_t bal = 0;  bal < 2; bal++)
      { ralent_do_test(2, (bool_t)bal);
        ralent_do_test(5, (bool_t)bal);
      }
    
    return 0;
  }
  
void ralent_do_test(int32_t ni, bool_t bal)
  { 
    fprintf(stderr, "--- testing {r2_align_enum} ni = %d ---\n", ni);

    fprintf(stderr, "... choosing ellipsoid center {ctr} ...\n");
    r2_t ctr[ni];   /* Central alignment vector. */
    double cmax = 5.000;
    r2_align_throw_ctr(ni, cmax, ctr, TRUE);

    fprintf(stderr, "... choosing basic ellipsoid radius {arad} ...\n");
    r2_t arad[ni];  /* Radius vector if basic ellipsoid {\RE}. */
    double rmax = 3.000;
    int32_t nvmin = 1;
    r2_align_throw_arad(ni, rmax, nvmin, arad, TRUE);    
    
    fprintf(stderr, "... Computing the dimension of the search ellipsoid {\\RF} ...\n");
    i2_t nv = r2_align_count_variable_coords (ni, arad);
    fprintf(stderr, "num of variable coords nv = (%d,%d)\n", nv.c[0], nv.c[1]);
    int32_t nd = r2_align_count_degrees_of_freedom(nv, bal);
    fprintf(stderr, "search dimensions nd = %d\n", nd);
    
    fprintf(stderr, "... Computing the main axes and radii of the search ellipsoid {\\RF} ...\n");
    r2_t U[ni*nd];
    double urad[nd];
    r2_align_compute_search_ellipsoid (ni, arad, bal, nd, U, urad);

    fprintf(stderr, "... Finding the largest dimension of {\\RF} ...\n");
    double ursup = 0.0;
    for (uint32_t k = 0;  k < nd; k++)  { if (urad[k] > ursup) { ursup = urad[k]; } }
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
            fprintf(stderr, "  dr = %.8f  de = %.8f\n", dr, de);
            /* Check if well within the ellipsoid and not too close to {ctr}: */
            if ((dr <= drmax) && (de >= demin)) { break; }
          }
        fprintf(stderr, "chose {popt} with %d attempts dr = %.8f de = %.8f\n", ntry, dr, de);
        r2_align_print_vector(stderr, ni, "popt", -1, popt);
        demand(dr <= 1.0 + 1.0e-8, "{popt} outside the ellipsoid");
      }
    else
      { /* Use the center: */
        for (uint32_t i = 0;  i < ni; i++) { popt[i] = ctr[i]; }
      }

    /* If {nd} is 2, define the axis vectors of the indep plot variables: */
    r2_t *u0 = (nd == 2 ? &(U[0*ni]) : NULL);
    r2_t *u1 = (nd == 2 ? &(U[1*ni]) : NULL);
    ralent_test_align_enum(ni, ctr, arad, bal, tol, u0, u1, popt);
    return;
  }
    
void ralent_choose_arad(int32_t ni, r2_t arad[])
  { 
    double rmax = 4.999;
    int32_t nvmin = 1;
    
    r2_align_throw_arad(ni, rmax, nvmin, arad, TRUE);
    return;
  }  
  
void ralent_choose_ctr(int32_t ni, r2_t ctr[])
  { 
    double cmax = 1.995;
    fprintf(stderr, "... choosing the center {ctr} ...\n");
    r2_align_throw_ctr(ni, cmax, ctr, TRUE);
    for (uint32_t i = 0;  i < ni; i++) { ctr[i] = (r2_t){{ 1.0, 2.0 }}; }
    r2_align_print_vector(stderr, ni, "ctr", -1, ctr);
    return;
  }
    
void ralent_test_align_enum(int32_t ni, r2_t ctr[], r2_t arad[], bool_t bal, double tol, r2_t u0[], r2_t u1[], r2_t popt[])
  {
    bool_t debug = FALSE;
    
    i2_t nv = r2_align_count_variable_coords (ni, arad);
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
    r2_align_enum(ni, &FD2, arad, bal, tol, psol, &F2sol);
    
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
