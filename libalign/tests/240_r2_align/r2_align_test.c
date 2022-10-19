#define PROG_NAME "r2_align_test"
#define PROG_DESC "test of {r2_align.h}"
#define PROG_VERS "1.0"

/* Last edited on 2022-02-28 13:11:45 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_align_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <r2.h>
#include <bool.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>
#include <float_image.h>

#include <r2_align.h>

int32_t main(int32_t argn, char **argv);

void ralt_test_rel_dist_sqr(int32_t ni);
  /* Tests {r2_align_rel_dist_sqr} for alignment vectors of {ni} elements.
    The parameter {ni} must be 2 or more. */

void ralt_test_compute_search_ellipsoid(int32_t ni);
  /* Tests {r2_align_compute_search_ellipsoid} for alignment vectors of {ni} elements. */

void ralt_choose_arad(int32_t ni, r2_t arad[]);
  /* Stores into {arad[0..ni-1]} a random reference radius vector. */

void ralt_choose_ctr(int32_t ni, r2_t ctr[]);
  /* Stores into {ctr[0..ni-1]} a random alignment vector to be the center of the basic 
    domain ellipsoid {\RE}. */

void ralt_plot_rel_dist_sqr(int32_t ni, r2_t ctr[], r2_t arad[]);
  /* Writes a  file "out/f2.dat" with a random 2D slice of the function
    {r2_align_rel_dist_sqr} over the ellipsoid with semi-axes {arad}. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srandom(4615*417);

    ralt_test_rel_dist_sqr(5);

    ralt_test_compute_search_ellipsoid(1);
    ralt_test_compute_search_ellipsoid(2);
    ralt_test_compute_search_ellipsoid(3);
    ralt_test_compute_search_ellipsoid(5);
    ralt_test_compute_search_ellipsoid(7);
    return 0;
  }
  
void ralt_test_rel_dist_sqr(int32_t ni)
  { 
    fprintf(stderr, "-- testing {r2_align_rel_dist_sqr} ni = %d ---\n", ni);

    r2_t arad[ni];  /* Search radius for each coordinate. */
    ralt_choose_arad(ni, arad);

    r2_t ctr[ni];  /* Initial alignment. */
    ralt_choose_ctr(ni, ctr);

    ralt_plot_rel_dist_sqr(ni, ctr, arad);
    return;
  }
    
void ralt_choose_arad(int32_t ni, r2_t arad[])
  { 
    fprintf(stderr, "... choosing the basic domain radius {arad} ...\n");
    double rmax = 4.999;
    double rmin = 1.500;
    r2_t zfrac = (r2_t){{ 0.25, 0.75 }};
    r2_align_throw_arad(ni, zfrac, rmin, rmax, arad);
    r2_align_print_vector(stderr, ni, "arad", -1, arad, TRUE);
    return;
  }  

void ralt_choose_ctr(int32_t ni, r2_t ctr[])
  { 
    fprintf(stderr, "... choosing the center alighnment {ctr} ...\n");
    r2_align_throw_ball_vector(ni, 0.0, 9.995, ctr);
    r2_align_print_vector(stderr, ni, "ctr", -1, ctr, FALSE);
    return;
  }

void ralt_test_compute_search_ellipsoid(int32_t ni)
  {
    bool_t debug = TRUE;
    
    fprintf(stderr, "--- testing {r2_align_compute_search_ellipsoid} ni = %d ---\n", ni);

    /* Choose the radius vector {arad} of the basic domain ellipsoid {\RE}: */
    r2_t arad[ni];  /* Search radius for each coordinate. */
    ralt_choose_arad(ni, arad);

    /* Compute the eigenvectors {U} and eigenvalues {urad} of the search ellipsoid {\RF}: */
    fprintf(stderr, "... computing the search ellipsoid basis {U,urad} ...\n");
    int32_t nd = r2_align_count_degrees_of_freedom(ni, arad);
    r2_t U[nd*ni];
    double urad[nd];
    r2_align_compute_search_ellipsoid (ni, arad, nd, U, urad);
    
    /* Validate {U,urad}: */
    fprintf(stderr, "... validating the search ellipsoid basis ...\n");
    r2_t t[ni]; /* A corner of the enclosing box of {\CF} */
    for (int32_t i = 0; i < ni; i++) { t[i] = (r2_t){{0.0, 0.0 }}; }
    for (int32_t k = 0; k < nd; k++)
      { fprintf(stderr, "  checking vector {U[%d]} and radius {urad[%d]} ...\n", k, k);
        
        r2_t *uk = &(U[k*ni]);
        if (debug) { r2_align_print_vector(stderr, ni, "    u", k, uk, FALSE); }
        if (debug) { fprintf(stderr, "    urad[%d] = %.8f\n", k, urad[k]); }
        
        if (k > 0)
          { fprintf(stderr, "    checking decreasing radius order...\n");
            demand(urad[k] <= urad[k-1], "radii out of order");
          }
        
        fprintf(stderr, "    checking if {U[%d]} is conformal with {arad}...\n", k);
        for (int32_t i = 0; i < ni; i++)
          { for (int32_t j = 0; j < 2; j++)
              { double rij = arad[i].c[j];
                if (rij == 0) { demand(uk[i].c[j] == 0, "{uk} is not conformal"); }
              }
          }

        fprintf(stderr, "    checking if {U[%d]} is balanced...\n", k);
        for (int32_t j = 0; j < 2; j++)
          { double sum = 0.0;
            for (int32_t i = 0; i < ni; i++) { sum += uk[i].c[j]; }
            if (debug) { fprintf(stderr, "      sum U[%d][*].c[%d] = %24.16e\n", k, j, sum); }
            demand(fabs(sum) < 1.0e-8, "not balanced");
          }
  
        fprintf(stderr, "    checking if {U[%d]} is normalized...\n", k);
        double sdot = r2_align_dot(ni, uk, uk);
        if (debug) { fprintf(stderr, "      dot(U[%d],U[%d]) = %.8f\n", k, k, sdot); }
        demand(fabs(sdot - 1.0) < 1.0e-8, "not normalized");
       
        if (k > 0) 
          { fprintf(stderr, "    checking if {U[%d]} is orthogonal to {U[0..%d]}...\n", k, k-1);
            for (int32_t r = 0; r < k; r++) 
              { r2_t *ur = &(U[r*ni]);
                double sdot = r2_align_dot(ni, uk, ur);
                if (debug) { fprintf(stderr, "      dot(U[%d],U[%d]) = %.8f\n", r, k, sdot); }
                demand(fabs(sdot) < 1.0e-8, "not orthogonal");
              }
          }
          
        fprintf(stderr, "    checking whether poles of {\\RF} are on boundary of {\\RE} ...\n");
        double sum2 = 0; 
        for (int32_t i = 0; i < ni; i++)
          { for (int32_t j = 0; j < 2; j++)
              { double skij = urad[k]*uk[i].c[j];
                double rij = arad[i].c[j];
                if (rij != 0) { double ekij = skij/rij; sum2 += ekij*ekij; }
                /* Add pole adjustment to diagonal adjustment {t}: */
                t[i].c[j] += skij;
              }
          }
        if (debug) { fprintf(stderr, "      rel dist = %.8f\n", sqrt(sum2)); }
        demand(fabs(sum2 - 1.0) < 1.0e-8, "{urad[k]*uk} not on boundary of {\\RE}");
      }
    fprintf(stderr, "  checking diagonal adjustment ...\n");
    double dr2 = r2_align_rel_dist_sqr (ni, t, NULL, arad);
    demand(fabs(dr2 - nd) < 1.0e-8, "diagonal mismatch");

    fprintf(stderr, "  search ellipsoid OK!\n\n");
  }
    
void ralt_plot_rel_dist_sqr(int32_t ni, r2_t ctr[], r2_t arad[])
  {
    char *fname = NULL;
    asprintf(&fname, "out/f%03d.dat", ni);
    FILE *wr = open_write(fname, TRUE);
    free(fname);

    auto double f2 (int32_t ni, r2_t p[]);
      /* The function to plot. */
    
    double step = 0.25;
    r2_align_plot_mismatch(wr, ni, ctr, arad, step, &f2);
    fclose(wr);

    fprintf(stderr, "done.\n");
    return;
    
    /* Internal implementations: */
    
    double f2 (int32_t ni, r2_t p[])
      { double fval = r2_align_rel_dist_sqr(ni, ctr, p, arad);
        return fval; 
      }

  }

