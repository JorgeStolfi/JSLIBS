#define PROG_NAME "test_hr2_pmap_generic_opt"
#define PROG_DESC "test of {hr2_pmap_generic_opt.h} and related modules"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-08 11:17:15 by stolfi */ 
/* Created on 2020-07-11 by J. Stolfi, UNICAMP */

#define test_hr2_pmap_generic_opt_COPYRIGHT \
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
#include <rmxn.h>
#include <rmxn_extra.h>
#include <bool.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>
#include <argparser.h>

#include <hr2_pmap_from_many_pairs.h>
#include <hr2_pmap_encode.h>
#include <hr2_pmap_opt.h>
#include <hr2_pmap_test_tools.h>
#include <hr2_pmap_opt_test_tools.h>

#include <hr2_pmap_generic_opt.h>

#define POINT_PERT 1.0e-6
  /* Absolute perturbation to apply to the data points. */
    
#define ht_GENERIC hr2_pmap_type_GENERIC
  /* Shorter name. */

/* The program takes no command line options. */

void hgot_test_hr2_pmap_generic_opt_quadratic__hr2_pmap_generic_opt_1D_plot
  ( sign_t sgn,
    bool_t verbose,
    char *outPrefix
  );
  /* Tests {hr2_pmap_generic_opt_quadratic}, and, if {outPrefix} is not
    {NULL}, also {hr2_pmap_generic_opt_1D_plot}, creating files with
    prefix {outPrefix}. */

int32_t main(int32_t argn, char **argv);

void hgot_do_test_opt_and_plot
  ( sign_t sgn,
    bool_t tight,
    bool_t ident,
    bool_t close,
    bool_t verbose,
    char *outPrefix
  );
  /* Tests {hr2_pmap_generic_opt_quadratic} for maps of the given sign.
     If {outPrefix} is not {NULL}, also tests {hr2_pmap_generic_opt_1D_plot}.
  
    The goal function for a map {M} is the discrepancy squared between a
    certain number {np} of random point pairs, mapped forward and
    backwards by {M}. The points are chosen so that the optimum map has
    the specified {sgn}.
    
    If {tight}, tests with the minimum number of points needed
    to determine the map, namely 4.
    
    If {ident}, the points are exactly related by the identity
    transformation (that is, each point is to be mapped to itself) if
    {sgn} is {+1} or by a Y coordinate negation if {sgn} is {-1};
    minus a small amount of noise.
    
    if {close} is true, puts the initial guess of the map very close to the
    optimal map, else starts with a random map.
    
    Also writes 1D plot data of the goal function for the initial and final maps.
    The tile names will be "{outPrefix}-{tag}-{stage}.txt" where
    {tag} will be "{xsgn}-tight{tight}-ident{ident}-close{close}",
    {xsgn} will be "m", "o", or "p" for {-1}, 0, and {+1},respectively, and 
    and {stage} will be "ini" or "opt". */

void hgot_do_test_1D_plot
  ( char *outPrefix,
    char *tag,
    char *stage,
    sign_t sgn,
    hr2_pmap_t *M,
    hr2_pmap_opt_func_t *f2,
    double urad,
    bool_t verbose
  );
  /* Plots the goal function {f2} in the neighborhood of the projective map {M},
    which must be of the given {sgn}.
    
    The plots will extend {urad} away from {M} in encoding parameter space.
    The tile names will be "{outPrefix}-{tag}-{stage}.txt".  */
  
hr2_pmap_t hgot_perturb_map(hr2_pmap_t *M, double pabs, double prel, sign_t sgn);
  /* Returns a slightly perturbed copy of map {M},
    using 

    The map {M} is assumed to be of the given {type} an {sgn},
    and the result will be too. */

void hgot_check_opt_pmap
  ( hr2_pmap_t *M,
    sign_t sgn,
    double epsy,
    hr2_pmap_opt_func_t *f2,
    double f2Stop,
    bool_t verbose
  );
  /* If {f2(M)} does not exceed {f2Stop}, does nothing.
    Else evaluates {f2} at various perturbed versions {N} of {M},
    and fails if {M} does not look optimum at all -- specifically,
    if {f2(N) < f2(M) - 0.5*f2Stop}. The norm of the perturbations
    will be {epsy}.  */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srand(4615*417);
    
    for (sign_t sgn = +1; sgn >= -1; sgn -= 2)
      { for (int32_t trial = 0; trial < 10; trial++)
          { bool_t verbose = (trial < 2);
            char *outPrefix = NULL;
            if (trial < 2) { asprintf(&outPrefix, "out/test-%03d", trial); }
            hgot_test_hr2_pmap_generic_opt_quadratic__hr2_pmap_generic_opt_1D_plot(sgn, verbose, outPrefix);
            free(outPrefix);
          }
      }
    return 0;
  }
  
void hgot_test_hr2_pmap_generic_opt_quadratic__hr2_pmap_generic_opt_1D_plot(sign_t sgn, bool_t verbose, char *outPrefix)
  {
    fprintf(stderr, "> --- %s sgn = %+2d ---\n", __FUNCTION__, sgn);
    
    for (int32_t kident = 1; kident >= 0; kident--)
      { bool_t ident = (kident == 1);
        for (int32_t ktight = 1; ktight >= 0; ktight--)
          { bool_t tight = (ktight == 1);
            for (int32_t kclose = 1; kclose >= 0; kclose--)
              { bool_t close = (kclose == 1);
                hgot_do_test_opt_and_plot
                  ( sgn, 
                    tight, ident, close, 
                    verbose, outPrefix
                  );
              }
          }
      }

    fprintf(stderr, "< --- %s ---\n", __FUNCTION__);
  }
  
void hgot_do_test_opt_and_plot
  ( sign_t sgn,
    bool_t tight,
    bool_t ident,
    bool_t close,
    bool_t verbose,
    char *outPrefix
  )
  {
    if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
    if (verbose)
      { fprintf(stderr, "  > --- %s sgn = %+d", __FUNCTION__, sgn);
        fprintf(stderr, " tight = %c ident = %c", "FT"[tight], "FT"[ident]);
        fprintf(stderr, " close = %c --- \n", "FT"[close]);
      }
    
    char *xsgn = (sgn == 0 ? "o" : (sgn < 0 ? "m" : "p"));
    char *tag = NULL;
    asprintf
      ( &tag, "%s-tight%c-ident%c-close%c", 
        xsgn, "FT"[tight], "FT"[ident], "FT"[close]
      );
    
    double pert = POINT_PERT;
    int32_t np;
    r2_t *p1 = NULL;
    r2_t *p2 = NULL;
    double *w = NULL;
    hr2_pmap_t MP;
    hr2_pmap_opt_test_choose_r2_point_pairs(ht_GENERIC, sgn, tight, ident, pert, verbose, &np, &p1, &p2, &w, &MP);
    if (verbose) 
      { fprintf(stderr, "    testing with %d point pairs\n", np);
        double f2MP = hr2_pmap_mismatch_sqr(&MP, np, p1, p2, w);
        hr2_pmap_opt_test_print_map("nominal p1->p2 map", &MP, f2MP);
      }
    
    hr2_pmap_t M;  /* Initial guess and result matrix. */
    if (close)
      { /* Start close to the correct map: */
        double pabs = 1.0e-8;
        double prel = 0.0001;
        M = hgot_perturb_map(&MP, pabs, prel, sgn);
      }
    else
      { /* Start with a random map: */
        hr2_test_throw_pmap(&M);
        hr2_pmap_set_sign(&M, sgn);
      }
      
    /* Compute max point-to-point error of the initial guess: */
    double maxDist = hr2_pmap_max_mismatch(&M, np, p1, p2); 
    if (verbose) { fprintf(stderr, "    max p1<-->p2 maping dist = %20.14f\n", maxDist); }
    
    auto double goal_func(hr2_pmap_t *P);
      /* Goal function: point mismatch squared. */
        
    auto bool_t ok_pred(hr2_pmap_t *P, double fP);
      /* Checks whether {P} is acceptable. Currently, that
        means {fP} is less than {f2Stop}. */
    
    int32_t maxIter = 20;
    double maxMod = 5.0*fmax(1.0, maxDist);
    double f2Stop = POINT_PERT*POINT_PERT;
    if (verbose) { fprintf(stderr, "    maxIter = %d maxMod = %20.14f  f2Stop = %20.14f = %12.4e\n", maxIter, maxMod, f2Stop, f2Stop); }
    
    double f2M = goal_func(&M);
    if (verbose) { hr2_pmap_opt_test_print_map("initial", &M, f2M); }
    if (outPrefix != NULL)
      { hgot_do_test_1D_plot(outPrefix, tag, "ini", sgn, &M, goal_func, maxMod, verbose); }
    
    hr2_pmap_generic_opt_quadratic(sgn, &goal_func, &ok_pred, maxIter, maxMod, &M, &f2M, verbose);
    
    if (verbose) { hr2_pmap_opt_test_print_map("final", &M, f2M); }
    if (outPrefix != NULL)
      { hgot_do_test_1D_plot(outPrefix, tag, "opt", sgn, &M, goal_func, 0.01*maxMod, verbose); }

    double epsy = 0.001;        /* Size of perturbations to encoded map. */
    hgot_check_opt_pmap(&M, sgn, epsy, &goal_func, f2Stop, verbose);

    free(tag);

    if (verbose) { fprintf(stderr, "  < --- %s ---\n", __FUNCTION__); }
    if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
    
    return;
    
    double goal_func(hr2_pmap_t *P)
      { double mm2 = hr2_pmap_mismatch_sqr(P, np, p1, p2, w);
        return mm2;
      }
     
    bool_t ok_pred(hr2_pmap_t *P, double fP)
      { if (verbose) { hr2_pmap_show_point_mismatch(P, np, p1, p2, w); }
        if (verbose) 
          { /* Check with the nominal map: */
            double d = sqrt(hr2_pmap_diff_sqr(P, &MP));
            fprintf(stderr, "    distance from nominal map = %20.14f = %12.4e\n", d, d);
          }
        return fP < f2Stop;
      }
  }

hr2_pmap_t hgot_perturb_map(hr2_pmap_t *M, double pabs, double prel, sign_t sgn)
  { 
    hr2_pmap_t N = (*M);
    rmxn_perturb_unif(3, 3, pabs, prel, &(N.dir.c[0][0]));
    r3x3_inv(&(N.dir), &(N.inv));
    hr2_pmap_set_sign(&N, sgn);
    return N;
  }

void hgot_do_test_1D_plot
  ( char *outPrefix,
    char *tag,
    char *stage,
    sign_t sgn,
    hr2_pmap_t *M,
    hr2_pmap_opt_func_t *f2,
    double urad,
    bool_t verbose
  ) 
  { 
    if (verbose) { fprintf(stderr, ". . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n"); }
    if (verbose) { fprintf(stderr, "    > --- %s ---\n", __FUNCTION__); }
    int32_t nu = 7;
    int32_t ns = 50;
    hr2_pmap_generic_opt_1D_plot(outPrefix, tag, stage, M, sgn, f2, nu, urad, ns);

    if (verbose) { fprintf(stderr, "    < --- %s ---\n", __FUNCTION__); }
    if (verbose) { fprintf(stderr, ". . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n"); }
  }

void hgot_check_opt_pmap
  ( hr2_pmap_t *M,
    sign_t sgn,
    double epsy,
    hr2_pmap_opt_func_t *f2,
    double f2Stop,
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "    > --- %s ---\n", __FUNCTION__); }
    if (verbose) { fprintf(stderr, "    f2Stop = %.18f  epsy = %.18f\n", f2Stop, epsy); }

    double f2M = f2(M);
    if (verbose) { fprintf(stderr, "    f2(M) = %.12f = %24.16e\n", f2M, f2M); }
    if (f2M <= f2Stop) 
      { /* Good enough: */ 
        if (verbose) { fprintf(stderr, "    M is good enough!\n"); }
      }
    else
      { if (verbose)
          { fprintf(stderr, "    !! M is not good enough: ");
            fprintf(stderr, "  f2M = %20.14f = %12.4e\n", f2M, f2M);
            fprintf(stderr, "  f2Stop = %20.14f = %12.4e\n", f2Stop, f2Stop);
            fprintf(stderr, "    probing maps around M\n");
          }

        double f2DropLim = 0.75*f2Stop;   /* Neighboring matrices cannot have {goal_func} this much lower than {M}. */
        bool_t ok = TRUE;
        hr2_pmap_t N_best;
        double f2N_best = f2M;
        int32_t nd = 30;
        double pabs = epsy;
        double prel = 0.001;
        for (int32_t kd = 0; kd < nd; kd++)
          { hr2_pmap_t N = hgot_perturb_map(M, pabs, prel, sgn);
            double f2N = f2(&N);
            if (f2M - f2N > f2DropLim)
              { if (verbose) 
                  { fprintf(stderr, "    oops, found a significantly better matrix:\n");
                    hr2_pmap_opt_test_print_map("N", &N, f2N);
                  }
                ok = FALSE;
                if (f2N < f2N_best) { f2N_best = f2N; N_best = N; }
              }
          }
        if (! ok) 
          { fprintf(stderr, "    !! M is definitely NOT OPTIMAL\n");
            hr2_pmap_opt_test_print_map("N_best", &N_best, f2N_best);
          }
        else
          { if (verbose) { fprintf(stderr, "    M seems to be optimal or nearly so\n"); } }
      }
   if (verbose) { fprintf(stderr, "    < --- %s ---\n", __FUNCTION__); }
  }
