/* See {lsq_robust.h} */
/* Last edited on 2024-11-07 00:49:48 by stolfi */

#define lsq_robust_C_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <rmxn.h>
#include <jsmath.h>
#include <gauss_elim.h>
#include <rn.h>

#include <lsq.h>
#include <lsq_array.h>
#include <lsq_robust.h>

void lsq_robust_fit
  ( int32_t nt,
    int32_t nx,
    double X[],
    double F[],
    double W[],
    int32_t maxiter,
    double U[],
    double P[],
    lsq_robust_report_t *report,
    bool_t verbose
  )
  {
    /* bool_t debug = (report != NULL); */
    bool_t debug = FALSE;

    int32_t nf = 1; /* Number of dependent variables (function samples). */

    double *A = rmxn_alloc(nx,nx); /* Moment matrix. */
    double *B = rn_alloc(nx); /* Right-hand side. */

    /* Compute the system matrices {A,B}: */
    lsq_array_compute_matrix(nt, nx, X, W, A);

    /* The matrix does not change: */
    rmxn_inv_full(nx, A, A);

    /* Solve first without correction: */
    lsq_array_compute_rhs(nt, nx, nf, X, F, W, B);
    rmxn_map_col(nx, nx, A, B, U);
    if (report != NULL) { report(0, nt, nx, NULL, F, U); }

    if (maxiter > 0)
      {
        /* Apply the outlier exclusion method: */
        double *Fa = rn_alloc(nt); /* Current approximation. */
        double *Fc = rn_alloc(nt); /* Deviations first, fake values later. */

        /* Estimated statistics: */
        double *Pc = (P != NULL ? P : rn_alloc(nt)); /* Probability of each data point to be inlier. */
        double pri_gud = NAN;  /* A priori probability of being inlier. */
        double avg_gud = NAN;  /* Mean residual of inliers. */
        double dev_gud = NAN;  /* Deviation of inliers. */

        double pri_bad = NAN;  /* A priori probability of being outlier. */
        double avg_bad = NAN;  /* Mean of outliers. */
        double dev_bad = NAN;  /* Deviation of outliers. */

        rn_all(nt, 0.5, Pc); /* A priori, inliers and outliers are equally likely: */
        for (int32_t iter = 1; iter <= maxiter; iter++)
          {
            if (debug) { fprintf(stderr, "    iteration %d\n", iter); }

            /* Set {Fa[0..n-1]} to current approximation: */
            rmxn_map_col(nt, nx, X, U, Fa);

            /* Set {Fc} to the residuals relative to the current approximation: */
            rn_sub(nt, F, Fa, Fc);
            if (debug)
              { fprintf(stderr, "      %3s   %24s  %24s %24s  %24s %24s\n", "it", "F", "Fa", "Fc", "W", "Pc");
                for (int32_t it = 0; it < nt; it++)
                  { fprintf(stderr, "      %3d  %24.16e  %24.16e %24.16e  %24.16e %24.16e\n", it, F[it], Fa[it], Fc[it], W[it], Pc[it]); }
                fprintf(stderr, "\n");
              }

            /* Recompute parameters of inlier residuals: */
            lsq_robust_compute_stats(nt, Fc, W, Pc, TRUE,  0.0,     &avg_gud, &dev_gud, &pri_gud, debug);

            /* Recompute parameters of outlier function samples: */
            lsq_robust_compute_stats(nt, F,  W, Pc, FALSE, dev_gud, &avg_bad, &dev_bad, &pri_bad, debug);

            assert(fabs(pri_bad + pri_gud - 1.0) < 0.0001);

            /* Recompute inlier/outlier probabilities {Pc[0..nt-1]} and adjusted data {Fc[0..nt-1]}: */
            for (int32_t k = 0; k < nt; k++)
              { /* Decide the probability {Pc[k]} of each data point {X[k,*},F[k,*]} being an inlier: */
                /* Grab the function sample {Fk} and current approximation {Fak}: */
                double Fk = F[k];
                double Fak = Fa[k];
                /* Compute its assumed average if inlier: */
                double avg_gud_k = Fak; /* Assumed average of inlier residual distribution at this point. */
                /* Compute the probability of {Fk} being inlier: */
                double Pk = lsq_robust_bayes(Fk, pri_gud, avg_gud_k, dev_gud, avg_bad, dev_bad, verbose);
                assert(isfinite(Pk) && (Pk >= 0.0) && (Pk <= 1.0));
                Pc[k] = Pk;
                /* Set {Fc[k]} to the expected value of the inlier part of {F[k]}: */
                Fc[k] = Pk*Fk + (1-Pk)*Fak;
              }

            /* Recompute coeffs using {Fc} instead of {F}: */
            lsq_array_compute_rhs(nt, nx, nf, X, Fc, W, B);
            rmxn_map_col(nx, nx, A, B, U);
            if (report != NULL) { report(iter, nt, nx, Pc, Fc, U); }
          }
        free(Fa);
        free(Fc);
        if (Pc != P) { free(Pc); }
        free(A);
        free(B);
      }
  }

void lsq_robust_compute_stats
  ( int32_t nt,
    double Y[],   /* Values to analyze (residuals or function samples, {nt} elements). */
    double W[],
    double P[],
    bool_t inlier,
    double dev_ref,
    double *avgP,
    double *devP,
    double *priP,
    bool_t verbose
  )
  {
    /* Compute the average and overall probability: */
    double sum_wpy = 0;
    double sum_wp = 0;
    double sum_w = 0;
    for (int32_t k = 0; k < nt; k++)
      { double Yk = Y[k];
        double Wk = (W == NULL ? 1.0 : W[k]);
        double Pk = (P == NULL ? 0.5 : (inlier ? P[k] : 1.0 - P[k]));
        sum_wpy += Wk*Pk*Yk;
        sum_wp += Wk*Pk;
        sum_w += Wk;
      }
    if (verbose) { fprintf(stderr, "    sum_wpy = %24.16e  sum_wp = %24.16e  sum_w = %24.16e\n", sum_wpy, sum_wp, sum_w); }
    sum_wpy += 0.5e-300; /* To avoid division by zero. */
    sum_wp += 0.5e-300; /* To avoid division by zero. */
    sum_w += 1.0e-300; /* To avoid division by zero. */
    double avg = sum_wpy/sum_wp; assert(isfinite(avg));
    double pri = sum_wp/sum_w; assert(isfinite(pri));

    /* Compute the variance and deviation: */
    double sum_wpd2 = 0;
    for (int32_t k = 0; k < nt; k++)
      { double Dk = Y[k] - (inlier ? 0.0 : avg);
        double Wk = (W == NULL ? 1.0 : W[k]);
        double Pk = (P == NULL ? 0.5 : (inlier ? P[k] : 1.0 - P[k]));
        sum_wpd2 += Wk*Pk*Dk*Dk;
      }
    if (verbose) { fprintf(stderr, "    sum_wpd2 = %24.16e\n", sum_wpd2); }
    double var = sum_wpd2/sum_wp; assert(isfinite(var));
    var += 1.0e-300; /* Make sure that the variance is nonzero. */
    assert(isfinite(var) && (var > 0));
     
    if (inlier)
      { /* Scale up the variance to compensate the shrink effect: */
        double alpha = 1.1; /* !!! Figure out the right factor !!! */
        var = alpha * var;
      }
    else
      { double beta = 2.0; /* Variance broadening factor. !!! Figure out !!! */
        double var_min = beta * dev_ref * dev_ref;
        if (var < var_min) { var = var_min; }
      }

    double dev = sqrt(var);

    if (verbose)
      { char* vname = ((char*[2]){ "function samples", "residuals" })[inlier];
        char* tag = ((char*[2]){ "outlier", "inlier" })[inlier];
        lsq_robust_debug_distr(stderr, vname, tag, avg, dev, pri);
      }
    /* Return: */
    assert(isfinite(avg));
    assert(isfinite(dev) && (dev >= 0));
    assert(isfinite(pri) && (pri >= 0) && (pri <= 1.0));
    (*avgP) = avg;
    (*devP) = dev;
    (*priP) = pri;
  }

double lsq_robust_bayes
  ( double F,
    double pri_gud,
    double avg_gud,
    double dev_gud,
    double avg_bad,
    double dev_bad,
    bool_t verbose
  )
  {
    /* Degenrate cases: */
    if ((dev_gud < 1.0e-200) || (dev_bad < 1.0e-200))
      { if (dev_gud < dev_bad) 
          { return 0.0; }
        else if (dev_gud > dev_bad) 
          { return 1.0; }
        else
          { return pri_gud; }
      }
    /* Compute the normalized discrepancies relative to the inlier and outlier means: */
    double d_gud = (F - avg_gud)/dev_gud; /* Normalized displacement of {F} from inlier average. */
    if (d_gud > 6.0) { return 0.0; }

    double d_bad = (F - avg_bad)/dev_bad; /* Normalized displacement of {F} from outlier average. */
    if (d_bad > 6.0) { return 1.0; }

    /* Apply Bayes's formula: */
    double pri_bad = 1 - pri_gud;
    double Pc_gud = exp(-d_gud*d_gud/2)/dev_gud; /* Prob of {F} assuming inlier, times {sqrt(2*PI)}. */
    double Pc_bad = exp(-d_bad*d_bad/2)/dev_bad; /* Prob of {F} assuming outlier, times {sqrt(2*PI)}. */
    double PP_gud = pri_gud*Pc_gud;
    double PP_bad = pri_bad*Pc_bad;
    return PP_gud/(PP_gud + PP_bad);
  }

void lsq_robust_debug_distr(FILE *wr, char *vname, char *tag, double avg, double dev, double pri)
  {
    fprintf(wr, "--- distribution of %s for %s data points ---\n", vname, tag);
    fprintf(wr, "mean = %+18.10e", avg);
    fprintf(wr, "  deviation = %17.10e", dev);
    fprintf(wr, "  prior prob = %9.7f\n", pri);
  }
