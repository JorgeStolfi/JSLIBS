/* See {btc_bubble_nl_opt_adjust_continuous_parameters.h} */
/* Last edited on 2017-03-13 21:06:32 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <rn.h>
#include <jsmath.h>
#include <sve_minn.h>

#include <btc_bubble_t.h>
#include <btc_bubble_nl_opt_gather_continuous_variable_parameters.h>
#include <btc_bubble_nl_opt_map_continuous_parameters.h>
#include <btc_bubble_nl_opt_unmap_continuous_parameters.h>
#include <btc_bubble_nl_opt_compute_dmax.h>
#include <btc_bubble_nl_opt_set_continuous_variable_parameters.h>
#include <btc_bubble_nl_opt_check_bubble_parms_in_range.h>
#include <btc_bubble_compute_basis.h>
#include <btc_bubble_fit_lsq.h>
#include <btc_bubble_parms_copy.h>
#include <btc_bubble_eval_rms_log_error.h>

#include <btc_bubble_nl_opt_adjust_continuous_parameters.h>

void btc_bubble_nl_opt_adjust_continuous_parameters
  ( int nd, 
    char* dt[], 
    double ap[], 
    double wt[],
    int nb, 
    btc_bubble_t bp_lo[], 
    btc_bubble_t bp[],    
    btc_bubble_t bp_hi[], 
    int hrad,
    int maxLSQIters, 
    int maxNLIters, 
    int id_ini,
    int id_fin,
    char* outPrefix, 
    double bval[]
  )
  {
    bool_t verbose = TRUE;
    
    if (maxNLIters > 0) 
      { 
        /* Identify the continuous parameters to optimize: */
        int npf = 0; /* Number of continuous parameters to optimize. */
        double* pf_lo;  /* Min values of the parameters. */
        double* pf;     /* Guessed values of the parameters. */
        double* pf_hi;  /* Max values of the parameters. */
        btc_bubble_nl_opt_gather_continuous_variable_parameters(nb, bp_lo, bp, bp_hi, &npf, &pf_lo, &pf, &pf_hi); 
        fprintf(stderr, "found %d adjustable continuous parameters\n", npf);

        if (npf > 0)
          { 
            /* For the non-linear optimization, we change the scale of all 
              parameters so that {pf[k]} in the interval {pf_lo[k] _ pf_hi[k]} is mapped 
              to {x[k]} in {[-1 _ 1]}. */

            /* Save the {x}-vector of input values  of initial guess, for bias. */
            double x_ini[npf]; /* Initial {x}-vector. */
            btc_bubble_nl_opt_map_continuous_parameters(npf, pf_lo, pf, pf_hi, x_ini);

            int n_evals = 0; /* Number of calls to {eval_Q}. */
            
            auto double eval_Q(int n, double x[]);
              /* Goal function for optimization.  Basically the RMS error in the 
                natural logs of model and actual prices.  The model prices
                are computed by unmapping {x} to obtain parameter values {pf},
                clipping them to the box {[pf_lo _ pf_hi]}
                then inserting {pf} into {bp}, building the table {bval},
                computing the coefficients {bp[].coef} by robust least squares,
                and adding a small bias term and an off-limits penalty. */

            double eval_Q(int n, double x[])
              {
                assert(n == npf);
                /* Unmap the {x} arguments and store them in {bp}, get distance {dBox} from box: */
                double dBox2 = btc_bubble_nl_opt_unmap_continuous_parameters(npf, x, pf_lo, pf, pf_hi);
                if (verbose)
                  { fprintf(stderr, "  --- eval %5d ------------------------------------------------\n", n_evals);
                    int k;
                    for (k = 0; k < npf; k++) { fprintf(stderr, "    pf[%02d] = %25.16e\n", k, pf[k]); }
                  }
                btc_bubble_nl_opt_set_continuous_variable_parameters(npf, pf, nb, bp_lo, bp, bp_hi);
                btc_bubble_nl_opt_check_bubble_parms_in_range(nb, bp_lo, bp, bp_hi);

                /* Compute the bubble basis and fit the coefficients {bp[].coef}: */
                btc_bubble_compute_basis(nd, nb, bp, hrad, bval); 
                btc_bubble_fit_lsq(nd, dt, ap, wt, nb, bp, bval, maxLSQIters, outPrefix); 

                /* Compute the RMS log error and add bias and penalty terms: */
                double Q = btc_bubble_eval_rms_log_error(nd, ap, id_ini, id_fin, nb, bp, bval);
                
                if (dBox2 > 0) 
                  { /* Add out-of-box penalty: */
                    double beta = 1.000; 
                    double penalty = beta*dBox2;
                    if (verbose) 
                      { fprintf(stderr, "    dBox = %25.16e  dBox2 = %25.16e  penalty = %25.16e\n", sqrt(dBox2), dBox2, penalty); }
                    Q += penalty;
                  }
                
                /* Add bias for initial guess: */
                double dIni2 = rn_dist_sqr(npf, x, x_ini);
                if (dIni2 > 0.0)
                  { double gamma = 1.0e-6;
                    double penalty = gamma*dIni2;
                    if (verbose) 
                      { fprintf(stderr, "    dIni = %25.16e  dIni2 = %25.16e  penalty = %25.16e\n", sqrt(dIni2), dIni2, penalty); }
                    Q += penalty;
                  }

                if (verbose) 
                  { fprintf(stderr, "    Q = %25.16e\n", Q);
                    fprintf(stderr, "  ----------------------------------------------------------------\n");
                  }
                
                n_evals++;
                return Q;
              }

            auto bool_t examine(int n, double x[], double Fx);
              /* Predicate that examines a candidate solution and 
                tells whether it is good enough. Then current version 
                always returns {FALSE}. */ 

            bool_t examine(int n, double x[], double Fx)
              { 
                return FALSE;
              }

            /* Compute the maximum search radius: */
            double dMax = btc_bubble_nl_opt_compute_dmax(npf, x_ini);
            bool_t dBox = FALSE;  /* Search in ball, not box. */

            /* Non-linear adjustment of the continuous parameters: */
            double x_opt[npf]; /* Trial argument vector. */
            int ip;
            for (ip = 0; ip < npf; ip++) { x_opt[ip] = x_ini[ip]; }
            double Q_opt = eval_Q(npf, x_opt);

            double rIni = 0.050;   /* Initial simplex radius. */
            double rMin = 0.005;   /* Minimum simplex radius. */
            double rMax = 0.200;   /* Maximum simplex radius. */
            double stop = 1.0e-4;  /* Stop if step is less than this. */

            sve_minn_iterate
              ( npf,        /* int n,          */
                eval_Q,     /* sve_goal_t *F,  */
                examine,    /* sve_pred_t *OK, */
                x_opt,      /* double x[],     */
                &Q_opt,     /* double *FxP,    */
                -1,         /* sign_t dir,     */
                dMax,       /* double dMax,    */
                dBox,       /* bool_t dBox,    */
                rIni,       /* double rIni,    */
                rMin,       /* double rMin,    */
                rMax,       /* double rMax,    */
                stop,       /* double stop,    */
                maxNLIters, /* int maxIters,   */
                verbose     /* bool_t debug    */
              );

            fprintf(stderr, "  tried %d continuous parameter combinations\n", n_evals);

            /* Store the optimum args {x_opt} back into {bp}: */
            assert(Q_opt < +INF);
            (void)btc_bubble_nl_opt_unmap_continuous_parameters(npf, x_opt, pf_lo, pf, pf_hi);
            if (verbose)
              { fprintf(stderr, "  --- optimum ----------------------------------------------------\n");
                int k;
                for (k = 0; k < npf; k++) 
                  { fprintf(stderr, "    pf[%02d] = %25.16e", k, pf[k]);
                    fprintf(stderr, "  in %25.16e", pf_lo[k]);
                    fprintf(stderr, " __ %25.16e\n", pf_hi[k]);
                  }
                fprintf(stderr, "    Q = %25.16e\n", Q_opt);
                fprintf(stderr, "  ----------------------------------------------------------------\n");
              }
                
            btc_bubble_nl_opt_set_continuous_variable_parameters(npf, pf, nb, bp_lo, bp, bp_hi);
            btc_bubble_nl_opt_check_bubble_parms_in_range(nb, bp_lo, bp, bp_hi);
          }

        free(pf_lo); free(pf); free(pf_hi);
      }
            
    /* Compute the basis for {bp} and adjust the linear coefficients {.coef}: */
    btc_bubble_compute_basis(nd, nb, bp, hrad, bval); 
    btc_bubble_fit_lsq(nd, dt, ap, wt, nb, bp, bval, maxLSQIters, outPrefix);
  }
