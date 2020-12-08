/* See {btc_bubble_nl_opt_adjust_parameters.h} */
/* Last edited on 2015-04-22 20:49:43 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsmath.h>

#include <btc_bubble_t.h>
#include <btc_bubble_nl_opt_gather_integer_variable_parameters.h>
#include <btc_bubble_nl_opt_gather_continuous_variable_parameters.h>
#include <btc_bubble_parms_copy.h>
#include <btc_bubble_nl_opt_set_integer_variable_parameters.h>
#include <btc_bubble_nl_opt_add_integer_parameter_bias.h>
#include <btc_bubble_nl_opt_check_bubble_parms_in_range.h>
#include <btc_bubble_nl_opt_adjust_continuous_parameters.h>
#include <btc_bubble_compute_basis.h>
#include <btc_bubble_fit_lsq.h>
#include <btc_bubble_eval_rms_log_error.h>

#include <btc_bubble_nl_opt_adjust_parameters.h>

void btc_bubble_nl_opt_adjust_parameters
  ( int nd, 
    char* dt[], 
    double ap[], 
    double wt[],
    int nb, 
    btc_bubble_t bp_lo[], /* Min values of parameters, or NIL. */
    btc_bubble_t bp[],    /* Guessed values of parameters. */
    btc_bubble_t bp_hi[], /* Max values of parameters, or NIL. */
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
        /* Perform non-linear optimization on the variable parameters in {bp}: */
        demand((bp_lo != NULL) && (bp_hi != NULL), "must specify the parameter bounds");
        
        /* Identify the integer parameters to optimize: */
        int npi = 0; /* Number of integer parameters to optimize. */
        int* pi_lo;  /* Min values of the parameters. */
        int* pi;     /* Guessed values of the parameters. */
        int* pi_hi;  /* Max values of the parameters. */
        btc_bubble_nl_opt_gather_integer_variable_parameters(nb, bp_lo, bp, bp_hi, &npi, &pi_lo, &pi, &pi_hi); 
        fprintf(stderr, "found %d adjustable integer parameters\n", npi);
        
        /* We must call {btc_bubble_nl_opt_adjust_continuous_parameters} even 
          if {npi} is zero, because there may be continuous parameters to adjust. */
        
        /* Keep track of the best parameter set seen: */
        btc_bubble_t bp_best[nb];  /* The best bubble parameter set found so far. */
        int pi_best[npi];          /* Best set of integer parameters found so far. */
        double Q_best = +INF;      /* Best goal function value found so far, or {+INF}. */
        
        /* Step through the integer parameters, remembering the optimum: */
        btc_bubble_t bp_try[nb]; /* Tentative bubble parameter set during NL optimization. */
        int pi_try[npi];      /* The trial parameter vector. */
        int k;
        for (k = 0; k < npi; k++) { pi_try[k] = pi_lo[k]; }
        int n_trials = 0;
        while (TRUE) 
          { 
            /* Try the integer parameter values {pi[0..npi-1]}, resetting the given cont ones: */
            if (verbose)
              { fprintf(stderr, "=== trial %3d ================================================\n", n_trials);
                for (k = 0; k < npi; k++) 
                  { fprintf(stderr, "  pi[%02d] = %5d = %s", k, pi_try[k], dt[pi_try[k]]);
                    fprintf(stderr, "  in %5d = %s", pi_lo[k], dt[pi_lo[k]]);
                    fprintf(stderr, " .. %5d = %s\n", pi_hi[k], dt[pi_hi[k]]);
                  }
              }
            btc_bubble_parms_copy(nb, bp, bp_try);
            btc_bubble_nl_opt_set_integer_variable_parameters(npi, pi_try, nb, bp_lo, bp_try, bp_hi);
            btc_bubble_nl_opt_adjust_continuous_parameters
              ( nd, dt, ap, wt, 
                nb, bp_lo, bp_try, bp_hi,
                hrad, 
                maxNLIters, maxLSQIters, id_ini, id_fin,
                outPrefix, 
                bval 
              );
            n_trials++;
            /* Assumes that {bp_try} is set to the optimum and {bval} is derived from it. */

            /* Compute the goal function for the adjusted continuous params: */
            double Q_try = btc_bubble_eval_rms_log_error(nd, ap, id_ini, id_fin, nb, bp_try, bval);
            if (verbose) { fprintf(stderr, "  Q_try = %25.16e\n", Q_try); }

            /* Add a slight bias to favor the given guess in case of ties or near-ties: */
            double alpha = 1.0e-6; /* Relative bias factor. */
            Q_try = btc_bubble_nl_opt_add_integer_parameter_bias(Q_try, npi, pi_try, pi, alpha);

            if ((Q_try < +INF) && (Q_try < Q_best))
              { /* Save this trial: */
                Q_best = Q_try;
                for (k = 0; k < npi; k++) { pi_best[k] = pi_try[k]; }
                btc_bubble_parms_copy(nb, bp_try, bp_best);
                if (verbose) { fprintf(stderr, "    updated!\n"); }
              }
              
            if (verbose)
              { fprintf(stderr, "==============================================================\n"); }

            /* Increment the integer parameter vector: */
            k = npi - 1;
            while ((k >= 0) && (pi_try[k] == pi_hi[k])) { pi_try[k] = pi_lo[k]; k--; }
            if (k < 0) { break; }
            assert(pi_try[k] < pi_hi[k]);
            pi_try[k]++;
          }
        
        if (verbose)
          { fprintf(stderr, "tried %d integer parameter combinations\n", n_trials);
            fprintf(stderr, "=== optimum ==================================================\n"); 
            for (k = 0; k < npi; k++) 
              { fprintf(stderr, "  pi[%02d] = %5d = %s", k, pi_best[k], dt[pi_best[k]]);
                fprintf(stderr, "  in %5d = %s", pi_lo[k], dt[pi_lo[k]]);
                fprintf(stderr, " .. %5d = %s\n", pi_hi[k], dt[pi_hi[k]]);
              }
            fprintf(stderr, "  Q = %25.16e\n", Q_best);
            fprintf(stderr, "==============================================================\n");
          }
            
        /* Copy the best parameters to {bp}, check range to be sure: */
        assert(Q_best < +INF);
        btc_bubble_parms_copy(nb, bp_best, bp);
        btc_bubble_nl_opt_check_bubble_parms_in_range(nb, bp_lo, bp, bp_hi);

        free(pi_lo); free(pi); free(pi_hi);
      }
      
    /* recompute the basis and fit the linear combination for {bp}: */
    btc_bubble_compute_basis(nd, nb, bp, hrad, bval); 
    btc_bubble_fit_lsq(nd, dt, ap, wt, nb, bp, bval, maxLSQIters, outPrefix);
  }
