/* See {btc_bubble_nl_opt_adjust_continuous_parameters.h} */
/* Last edited on 2025-03-19 12:46:24 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <rn.h>
#include <jsmath.h>
#include <sve_minn.h>
#include <sve_minn_iterate.h>

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
#include <btc_bubble_nl_opt_do_adjust_continuous_parameters.h>

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
          { btc_bubble_nl_opt_do_adjust_continuous_parameters
              ( npf, pf_lo, pf, pf_hi );
          }
        free(pf_lo); free(pf); free(pf_hi);
      }
            
    /* Compute the basis for {bp} and adjust the linear coefficients {.coef}: */
    btc_bubble_compute_basis(nd, nb, bp, hrad, bval); 
    btc_bubble_fit_lsq(nd, dt, ap, wt, nb, bp, bval, maxLSQIters, outPrefix);
  }
