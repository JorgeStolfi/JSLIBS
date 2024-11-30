#ifndef btc_bubble_nl_opt_do_adjust_continuous_parameters_H
#define btc_bubble_nl_opt_do_adjust_continuous_parameters_H

/* Non-linear optimization of the continuous parameters of a BTC price bubble. */
/* Last edited on 2024-11-08 18:17:14 by stolfi */

#include <btc_bubble_t.h>

void btc_bubble_nl_opt_do_adjust_continuous_parameters
  ( int npf,         /* Number of continuous parameters to optimize. */
    double pf_lo[],  /* Min values of the parameters. */
    double pf[],     /* Guessed values of the parameters. */
    double pf_hi[]   /* Max values of the parameters. */
  );
  /* Like {btc_bubble_nl_opt_adjust_continuous_parameters}, but when {maxLSQIters}
    and {npf} are known and positive,
    and the lower limit {pf_lo[0..npf-1]}, initial guess {pf[0..npf-1]},
    and upper limit {pf_hi[0..npf-1]} of the continuously adjustable parameters
    have been determined. */
    
#endif
