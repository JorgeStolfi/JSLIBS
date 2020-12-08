#ifndef btc_bubble_nl_opt_adjust_parameters_H
#define btc_bubble_nl_opt_adjust_parameters_H

/* Non-linear optimization of BTC price bubbe parameters. */
/* Last edited on 2015-04-22 20:48:43 by stolfilocal */

#include <btc_bubble_t.h>

void btc_bubble_nl_opt_adjust_parameters
  ( int nd, 
    char* dt[], 
    double ap[], 
    double wt[],
    int nb, 
    btc_bubble_t bp_lo[], /* Min values of parameters, or NIL. */ 
    btc_bubble_t bp[],    /* (IN/OUT) Guessed values of parameters. */    
    btc_bubble_t bp_hi[], /* Max values of parameters, or NIL. */
    int hrad,
    int maxLSQIters, 
    int maxNLIters, 
    int id_ini,
    int id_fin,
    char* outPrefix, 
    double bval[]         /* (OUT) Bubble basis for {bp}. */
  );
  /* Adjusts any adjustable bubble parameters {bp[0..nb-1]} as to get
    the best fit between the modeled price series with those parameters
    and the given price series {ap[0..nd-1]}.
      
    If {maxNLIters} is positive, some parameters of some bubbles (other
    than the coefficients {.coef}) may be adjusted by non-linear
    optimization. The parameters to be adjusted are found by comparing
    the values of certain parameters in bracketing sets {bp_lo} and
    {bp_hi} (see
    {btc_bubble_nl_opt_gather_integer_variable_parameters}).
    
    More precisely, for each parameter that can be adjusted, its value
    {bp[jb].{XX}} in each bubble {jb} must be bracketed by its values
    {bp_lo[jb].{XX}} and {bp_hi[jb].{XX}}. If that interval has nonzero
    width, then that parameter is considered for non-linear
    optimization, within that range.
    
    When a date parameter ({.id_fin_up} or {.id_ini_dn}) needs
    adjusting, its value is stepped through the dates in the specified
    interval. Beware of the combinatorial explosion when trying to
    adjust too many date parameters at the same time. The rates {.rt_up}
    and {.rt_dn} that need fitting are adjusted with {maxNLIters}
    iterations of the edge-divided simplex optimization method.
    
    For any parameter that gets adjusted, the corresponding value
    {bp[jb].{XX}} in {bp} is used as an initial guess (if the parameter
    is continuous) and as the favored solution if there is no clear
    optimum. The best setting of the adjusted parameters, including the
    coefficients, are stroed back into {bp[0..nb-1]}.
   
    Parameters that are never adjusted, such as {.id_ini_sg} and
    {.id_fin_sg}, must be equal in all three parameter sets
    {bp_lo,bp,bp_hi}. The fields {.coef}, {.tag}, and {.color} of
    {bp_lo[jb]} and {bp_hi[jb]}, for all bubbles {jb}, are ignored.
    
    If {maxNLIters} is zero, the parameters in {bp} (other than {.coef})
    are not changed. The parameters {bp_lo,bp_hi} can be NULL in that case. 
    
    In any case, the procedure stores in {bval[0..nd*nb-1]} the the
    bubble function basis corresponding to the parameters {bp[0,,nb-1]},
    smoothed with Hann window with radius {hrad}; see
    {btc_bubble_compute_basis}.  Then it adjusts the coefficients
    {bp[0..nb-1].coef} with {btc_bubble_fit_lsq}, 
    using {maxLSQIters} iterations.  This is done before exiting and 
    internally at every trial set of parameters.
    
    The goal function for non-linear optimization is the RMS error
    bewteen the modeled series and the given series, in log scale,
    between the samples {id_ini} and {id_fin} inclusive. */

#endif
