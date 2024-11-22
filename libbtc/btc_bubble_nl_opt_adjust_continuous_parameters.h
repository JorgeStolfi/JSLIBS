#ifndef btc_bubble_nl_opt_adjust_continuous_parameters_H
#define btc_bubble_nl_opt_adjust_continuous_parameters_H

/* Non-linear optimization of the continuous parameters of a BTC price bubble. */
/* Last edited on 2024-11-08 18:18:52 by stolfi */

#include <btc_bubble_t.h>

void btc_bubble_nl_opt_adjust_continuous_parameters
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
    double bval[]         /* (OUT) Bubble basis computed with parameters {bp}. */
  );
  /* Adjusts any adjustable continuous bubble parameters in {pb[0..nb-1]}
    so as to best fit the price series {ap[0..nd-1]}. 
     
    If {maxNLIters} is positive, some parameters of some bubbles (other
    than the coefficients {.coef}) may also be adjusted by non-linear
    optimization, specifically with {maxNLIters} iterations of the
    edge-divided simplex optimization method. The parameters to be
    adjusted are found by comparing the values of certain parameters in
    bracketing sets {bp_lo} and {bp_hi} (see
    {btc_bubble_nl_opt_gather_continuous_variable_parameters}).
    
    More precisely, for each parameter that can be adjusted, its value
    {bp[jb].{XX}} in each bubble {jb} must be bracketed by its values
    {bp_lo[jb].{XX}} and {bp_hi[jb].{XX}}. If that interval has nonzero
    width, then that parameter is considered for non-linear optimization
    
    For any parameter that gets adjusted, the corresponding value
    {bp[jb].{XX}} in {bp} is used as an initial guess, and as the
    favored solution if there is no clear optimum. The best setting of
    the adjusted parameters, including the coefficients, are stored back
    into {bp[0..nb-1]}.
    
    Parameters that are never adjusted, such as {.id_ini_sg} and
    {.id_fin_sg}, must have {bp_lo,bp,bp_hi} all equal. For adjstable
    integer parameters, like {.id_fin_up} or {.id_ini_dn}, their values
    in {bp} must be bracketed by the correspoding values in
    {bp_lo,bp_hi}, and are not changed by the procedure. The fields
    {.coef}, {.tag}, and {.color} of {bp_lo[jb]} and {bp_hi[jb]}, for
    all bubbles {jb}, are ignored.
   
    If {maxNLIters} is zero, the parameters in {bp} (other than {.coef})
    are not changed. The parameters {bp_lo,bp_hi} can be NULL in that case. 
    
    In any case, the procedure stores in {bval[0..nd*nb-1]} the the
    bubble function basis, smoothed with Hann window with radius {hrad};
    see {btc_bubble_compute_basis}. Then it adjusts the coefficients
    {bp[0..nb-1].coef} with robust least-squares fitting (see
    {btc_bubble_fit_lsq}), with given prices {ap[0..nd-1]} and weights
    {wt[0..nd-1]}, using {maxLSQIters} iterations. This is done before
    exiting and internally at every trial set of parameters.
    
    The goal function for non-linear optimization is the RMS error
    bewteen the modeled series and the given series, in log scale,
    between the samples {id_ini} and {id_fin} inclusive. */
    

#endif
