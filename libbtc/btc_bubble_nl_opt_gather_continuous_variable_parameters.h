#ifndef btc_bubble_nl_opt_gather_continuous_variable_parameters_H
#define btc_bubble_nl_opt_gather_continuous_variable_parameters_H

/* Collecting adjustable integer paramters of a BTC price bubble model. */
/* Last edited on 2015-04-20 13:40:36 by stolfilocal */

#include <btc_bubble_t.h>

void btc_bubble_nl_opt_gather_continuous_variable_parameters
  ( int nb, 
    btc_bubble_t bp_lo[], 
    btc_bubble_t bp[], 
    btc_bubble_t bp_hi[], 
    int* npfP, 
    double** pf_loP, 
    double** pfP, 
    double** pf_hiP
  );
  /* Scans the bubble parameter sets {{bp_lo,bp,bp_hi}[0..nb-1]} for
    continuous ({double}) parameters that are to be adjusted, and returns their values
    in the vectors {{pf_lo,pf,pf_hi}[0..npf-1]}.
    
    More precisely, if {bp_lo[ib].{XX} < bp_hi[ib].{XX}} for some bubble
    index {ib} in {0..nb-1} and some adjustable continuous field
    {.{XX}}, then it considers that parameter of that bubble adjustable,
    and saves those limit values in {pf_lo[ip],pf_hi[ip]}, for
    successive values of {ip}. Also saves the value of {bp[ib].{XX}}
    (which must always be in the specified range) into {pf[ip]}.
    
    Currently considers only the parameters {.rt_up} and {.rt_dn}, in
    that order. The parameters {.id_fin_up}, {.id_ini_dn}, {.id_ini_sg},
    and {.id_fin_sg} in {bp} must be bracketed by the values in
    {bp_lo,bp_hi}, but are otherwise ignored.
    
    The arrays {pf_lo,pf,pf_hi} are allocated by the procedure, and
    their addresses are returned in {*pf_loP,*pfP,*pf_hiP}. The number
    {npf} of adjustable parameters found is returned in {*npfP}. */
    

#endif
