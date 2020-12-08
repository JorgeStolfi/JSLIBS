#ifndef btc_bubble_nl_opt_gather_integer_variable_parameters_H
#define btc_bubble_nl_opt_gather_integer_variable_parameters_H

/* Collecting adjustable integer paramters of a BTC price bubble model. */
/* Last edited on 2015-04-20 13:38:18 by stolfilocal */

#include <btc_bubble_t.h>

void btc_bubble_nl_opt_gather_integer_variable_parameters
  ( int nb, 
    btc_bubble_t bp_lo[], 
    btc_bubble_t bp[], 
    btc_bubble_t bp_hi[], 
    int* npiP, 
    int** pi_loP, 
    int** piP, 
    int** pi_hiP
  );
  /* Scans the bubble parameter sets {{bp_lo,bp,bp_hi}[0..nb-1]} for
    integer parameters that are to be adjusted, and returns their values
    as integer vectors {{pi_lo,pi,pi_hi}[0..npi-1]}.
    
    More precisely, if {bp_lo[ib].{XX} < bp_hi[ib].{XX}} for some bubble
    index {ib} in {0..nb-1} and some adjustable integer field {.{XX}},
    then it considers that parameter of that bubble adjustable, and
    saves those limit values in {pi_lo[ip],pi_hi[ip]}, for successive
    values of {ip}. Also saves the value of {bp[ib].{XX}} (which must
    always be in the specified range) into {pi[ip]}.
    
    Currently considers only the parameters {.id_fin_up} and
    {.id_ini_dn}, in that order. The parameters {.id_ini_sg} and
    {.id_fin_sg} must be identical in all three parameter sets.
    
    The arrays {pi_lo,pi,pi_hi} are allocated by the procedure, and
    their addresses are returned in {*pi_loP,*piP,*pi_hiP}. The number
    {npi} of adjustable parameters found is returned in {*npiP}. */
    

#endif
