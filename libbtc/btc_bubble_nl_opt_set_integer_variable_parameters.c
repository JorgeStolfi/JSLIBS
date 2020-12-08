/* See {btc_bubble_nl_opt_set_integer_variable_parameters.h} */
/* Last edited on 2015-04-29 23:59:08 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <btc_bubble_t.h>

#include <btc_bubble_nl_opt_set_integer_variable_parameters.h>


void btc_bubble_nl_opt_set_integer_variable_parameters
  ( int npi, 
    int pi[], 
    int nb, 
    btc_bubble_t bp_lo[], 
    btc_bubble_t bp[], 
    btc_bubble_t bp_hi[]
  )
  {
    demand((bp_lo != NULL) && (bp_hi != NULL), "min and max parameter sets must be given");

    int ip = 0; /* Number of parameters already set. */
    
    int ib;
    for (ib = 0; ib < nb; ib++)
      { btc_bubble_t* blo = &(bp_lo[ib]);
        btc_bubble_t* b   = &(bp[ib]);
        btc_bubble_t* bhi = &(bp_hi[ib]);
        if (blo->id_fin_up < bhi->id_fin_up) { assert(ip < npi); b->id_fin_up = pi[ip]; ip++; }
        if (blo->wd_plat < bhi->wd_plat) { assert(ip < npi); b->wd_plat = pi[ip]; ip++; }
      }
      
    assert(ip == npi);
  }

