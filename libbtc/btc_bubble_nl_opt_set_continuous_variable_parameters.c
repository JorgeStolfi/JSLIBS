/* See {btc_bubble_nl_opt_set_continuous_variable_parameters.h} */
/* Last edited on 2015-04-21 01:52:57 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <btc_bubble_t.h>

#include <btc_bubble_nl_opt_set_continuous_variable_parameters.h>


void btc_bubble_nl_opt_set_continuous_variable_parameters
  ( int npf, 
    double pf[], 
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
        if (blo->rt_up < bhi->rt_up) { assert(ip < npf); b->rt_up = pf[ip]; ip++; }
        if (blo->rt_dn < bhi->rt_dn) { assert(ip < npf); b->rt_dn = pf[ip]; ip++; }
      }
      
    assert(ip == npf);
  }

