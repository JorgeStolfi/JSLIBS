/* See {btc_bubble_nl_opt_gather_integer_variable_parameters.h} */
/* Last edited on 2024-12-05 10:22:55 by stolfi */

#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <btc_bubble_t.h>
#include <btc_is_trivial_int_range.h>

#include <btc_bubble_nl_opt_gather_integer_variable_parameters.h>

void btc_bubble_nl_opt_gather_integer_variable_parameters
  ( int nb, 
    btc_bubble_t bp_lo[], 
    btc_bubble_t bp[], 
    btc_bubble_t bp_hi[], 
    int* npiP, 
    int** pi_loP, 
    int** piP, 
    int** pi_hiP
  )
  {
    demand((bp_lo != NULL) && (bp_hi != NULL), "min and max parameter sets must be given");

    int npi_pbb = 2; /* Max adjustable integer parameters per bubble. */
    
    /* Allocate output vectors: */
    int npi_max = npi_pbb*nb; /* Max adjustable integer parameter in all. */
    int* pi_lo = notnull(malloc(npi_max*sizeof(int)), "no mem");
    int* pi    = notnull(malloc(npi_max*sizeof(int)), "no mem");
    int* pi_hi = notnull(malloc(npi_max*sizeof(int)), "no mem");

    int npi = 0;
    
    auto void collect(int vlo, int v, int vhi, int ib, char* vname);
      /* If the interval {vlo..vhi} is not trivial, appends {vlo,v,vhi}
        to {{pi_lo,pi,pi_hi}[0..npi-1]}, incrementing {npi}.
        In any case, checks that {v} is in {vlo..vhi}.
        The bubble index {ib} and the field name {*vname} 
        are printed if error. */
      
    int ib;
    for (ib = 0; ib < nb; ib++)
      { btc_bubble_t* blo = &(bp_lo[ib]);
        btc_bubble_t* b   = &(bp[ib]);
        btc_bubble_t* bhi = &(bp_hi[ib]);
        
        (void)btc_is_trivial_int_range(blo->id_ini_sg, b->id_ini_sg, bhi->id_ini_sg, ib, ".id_ini_sg", TRUE);
        collect(blo->id_fin_up, b->id_fin_up, bhi->id_fin_up, ib, ".id_fin_up");
        collect(blo->wd_plat, b->wd_plat, bhi->wd_plat, ib, ".wd_plat");
        (void)btc_is_trivial_int_range(blo->id_fin_sg, b->id_fin_sg, bhi->id_fin_sg, ib, ".id_fin_sg", TRUE);
      }
      
    (*npiP) = npi;
    (*pi_loP) = realloc(pi_lo, npi*sizeof(int));
    (*piP) = realloc(pi, npi*sizeof(int));
    (*pi_hiP) = realloc(pi_hi, npi*sizeof(int));
    
    /* INTERNAL IMPLEMENTATIONS */

    void collect(int vlo, int v, int vhi, int ib, char* vname)
      { if (! btc_is_trivial_int_range(vlo, v, vhi, ib, vname, FALSE))
          { assert(npi < npi_max);
            pi_lo[npi] = vlo;
            pi[npi] = v;
            pi_hi[npi] = vhi;
            npi++;
          }
      }
        
  }
        
