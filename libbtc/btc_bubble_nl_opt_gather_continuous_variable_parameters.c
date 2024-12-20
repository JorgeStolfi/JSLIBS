/* See {btc_bubble_nl_opt_gather_continuous_variable_parameters.h} */
/* Last edited on 2024-12-05 10:22:53 by stolfi */

#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <btc_bubble_t.h>
#include <btc_is_trivial_double_range.h>

#include <btc_bubble_nl_opt_gather_continuous_variable_parameters.h>

void btc_bubble_nl_opt_gather_continuous_variable_parameters
  ( int nb, 
    btc_bubble_t bp_lo[], 
    btc_bubble_t bp[], 
    btc_bubble_t bp_hi[], 
    int* npfP, 
    double** pf_loP, 
    double** pfP, 
    double** pf_hiP
  )
  {
    demand((bp_lo != NULL) && (bp_hi != NULL), "min and max parameter sets must be given");

    int npf_pbb = 2; /* Max adjustable integer parameters per bubble. */
    
    /* Allocate output vectors: */
    int npf_max = npf_pbb*nb; /* Max adjustable integer parameter in all. */
    double* pf_lo = notnull(malloc(npf_max*sizeof(double)), "no mem");
    double* pf    = notnull(malloc(npf_max*sizeof(double)), "no mem");
    double* pf_hi = notnull(malloc(npf_max*sizeof(double)), "no mem");

    int npf = 0;
    
    auto void collect(double vlo, double v, double vhi, int ib, char* vname);
      /* If the interval {vlo..vhi} is not trivial, appends {vlo,v,vhi}
        to {{pf_lo,pf,pf_hi}[0..npf-1]}, incrementing {npf}.
        In any case, checks that {v} is in {vlo..vhi}.
        The bubble index {ib} and the field name {*vname} 
        are printed if error. */
      
    int ib;
    for (ib = 0; ib < nb; ib++)
      { btc_bubble_t* blo = &(bp_lo[ib]);
        btc_bubble_t* b   = &(bp[ib]);
        btc_bubble_t* bhi = &(bp_hi[ib]);
        
        collect(blo->rt_up, b->rt_up, bhi->rt_up, ib, ".rt_up");
        collect(blo->rt_dn, b->rt_dn, bhi->rt_dn, ib, ".rt_dn");
      }
      
    (*npfP) = npf;
    (*pf_loP) = realloc(pf_lo, npf*sizeof(double));
    (*pfP) = realloc(pf, npf*sizeof(double));
    (*pf_hiP) = realloc(pf_hi, npf*sizeof(double));
    
    /* INTERNAL IMPLEMENTATIONS */

    void collect(double vlo, double v, double vhi, int ib, char* vname)
      { if (! btc_is_trivial_double_range(vlo, v, vhi, ib, vname, FALSE))
          { assert(npf < npf_max);
            pf_lo[npf] = vlo;
            pf[npf] = v;
            pf_hi[npf] = vhi;
            npf++;
          }
      }
        
  }
        
