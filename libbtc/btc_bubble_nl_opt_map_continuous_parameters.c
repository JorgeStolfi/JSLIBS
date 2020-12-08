/* See {btc_bubble_nl_opt_map_continuous_parameters.h} */
/* Last edited on 2015-04-20 17:21:01 by stolfilocal */

#define _GNU_SOURCE

#include <bool.h>
#include <affirm.h>

#include <btc_bubble_nl_opt_map_continuous_parameters.h>

void btc_bubble_nl_opt_map_continuous_parameters
  ( int npf, 
    double pf_lo[], 
    double pf[], 
    double pf_hi[],
    double x[]
  )
  {
    int ip;
    for (ip = 0; ip < npf; ip++)
      { double vlo = pf_lo[ip];
        double v = pf[ip];
        double vhi = pf_hi[ip];
        demand(vlo < vhi, "invalid range");
        x[ip] = (v - vlo)/(vhi - vlo);
      }
  }

