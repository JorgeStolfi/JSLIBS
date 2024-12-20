/* See {btc_bubble_nl_opt_unmap_continuous_parameters.h} */
/* Last edited on 2024-12-05 10:23:05 by stolfi */

#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <btc_bubble_nl_opt_unmap_continuous_parameters.h>

double btc_bubble_nl_opt_unmap_continuous_parameters
  ( int npf,
    double x[], 
    double pf_lo[], 
    double pf[], 
    double pf_hi[]
  )
  {
    double dBox2 = 0;
    int ip;
    for (ip = 0; ip < npf; ip++)
      { 
        double xi = x[ip];
        double vlo = pf_lo[ip];
        double v;
        double vhi = pf_hi[ip];
        demand(vlo < vhi, "invalid range");
        if (xi <= -1.0)
          { double dx = xi + 1.0; dBox2 += dx*dx; v = vlo; }
        else if (xi >= +1.0)
          { double dx = xi - 1.0; dBox2 += dx*dx; v = vhi; }
        else
          { v = 0.5*((1.0 - xi)*vlo + (xi + 1.0)*vhi); 
            assert((vlo <= v) && (v <= vhi));
          }
        pf[ip] = v;
      }
    return dBox2;
  }
