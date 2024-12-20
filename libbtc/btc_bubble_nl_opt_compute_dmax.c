/* See {btc_bubble_nl_opt_compute_dmax.h} */
/* Last edited on 2024-12-05 10:22:48 by stolfi */

#include <math.h>

#include <bool.h>
#include <affirm.h>

#include <btc_bubble_nl_opt_compute_dmax.h>

double btc_bubble_nl_opt_compute_dmax(int npf, double x_ini[])
  {
    double sum_d2 = 0;
    int ip;
    for (ip = 0; ip < npf; ip++)
      { double xi = x_ini[ip];
        double d = (xi > 0 ? xi + 1 : 1 - xi);
        sum_d2 += d*d;
      }
    return sqrt(sum_d2);
  }
