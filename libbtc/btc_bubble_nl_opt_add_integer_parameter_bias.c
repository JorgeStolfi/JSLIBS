/* See {btc_bubble_nl_opt_add_integer_parameter_bias.h} */
/* Last edited on 2024-12-05 10:41:53 by stolfi */

#include <btc_bubble_nl_opt_add_integer_parameter_bias.h>

double btc_bubble_nl_opt_add_integer_parameter_bias(double Q, int npi, int pi_a[], int pi_b[], double alpha)
  { int ip;
    double sum_d2 = 0;
    for (ip = 0; ip < npi; ip++)
      { double d = (double)(pi_a[ip] - pi_b[ip]);
        sum_d2 += d*d;
      }
    return Q + alpha*sum_d2;
  }
    
