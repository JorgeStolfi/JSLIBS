/* Last edited on 2024-11-23 06:07:05 by stolfi */
/* See {interp_spline_B.h}. */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>

#include <interp_spline_B.h>
 
uint32_t interp_spline_B_compute_num_samples(int32_t ord)
  {
    demand(ord >= -1, "invalid {ord}");
    return (uint32_t)(ord+2);
  }
 
void interp_spline_B_get_weights(double z, int32_t ord, uint32_t nw, double wt[])
  {
    demand(ord >= -1, "invalid {ord}");
    demand(nw == ord+2, "bad window width {nw}");
    wt[0] = 1.0;
    if (ord >= 0)
      { /* Shift {z} by {1/2} depending on parity of {ord}: */
        z = z - 0.5*(ord + 1);

        /* Get the fractional part {u} of {z} and its complement {v}: */
        double u = z - floor(z);
        assert(u >= 0.0);
        assert(u <= 1.0);
        double v = 1.0 - u;
        
        /* Compute the interpolation weights {wt[0..nw-1]} recursively: */
        int32_t deg = ord+1;
        for (int32_t k = 1;  k <= deg; k++)
          { wt[k] = (u/k)*wt[k-1];
            for (int32_t r = k-1; r > 0; r--)
              { wt[r] = ((k-r+u)/k)*wt[r-1] + ((r+v)/k)*wt[r]; }
            wt[0] = (v/k)*wt[0];
          }
      }
  }
