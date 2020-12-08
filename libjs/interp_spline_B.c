/* Last edited on 2013-10-28 20:21:51 by stolfilocal */
/* See {interp_spline_B.h}. */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>

#include <interp_spline_B.h>
 
int interp_spline_B_compute_num_samples(int ord)
  {
    return ord+2;
  }
 
void interp_spline_B_get_weights(double z, int ord, int nw, double wt[])
  {
    assert(ord >= -1);
    assert(nw == ord+2);
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
        int deg = ord+1;
        int r, k;
        for (k = 1; k <= deg; k++)
          { wt[k] = (u/k)*wt[k-1];
            for (r = k-1; r > 0; r--)
              { wt[r] = ((k-r+u)/k)*wt[r-1] + ((r+v)/k)*wt[r]; }
            wt[0] = (v/k)*wt[0];
          }
      }
  }
