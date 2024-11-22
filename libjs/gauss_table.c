/* See gauss_table.h */
/* Last edited on 2024-11-19 06:08:08 by stolfi */

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <values.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <gauss_table.h>
#include <gauss_bell.h>

double *gauss_table_make(uint32_t n, double avg, double dev, bool_t normSum, bool_t folded)
  {
    demand(dev >= 0, "invalid standard deviation");

    /* Allocate table: */
    double *w = talloc(n, double);

    /* Fill the table: */
    if (folded)
      { for (int32_t i = 0; i < n; i++) 
          { w[i] = gauss_table_folded_bell((double)i - avg, dev, n); }
      }
    else
      { for (int32_t i = 0; i < n; i++) 
          { w[i] = gauss_bell_eval((double)i, avg, dev); }
      }
      
    if (normSum)
      { /* Normalize to unit sum: */
        double sum = 0.0; 
        for (int32_t i = 0; i < n; i++) { sum += w[i]; }
        demand(sum > 0, "cannot normalize a zero-sum table");
        for (int32_t i = 0; i < n; i++) { w[i] /= sum; }
      }
    else if (folded)
      { /* Normalize to unity at {avg}: */
        double wmax = gauss_table_folded_bell(0.0, dev, n);
        assert(wmax > 0.0);
        for (int32_t i = 0; i < n; i++) { w[i] /= wmax; }
      }
    return w;
  }
  
double gauss_table_folded_bell(double z, double dev, uint32_t n)
  {
    assert(dev >= 0);
    if ((n > 0) && (dev > gauss_table_BIG_DEV*n))
      { /* The folded bell is flat to 10^-12 or more: */
        return 1.0;
      }
    else
      { if ((n > 0) && ((z < 0) || (z >= n)))
          { /* Reduce {z} to the range {[0_n)}: */
            z = z - n*floor(z/n);
            while (z >= n) { z -= n; }
            while (z < 0) { z += n; }
            assert(z >= 0);
            assert(z < n);
          }
          
        /* Check for degenerate bell: */
        if (dev == 0.0) { return (z == 0.0 ? 1.0 : 0.0); }

        /* !!! Find a faster formula !!! */

        /* Add all fold-over terms that are {10^{-16}} or more: */
        int32_t kmax = (n == 0 ? 0 : (int32_t)ceil((gauss_bell_BIG_ARG*dev + z)/n));
        assert(kmax >= 0);
        double sum = 0;
        
        /* Add all terms {z ± kmax*n} to {z ± n}, from small to large: */
        int32_t k = kmax;
        while (k > 0)
          { double um = (z - k*(int32_t)n)/dev; 
            double up = (z + k*(int32_t)n)/dev; 
            sum += exp(-um*um/2) + exp(-up*up/2);
            k--;
          }
        /* Add term for {z}: */
        double u = z/dev;
        sum += exp(-u*u/2);
        return sum;
      }
  }
