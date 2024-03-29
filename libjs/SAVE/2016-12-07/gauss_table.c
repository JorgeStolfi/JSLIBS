/* See gauss_table.h */
/* Last edited on 2011-05-29 10:44:42 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <values.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <gauss_table.h>

double *gauss_table_make(int n, double avg, double dev, bool_t normSum)
  {
    demand(n >= 0, "invalid table size");
    demand(dev >= 0, "invalid standard deviation");

    /* Allocate table: */
    double *w = notnull(malloc(n*sizeof(double)), "no mem");

    /* Fill the table: */
    int i;
    for (i = 0; i < n; i++) 
      { w[i] = gauss_table_folded_bell((double)i - avg, dev, n); }
      
    if (normSum)
      { /* Normalize to unit sum: */
        double sum = 0.0; 
        for (i = 0; i < n; i++) { sum += w[i]; }
        demand(sum > 0, "cannot normalize a zero-sum table");
        for (i = 0; i < n; i++) { w[i] /= sum; }
      }
    else
      { /* Normalize to unit at zero: */
        double wmax = gauss_table_folded_bell(0.0, dev, n);
        assert(wmax > 0.0);
        for (i = 0; i < n; i++) { w[i] /= wmax; }
      }
    return w;
  }
  
double gauss_table_bell(double z, double dev)
  {
    assert(dev >= 0);
    if (dev == 0.0)
      { return (z == 0 ? 1.0 : 0.0); }
    else if (z >= gauss_table_BIG_ARG*dev)
      { return 0.0; }
    else if (z < gauss_table_TINY_ARG*dev)
      { return 1.0; }
    else
      { double u = z/dev;
        return exp(-u*u/2); 
      }
  }
  
double gauss_table_folded_bell(double z, double dev, int n)
  {
    assert(dev >= 0);
    if ((n > 0) && (dev > 1.5*n))
      { /* The folded bell is flat to 10^-12 or more: */
        return 1.0;
      }
    else
      { if ((n > 0) && ((z < 0) || (z >= n)))
          { /* Reduce {z} to the range {[0_n)}: */
            z = z - n*(int)floor(z/n);
            while (z >= n) { z -= n; }
            assert(z >= 0);
            assert(z < n);
          }
          
        /* Check for degenerate bell: */
        if (dev == 0.0) { return (z == 0.0 ? 1.0 : 0.0); }

        /* !!! Find a faster formula !!! */

        /* Add all fold-over terms that are {10^{-16}} or more: */
        int kmax = (n == 0 ? 0 : (int)ceil ((gauss_table_BIG_ARG*dev + z)/n));
        assert(kmax >= 0);
        double sum = 0;
        
        /* Add all terms for {z - kmax*n} to {z + kmax*n}: */
        int k = kmax;
        while (k > 0)
          { double um = (z - k*n)/dev; 
            double up = (z + k*n)/dev; 
            sum += exp(-um*um/2) + exp(-up*up/2);
            k--;
          }
        double uz = z/dev;
        sum += exp(-uz*uz/2);
        return sum;
      }
  }
