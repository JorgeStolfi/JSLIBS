/* See conv_filter.h */
/* Last edited on 2024-11-23 04:42:18 by stolfi */

#define conv_filter_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#include <stdint.h>

#include <affirm.h>
#include <ix_reduce.h>

#include <conv_filter.h>

void conv_filter
  ( uint64_t nx, double x[],
    ix_reduce_mode_t ixred,
    int64_t skip,
    int64_t step,
    uint64_t nw, double w[],
    uint64_t ny, double y[]
  )
  { demand((nw % 2) == 1, "weight table size must be odd");
    uint64_t hw = (nw - 1)/2;
    for (uint64_t i = 0;  i < ny; i++)
      { double sum_w = 0;
        double sum_wx = 0;
        int64_t kx = (int64_t)skip + step*(int64_t)i; 
        /* {x[kx]} is the sample aligned with window center to compute {y[i]}. */
        for (uint64_t j = 0;  j < nw; j++)
          { /* Identfy the matching sample {x[r]} for weight {w[j]}: */
            int64_t r = kx + (int64_t)j - (int64_t)hw;
            /* If {r} is outside the range {0..nx}, try fixing it as per {ixred}: */
            if ((r < 0) || (r >= nx)) { r = ix_reduce(r, nx, ixred); }
            /* Depending on {ixred}, {r} may still be out of range: */
            if ((r >= 0) && (r < nx))
              { double wj = w[j]; 
                double xr = x[r];
                sum_w += wj; 
                sum_wx += wj*xr;
              }
          }
        y[i] = sum_wx/sum_w;
      }
  }
