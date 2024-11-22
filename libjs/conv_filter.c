/* See conv_filter.h */
/* Last edited on 2024-11-15 19:12:11 by stolfi */

#define conv_filter_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#include <stdint.h>
#include <affirm.h>

#include <conv_filter.h>

void conv_filter(int32_t nx, double x[], int32_t skip, int32_t step, int32_t nw, double w[], int32_t ny, double y[])
  { demand((nw % 2) == 1, "weight table size must be odd");
    int32_t hw = (nw - 1)/2;
    demand(step > 0, "invalid {step}");
    for (int32_t ky = 0; ky < ny; ky++)
      { double sum_w = 0, sum_wx = 0;
        int32_t kx = skip + step*ky; /* Sample of {x} aligned with {y[ky]}. */
        demand((kx >= 0) && (kx < nx), "invalid index into old seq");
        for(int32_t j = 0; j < nw; j++)
          { int32_t ix = kx + j - hw;
            if ((ix >= 0) && (ix < nx))
              { double wj = w[j]; 
                double xi = x[ix];
                sum_w += wj; 
                sum_wx += wj*xi;
              }
          }
        y[ky] = sum_wx/sum_w;
      }
  }
