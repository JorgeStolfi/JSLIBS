/* See conv_filter.h */
/* Last edited on 2014-07-26 23:04:58 by stolfilocal */

#define conv_filter_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE

#include <affirm.h>

#include <conv_filter.h>

void conv_filter(int nx, double x[], int skip, int step, int nw, double w[], int ny, double y[])
  { demand((nw % 2) == 1, "weight table size must be odd");
    int hw = (nw - 1)/2;
    demand(step > 0, "invalid {step}");
    int ky;
    for (ky = 0; ky < ny; ky++)
      { double sum_w = 0, sum_wx = 0;
        int kx = skip + step*ky; /* Sample of {x} aligned with {y[ky]}. */
        demand((kx >= 0) && (kx < nx), "invalid index into old seq");
        int j;
        for(j = 0; j < nw; j++)
          { int ix = kx + j - hw;
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
