/* Last edited on 2013-10-25 22:06:29 by stolfilocal */
/* See {interp_spline.h}. */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <ix.h>
#include <affirm.h>

#include <interp_spline.h>
#include <interp_spline_I.h>
#include <interp_spline_B.h>
#include <interp_spline_O.h>

int interp_spline_compute_num_samples(int ord, interp_spline_kind_t knd)
   { switch(knd)
      { 
        case interp_spline_kind_B:
          return interp_spline_B_compute_num_samples(ord);
        case interp_spline_kind_I:
          return interp_spline_I_compute_num_samples(ord);
        case interp_spline_kind_O:
          return interp_spline_O_compute_num_samples(ord);
        default:
          demand(FALSE, "invalid {knd}");
          return 0;
      }
  }
  
void interp_spline_get_indices(double z, int ns, ix_reduction_t red, int nw, int ix[])
  {
    bool_t debug = FALSE;

    /* Get the raw index of the first source data sample: */
    int k = (int)floor(z - 0.5*(nw - 1));
    
    /* Compute the reduced indices {ix[0..nw-1]}: */
    int i;
    for (i = 0; i < nw; i++) { ix[i] = (int)ix_reduce(k + i, ns, red); }

    if (debug)
      { fprintf(stderr, "z = %7.4f", z);
        for (i = 0; i < nw; i++) { fprintf(stderr, "  ix[%d] = %d\n", i, ix[i]); }
        fprintf(stderr, "\n");
      }
  }
  
void interp_spline_get_weights(double z, int ord, interp_spline_kind_t knd, int nw, double wt[])
  {
    switch(knd)
      { 
        case interp_spline_kind_B:
          interp_spline_B_get_weights(z, ord, nw, wt);
          break;
        case interp_spline_kind_I:
          interp_spline_I_get_weights(z, ord, nw, wt);
          break;
        case interp_spline_kind_O:
          interp_spline_O_get_weights(z, ord, nw, wt);
          break;
        default:
          demand(FALSE, "invalid {knd}");
      }
  }
