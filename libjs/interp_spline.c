/* Last edited on 2024-11-23 06:06:26 by stolfi */
/* See {interp_spline.h}. */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <ix_reduce.h>
#include <affirm.h>

#include <interp_spline.h>
#include <interp_spline_I.h>
#include <interp_spline_B.h>

uint32_t interp_spline_compute_num_samples(int32_t ord, interp_spline_kind_t knd)
   { switch(knd)
      { 
        case interp_spline_kind_B:
          return interp_spline_B_compute_num_samples(ord);
        case interp_spline_kind_I:
          return interp_spline_I_compute_num_samples(ord);
        default:
          demand(FALSE, "invalid {knd}");
          return 0;
      }
  }
  
void interp_spline_get_indices(double z, uint32_t ns, ix_reduce_mode_t red, uint32_t nw, int32_t ix[])
  {
    bool_t debug = FALSE;

    /* Get the raw index of the first source data sample: */
    int32_t k = (int32_t)floor(z - 0.5*(nw - 1));
    
    /* Compute the reduced indices {ix[0..nw-1]}: */
    for (int32_t i = 0;  i < nw; i++) { ix[i] = (int32_t)ix_reduce(k + i, ns, red); }

    if (debug)
      { fprintf(stderr, "z = %7.4f", z);
        for (uint32_t i = 0;  i < nw; i++) { fprintf(stderr, "  ix[%d] = %d\n", i, ix[i]); }
        fprintf(stderr, "\n");
      }
  }
  
void interp_spline_get_weights(double z, int32_t ord, interp_spline_kind_t knd, uint32_t nw, double wt[])
  {
    switch(knd)
      { 
        case interp_spline_kind_B:
          interp_spline_B_get_weights(z, ord, nw, wt);
          break;
        case interp_spline_kind_I:
          interp_spline_I_get_weights(z, ord, nw, wt);
          break;
        default:
          demand(FALSE, "invalid {knd}");
      }
  }
