/* See {nmsim_firing_func_lramp.h} */
/* Last edited on 2019-01-02 23:07:56 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <nmsim_firing_func.h>

#include <nmsim_firing_func_lramp.h>

void nmsim_firing_func_lramp_eval
  ( double V_M,
    double V_D,
    double V,
    double *prP,
    double *dprP
  )
  {
    double z = (V - V_M)/V_D;
    if (fabs(z) >= 1.0)
      { if (prP != NULL) { *(prP) = (z < 0 ? 0.0 : 1.0); }
        if (dprP != NULL) { (*dprP) = 0.0; }
      }
    else
      { if (prP != NULL) { (*prP) = (1 + z)/2; }
        if (dprP != NULL) { (*dprP) = 0.5/V_D; }
      }
  }

double nmsim_firing_func_lramp_eval_inv
  ( double V_M,
    double V_D,
    double pr
  )
  {
    /* Assumes that {pr} is in {[0_1]}. */
    double z;
    if (pr <= 0)
      { z = -INF; }
    else if (pr >= 1)
      { z = +INF; }
    else
      { z = 2*(pr - 0.5); }
      
    /* Clip {z} to the significant range: */
    if (z < -1.0) { z = -1.0; }
    if (z > +1.0) { z = +1.0; }

    /* Convert to potential: */
    double V = V_M + V_D * z;
    return V;
  }
