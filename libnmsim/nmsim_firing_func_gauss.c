/* See {nmsim_firing_func_gauss.h} */
/* Last edited on 2019-01-02 23:54:15 by jstolfi */

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

#include <nmsim_firing_func_gauss.h>

#define JSM_SQRT2PI  (2.50662827463100050240) 
  /* Constant sqrt(2*pi) */

#define nmsim_firing_func_gauss_Z_SAT (6.0)
 /* Value of {z} parameter for which the Gaussian Phi saturates. */

void nmsim_firing_func_gauss_eval
  ( double V_M,
    double V_D,
    double V,
    double *prP,
    double *dprP
  )
  {
    double z = M_SQRT1_2*(V - V_M)/V_D;
    if (fabs(z) >= nmsim_firing_func_gauss_Z_SAT)
      { if (prP != NULL) { *(prP) = (z < 0 ? 0.0 : 1.0); }
        if (dprP != NULL) { (*dprP) = 0.0; }
      }
    else
      { if (prP != NULL) { (*prP) = (erf(z) + 1)/2; }
        if (dprP != NULL) { (*dprP) = exp(-z*z)/V_D/JSM_SQRT2PI; }
      }
  }

double nmsim_firing_func_gauss_eval_inv
  ( double V_M,
    double V_D,
    double pr
  )
  {
    /* Assumes that {pr} is in {[0_1]}. */
    double zsat = nmsim_firing_func_gauss_Z_SAT;
    double z;
    if (pr <= 0)
      { z = -INF; }
    else if (pr >= 1)
      { z = +INF; }
    else
      { z = erf_inv(2*pr - 1); }
      
    /* Clip {z} to the significant range: */
    if (z < -zsat) { z = -zsat; }
    if (z > +zsat) { z = +zsat; }

    /* Convert to potential: */
    double V = V_M + V_D * z*M_SQRT2;
    return V;
  }
