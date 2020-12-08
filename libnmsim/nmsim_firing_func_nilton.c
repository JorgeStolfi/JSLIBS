/* See {nmsim_firing_func_nilton.h} */
/* Last edited on 2020-12-04 21:17:17 by jstolfi */

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

#include <nmsim_firing_func_nilton.h>

void nmsim_firing_func_nilton_eval
  ( double V_M,
    double V_D,
    double V,
    double *prP,
    double *dprP
  )
  {
    double V_min = V_M - V_D;
    double V_max = V_M + V_D;
    if ((prP == NULL) && (dprP == NULL)) { return; }
    if (V <= V_min)
      { if (prP != NULL) { (*prP) = 0.0; }
        if (dprP != NULL) { (*dprP) = 0.0; }
      }
    else if (V >= V_max)
      { if (prP != NULL) { (*prP) = 1.0; }
        if (dprP != NULL) { (*dprP) = 0.0; }
      }
    else
      { double z = (V - V_min)/(V_max - V_min);
        if (prP != NULL) { (*prP) = pow(z, nmsim_firing_func_nilton_EXP); }
        if (dprP != NULL) { (*dprP) = 0.4*(*prP)/z/(V_max - V_min); }
      }
  }

double nmsim_firing_func_nilton_eval_inv
  ( double V_M,
    double V_D,
    double pr
  )
  {
    /* Assumes that {pr} is in {[0_1]}. */
    double V_min = V_M - V_D;
    double V_max = V_M + V_D;
    if (pr <= 0.0)
      { return V_min; }
    else if (pr >= 1.0)
      { return V_max; }
    else
      { double z = pow(pr, 1.0/nmsim_firing_func_nilton_EXP);
        double V = V_min + z*(V_max - V_min);
        return V;
      }
  }
