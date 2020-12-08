#ifndef nmsim_firing_func_lramp_H
#define nmsim_firing_func_lramp_H

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_firing_func.h>

/* The linear ramp firing function {\Phi^L}. */
/* Last edited on 2019-01-02 23:04:17 by jstolfi */
 
void nmsim_firing_func_lramp_eval
  ( double V_M,
    double V_D,
    double V, 
    double *prP,
    double *dprP
  );
  /* Evaluates the linear firing function {\Phi^L} at {V};
    that rises linearly from 0 at {V_M - V_D} to 1 at {V_M + V_D}. */

double nmsim_firing_func_lramp_eval_inv
  ( double V_M,
    double V_D,
    double pr
  );
  /* Evaluates the inverse of {\Phi^L} at {pr}. */

#endif
