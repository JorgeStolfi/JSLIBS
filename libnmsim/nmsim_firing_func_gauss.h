#ifndef nmsim_firing_func_gauss_H
#define nmsim_firing_func_gauss_H

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_firing_func.h>

/* The Gaussian firing function {\Phi^G}. */
/* Last edited on 2019-01-02 23:04:02 by jstolfi */
 
void nmsim_firing_func_gauss_eval
  ( double V_M,
    double V_D,
    double V, 
    double *prP,
    double *dprP
  );
  /* Evaluates the Gaussian firing function {\Phi^G} at {V};
    with is the integral of the Gaussian distribution with 
    mean {V_M} and deviation {V_D}. */

double nmsim_firing_func_gauss_eval_inv
  ( double V_M,
    double V_D,
    double pr
  );
  /* Evaluates the inverse of {\Phi^G} at {pr}. */

#endif
