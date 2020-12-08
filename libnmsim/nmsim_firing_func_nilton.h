#ifndef nmsim_firing_func_nilton_H
#define nmsim_firing_func_nilton_H

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_firing_func.h>

/* Nilton's firing function {\Phi^N}. */
/* Last edited on 2020-12-03 11:24:44 by jstolfi */
 
#define nmsim_firing_func_nilton_EXP (0.4)
  /* Exponent of Nilton's firing function. */

void nmsim_firing_func_nilton_eval
  ( double V_M,
    double V_D,
    double V, 
    double *prP,
    double *dprP
  );
  /* Evaluates Nilton's firing function {\Phi^N} at {V};
    namely {((V - V_min)/(V_max - V_min))^E} for {V_min < V < V_max},
    where {V_min = V_M-V_D}, {V_max = V_M+V_D}, and 
    {E = nmsim_firing_func_nilton_EXP}. */

double nmsim_firing_func_nilton_eval_inv
  ( double V_M,
    double V_D,
    double pr
  );
  /* Evaluates the inverse of {\Phi^N} at {pr}. */

#endif
