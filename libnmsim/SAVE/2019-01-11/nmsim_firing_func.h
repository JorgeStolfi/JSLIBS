#ifndef nmsim_firing_func_H
#define nmsim_firing_func_H
#define _GNU_SOURCE
 
/* Firing functions for the Galves-Löcherbach neuron model. */
/* Last edited on 2016-08-11 11:05:08 by stolfilocal */

typedef struct nmsim_firing_func_t nmsim_firing_func_t;
  /* A descriptor of of a firing function {Phi}. */

void nmsim_firing_func_eval
  ( struct nmsim_firing_func_t *Phi,
    double V, 
    double *pr,
    double *dpr
  );
  /* Evaluates the firing function {Phi} at the argument {V} (which
    should have been multiplied by the firing gain factor).
    
    If {dpr} is {NULL}, the derivative will not be computed. */

nmsim_firing_func_t *nmsim_firing_func_gauss(double V_M, double D_M);
  /* Returns a descriptor for a Gaussian-integral firing
    function {Phi} with midpoint potential {V_M} and slope
    {D_M} at that potential. */

#define nmsim_firing_func_INFO \
  "???"

#endif
