#ifndef nmsim_firing_func_H
#define nmsim_firing_func_H

#define _GNU_SOURCE
#include <stdint.h>
 
/* Firing functions for the Galves-Löcherbach neuron model. */
/* Last edited on 2017-08-02 20:30:13 by jstolfi */

typedef struct nmsim_firing_func_t nmsim_firing_func_t;
  /* A descriptor of of a firing function {Phi}. */
  
struct nmsim_firing_func_t *nmsim_firing_func_new
  ( char *class, 
    double deg, 
    double V_M, 
    double D_M
  );
  /* Returns a descriptor for a firing function of the given {class},
    with midpoint potential {V_M} (in millivoltds) and slope {1/D_M}
    at that potential, and degree {deg}. The current classes are
    described in {nmsim_firing_func_INFO}. */

void nmsim_firing_func_eval
  ( struct nmsim_firing_func_t *Phi,
    double V, 
    double *prP,
    double *dprP
  );
  /* Evaluates the firing function {Phi} at the argument {V} (which
    should have been multiplied by the firing modulator).
    
    Stores the firing probability {Phi(V)} into {*prP}, and the derivative 
    of {Phi} at that point into {*dprP}.
    
    If {dprP} is {NULL}, the derivative will not be computed. */

double nmsim_firing_func_eval_inv(struct nmsim_firing_func_t *Phi, double pr);
  /* The potential {V} such that {Phi)V) = pr}.  
    
    Fails if {pr < 0} or {pr > 1}.  Otherwise returns 
    {-INF} or {+INF} if {Phi(V)} is greater than or less than {pr},
    respectively, for all finite {V}. */

#define nmsim_firing_func_INFO \
  "??? The available firing functions are \"gauss\" (integral of Gaussian bell) ... ???"

#endif
