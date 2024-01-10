#ifndef nmsim_firing_func_gauss_H
#define nmsim_firing_func_gauss_H

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_firing_func.h>
 
nmsim_firing_func_t *nmsim_firing_func_gauss_new(double V_M, double D_M);
  /* Returns a descriptor for a Gaussian-integral firing
    function {Phi} with midpoint potential {V_M} and slope
    {1/D_M} at that potential. */

#endif
