#ifndef nmsim_test_H
#define nmsim_test_H
 
/* Test tools for the {nmsim} library. */
/* Last edited on 2019-03-28 15:57:48 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include <nmsim_firing_func.h>

void nmsim_test_firing_func
  ( nmsim_firing_func_t *Phi,
    double r,
    int32_t ns
  );
  /* Generates a plot data file of the firing function {Phi} and its derivative.
    Also checks the derivative numerically, and the inverse.
    
    The plotted range of the potential is from {Phi.V_M - r*Phi.V_D} to {Phi.V_M + r*Phi.V_D}.
    
    The file name is "out/Phi_{c}_M{MMM.MM}_D{DDD.DD}.txt" where {c} is the
    function's class {Phi->class} ('G', 'L', etc.), {MMM.MM} is the midpoint potential
    {Phi->V_M} (with two decimals and sign, zero-padded), and {DDD.DD} is 
    the deviation {Phi->V_D} (with two decimals, unsigned, zero-padded). */

double *nmsim_test_NAN_vector(size_t n);
   /* Returns a newly allocated vector of {double}s, all set to {NAN}. */

#endif


