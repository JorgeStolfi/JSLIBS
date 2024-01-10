#ifndef nmsim_test_H
#define nmsim_test_H
 
/* Test tools for the {nmsim} library. */
/* Last edited on 2017-08-04 11:47:34 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_neuron_parms.h>

void nmsim_test_plot_firing_func(char *prefix, nmsim_neuron_parms_t *parms);
  /* Writes values of the firing function {parms->Phi}
    and its derivative, for various potentials, to file 
    "{prefix}-Phi.txt". */

#endif


