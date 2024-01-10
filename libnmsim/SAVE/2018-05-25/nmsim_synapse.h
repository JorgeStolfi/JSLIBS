#ifndef nmsim_synapse_H
#define nmsim_synapse_H
 
/* Synapses in a neuron-level network. */
/* Last edited on 2017-08-03 02:22:24 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <vec.h>

#include <nmsim_basic.h>

typedef struct nmsim_synapse_t 
  { nmsim_neuron_ix_t pre; /* Index of pre-synaptic neuron. */
    nmsim_neuron_ix_t pos; /* Index of post-synaptic neuron. */
    double swt;            /* Resting weight of synapse. */
  } nmsim_synapse_t;
  /* Specifies a chemical synapse from neuron {pre} to neuron {pos}. */

vec_typedef(nmsim_synapse_vec_t, nmsim_synapse_vec, nmsim_synapse_t);
  /* Type of an extensible vector of synapses. */

#endif
