#ifndef nmsim_bundle_parms_H
#define nmsim_bundle_parms_H

/* Parameters for a bundle of connections two populations. */ 
/* Last edited on 2018-03-24 05:24:19 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_basic.h>

typedef struct nmsim_bundle_parms_t
  { double K;     /* Average number of synapses per source neuron. */
    double W_dev; /* Deviation of total rested synapse gain. */
    double W_avg; /* Average total rested synapse gain. */
  } nmsim_bundle_parms_t;
  /* Parameters of a set of synapses from a homogeneous population {A}
    of neurons to another homogeneous population {B} (possibly the
    same as {A}). If the population {B} has {N} neurons, each neuron
    of {A} has a synapse to each neuron of {B} with independent
    probability {K/N}. The strength of that synapse, in its fully
    rested state, is a random variable with log-normal distribution,
    mean {W_avg/K}, and deviation {W_dev/K}. */

#endif
