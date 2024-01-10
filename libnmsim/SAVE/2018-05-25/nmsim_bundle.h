#ifndef nmsim_bundle_H
#define nmsim_bundle_H

/* Parameters for a bundle of connections two populations. */ 
/* Last edited on 2017-07-21 23:35:33 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_neuron_parms.h>
#include <nmsim_cohort.h>
#include <nmsim_pop.h>

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
