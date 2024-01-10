#ifndef nmsim_network_H
#define nmsim_network_H
 
/* Types and functions for the Galves-Löcherbach neuron model. */
/* Last edited on 2019-01-10 04:28:21 by jstolfi */

#define _GNU_SOURCE

#include <nmsim_neuron_elem.h>
#include <nmsim_cohort.h>
#include <nmsim_pop.h>

typedef struct nmsim_network_bundle_parms_t
  { double W_dev; /* Deviation of total rested synapse gain. */
    double W_avg; /* Average total rested synapse gain. */
  } nmsim_network_bundle_parms_t;
  /* Parameters of a set of synapses from a homogeneous population {A}
    of neurons to another homogeneous population {B} (possibly the
    same as {A}). If the population {B} has {N} neurons, each neuron
    of {A} has a synapse to each neuron of {B} with independent
    probability {K/N}. The strength of that synapse, in its fully
    rested state, is a random variable with log-normal distribution,
    mean {W_avg/K}, and deviation {W_dev/K}. */

typedef struct nmsim_network_t
  { int32_t np;  /* Number of populations. */
    nmsim_neuron_class_t *parms; /* Neuron parameters in each population ({np} elements). */
    nmsim_network_bundle_parms_t *conn; /* Connection between populations ({np*np} elements). */
  } nmsim_network_t;
  /* Describes a network with {np} homogeneous populations of neurons,
    with random connections between each pair of populations.
    
    The data about population {r} is {.parms[r]}, for {r} in {0..np-1}.
    The parameters of the distribution of synapses from population {r}
    to population {s} are {.conn[r + s*np]}, for {r,s} in {0.np-1}. Note
    that a population usually has connections to itself too.

/* MEAN FIELD SIMULATION */

void nmsim_network_mf_evolve
  (
  );
  /* Simulates the evolution of a network with one or more homogeneous
    populations of neurons with stochastic states {net.pso[k]} and
    neuron parameters {parm} from some discrete time {t} to time
    {t+1}.

    In addition to synaptic inputs, each neuron {i} in the 
    population {k} is also assumed to have received, between times
    {t} and {t+1}, an external input signal {I[i,t]},
    which is a normally distributed variable with mean
    and deviation given by {net.exin(k,t,&I_avg,I_dev)}

    ??? Should the input signal be modulated by the input gain
    factor {G[i,t]} too? ???  */

#endif
