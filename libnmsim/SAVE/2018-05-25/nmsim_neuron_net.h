#ifndef nmsim_neuron_net_H
#define nmsim_neuron_net_H
 
/* Neuron-level modeling of networks of Galves-Löcherbach neurons. */
/* Last edited on 2018-04-11 12:47:15 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <nmsim_basic.h>
#include <nmsim_synapse.h>
#include <nmsim_neuron_parms.h>
#include <nmsim_neuron_state.h>

/* !!! remove state from network description !!! */
  
#define nmsim_neuron_net_MAX_NEURONS (((uint64_t)1) << 36)
  /* Max number of neurons in a network. */

#define nmsim_neuron_net_MAX_SYNAPSES (((uint64_t)1) << 36)
  /* Max number of synapses in a network. */

typedef struct nmsim_neuron_net_t
  { nmsim_neuron_count_t nn;       /* Number of neurons. */
    nmsim_synapse_count_t ns;      /* Number of synapses. */
    /* These arrays are indexed with neuron index {0..nn-1}: */
    nmsim_neuron_parms_t **parms;  /* Parameters of each neuron. */
    nmsim_neuron_count_t *deg_in;  /* Number of input synapses of each neuron. */
    nmsim_neuron_count_t *deg_ot;  /* Number of output synapses of each neuron. */
    nmsim_neuron_state_t *state;   /* Current state of each neuron. */
    /* This array is indexed with synapse index {0..ns-1}: */
    nmsim_synapse_vec_t syn;       /* List of synapses. */
  } nmsim_neuron_net_t;
  /* Description of a network with {nn} Galves-Loecherbach neurons,
    with {ns} specific synapses and neuron and synapse parameters.
    
    The pointers {.parms[i]} and {.parms[j]} for distinct neurons {i,j}
    may point to the same parameter record.
    
    The parameters of synapse {s}, including its origin and destination,
    are in {.syn[s]}, for {s} in {0..ns-1}. */
    
/* NETWORK CREATION */

nmsim_neuron_net_t *nmsim_neuron_net_new(nmsim_neuron_count_t nn);
  /* Allocates an {nmsim_neuron_net_t} structure with {nn}
    neurons and initially with zero synapses.
    All parameter pointers are {NULL},
    and all neuron states are undefined. */

void nmsim_neuron_net_add_random_synapses
  ( nmsim_neuron_net_t *net,
    nmsim_neuron_ix_t ia0,
    nmsim_neuron_ix_t ia1,
    nmsim_neuron_ix_t ib0,
    nmsim_neuron_ix_t ib1,
    double pr,
    double swt_tot_avg,
    double swt_tot_dev
  );
  /* Appends to the vector {net->syn[0..ns-1]} a number of random synapses
    from neurons {ia_0..ia_1} to neurons {ib_0..ib_1}.
    
    On input, the variable {net->ns} must contain the number {ns}
    of previously defined synapses. On output, it will be incremented
    with the number of synapses created by the function call. The
    vector {net->syn} must be allocated before the call (even if
    with zero size) and will be extended by the function if and when needed.
    
    The parameter {pr} is basically the probability of adding a synapse
    between any neuron {ia} in {ia0..ia1} and any neuron {ib} in {ib0..ib1}.
    In particular, if {pr} is zero, no synapses are added and the call is
    a no-op.  
    
    If {pr} is positive, usually at most one synapse is added, with
    probability {pr}, for each of those pairs {ia,ib}. Then the number
    of synapses added by the call is a random variable with binomial
    distribution, max value {na*nb} and mean value {pr*na*nb} where
    {na = ia1-ia0+1} and {nb = ib1-ib0+1}. On average, the call will
    increment the output degree of each neuron in {ia0..ia1} by
    {pr*nb} synapses, and the input degree of each neuron in
    {ib0..ib1} by {pr*na} synapses.
    
    Exceptionally, if this process gives zero new synapses to a neuron
    {ib}, but {pr} is not zero, one synapse is added anyway, from a
    neuron {ia} randomly chosen among {ia0..ia1}. Thus, if {pr} is not
    zero, the procedure always adds at least {nb} synapses. This
    exceptional case is unlikely to occur if {pr*na} is substantially
    greater than 1.
    
    The weights of the added synapses will have a log-normal
    distribution. Specifically, the log {zwt} of the absolute weight
    will be a random variable with normal distribution.  The
    sign of the weight will be the same as that of {swt_tot_avg}.
    
    The mean {zavg} and deviation {zdev} of the log-weight {zwt} will
    be chosen so that the SUM of the weights of the synapses entering
    each neuron {ib} in {ib0..ib1} is a random variable with mean
    {swt_tot_avg} and deviation {swt_tot_dev}. In particular, if
    {swt_tot_dev} is zero, every new synapse entering {ib} will have
    the same weight {swt = swt_tot_avg/din} where {din} is the number
    of new input synapses added to {ib}. */

/* OUTPUT */

void nmsim_neuron_net_write(FILE *wr, nmsim_neuron_net_t *net);
  /* Writes a description of the neuron net {net} to file {wr}. */

/* NEURON-LEVEL SIMULATION */  

void nmsim_neuron_net_tot_inputs(nmsim_neuron_net_t *net, bool_t X[], double J[], double dV[]);
  /* Computes the total voltage increments {dV[0..nn-1]} for each
    neuron {i} of network {net} accumulated during a simulation step,
    given the firing indicators {X[0..nn-1]} and the external inputs
    {J[0..nn-1]} of all neurons during that step. */

void nmsim_neuron_net_evolve
  ( void
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

    ??? Should the input signal be modulated by the input modulator
    {G[i,t]} too? ???  */

void nmsim_neuron_net_free(nmsim_neuron_net_t *net);
  /* Releases the storage taken by {net}, including synapse 
    lists, neuron property and state tables, and the top-level
    descriptor itself. */

#endif

