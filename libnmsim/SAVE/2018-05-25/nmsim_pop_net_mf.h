#ifndef nmsim_pop_net_mf_H
#define nmsim_pop_net_mf_H
 
/* Mean-field simulation of multi-population networks of Galves-Löcherbach neurons. */
/* Last edited on 2018-04-11 10:34:13 by jstolfi */

#define _GNU_SOURCE

#include <nmsim_basic.h>
#include <nmsim_neuron_parms.h>
#include <nmsim_bundle_parms.h>
#include <nmsim_pop_net_parms.h>
#include <nmsim_pop_mf_state.h>

void nmsim_pop_net_mf_evolve
  ( nmsim_pop_net_parms_t *parms,  /* Network description. */
    nmsim_pop_mf_state_t state[],  /* Statistical state of each population. */
    int64_t nt                     /* Number of time steps. */
  );
  /* Simulates the evolution of a network with one or more homogeneous
    populations of neurons with stochastic states {net.pso[k]} and
    neuron parameters {parms} from some discrete time {t} to time
    {t+1}.

    In addition to synaptic inputs, each neuron {i} in the 
    population {k} is also assumed to have received, between times
    {t} and {t+1}, an external input signal {I[i,t]},
    which is a normally distributed variable with mean
    and deviation given by {net.exin(k,t,&I_avg,I_dev)}

    ??? Should the input signal be modulated by the input gain
    factor {G[i,t]} too? ???  */

#endif
