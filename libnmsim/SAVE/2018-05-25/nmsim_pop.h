#ifndef nmsim_pop_H
#define nmsim_pop_H

/* Mean-field state and tools for a population of homogeneous neurons. */
/* Last edited on 2017-07-24 21:36:02 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_neuron_parms.h>
#include <nmsim_cohort.h>
#include <nmsim_bundle.h>

typedef struct nmsim_pop_mf_state_t
  { int32_t nc;  /* Number of cohorts individually represented. */
    nmsim_cohort_mf_state_t *cs; /* State of neurons in each cohort ({nc+1} elements). */
  } nmsim_pop_mf_state_t;
  /* Stochastic state of a homogeneous population of neurons in
    mean-field simulation. The neurons with firing age {k} have state
    {cs[k]}, for {k} in {0..nc-1}. All neurons with age {nc} or greater
    are lumped together and assumed to have state {cs[nc].} */

void nmsim_pop_mf_evolve
  ( nmsim_pop_mf_state_t *pso,
    double DV_avg,
    double DV_dev,
    nmsim_neuron_parms_t *parms,
    nmsim_pop_mf_state_t *psn
  );
  /* Simulates the evolution of a homogeneous population of neurons
    with stochastic state {pso} and neuron parameters {parms} from some
    discrete time {t} to time {t+1}.

    It is assumed that the total input that each neuron {i} in the
    population received from all its input synapses betwen time {t}
    and time {t+1} is a normal random variable with mean {DV_avg} and
    deviation {DV_dev}. This potential increment includes the pulses
    {X[j,t]} of each input neuron {j} in the network, modulated by its
    output modulator {H[j,t]} and by the rested synaptic gain
    {wfix[j,i]}; but not yet modulated by the input gain {G[i,t]} of
    neuron {i}.

    In addition to synaptic inputs, the potential increment {DV_avg,DV-dev} 
    is also assumed include any external input signal {I[i,t]}
    received by each neuron in that population between times {t} and {t+1}.

    ??? Should the input signal be modulated by the input modulator
    {G[i,t]} too? ??? */

#endif
