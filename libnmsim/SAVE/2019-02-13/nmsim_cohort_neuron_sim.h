#ifndef nmsim_cohort_neuron_sim_H
#define nmsim_cohort_neuron_sim_H

/* Cohort evolution in group-level network simulation. */
/* Last edited on 2019-02-13 17:35:46 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_cohort_neuron_state.h>
#include <nmsim_class_neuron.h>

#include <nmsim_cohort_neuron_sim.h>

void nmsim_cohort_neuron_sim_step
  ( nmsim_cohort_neuron_state_t *cso, 
    double DV_avg, 
    double DV_dev, 
    nmsim_class_neuron_t *nclass, 
    nmsim_cohort_neuron_state_t *csn_fire, 
    nmsim_cohort_neuron_state_t *csn_fail
  ); 
  /* Simulates the evolution of neurons in a specific cohort {S[k,t]} of
    a homogeneous population, from some discrete time {t} to the next
    time {t+1}. The cohort is assumed to have some age {k} and the
    stochastic state {cso} at time {t}. A fraction {cso->rho} of those
    neurons is assumed to fire, and they will join the cohort {S[0,t+1]}
    of age zero at time {t+1}. The new state of those neurons that fire is
    returned in {csn_fire}. The neurons that fail to fire will become
    the cohort {S[k+1,t+1]}, with age {k+1} at time {t+1}, whose new state
    is stored by the procedure into {csn_fail}.

    It is assumed that the total input that each neuron {i} in {S[k,t]}
    received from all its input synapses betwen time {t} and time {t+1}
    is a normal random variable with mean {DV_avg} and deviation
    {DV_dev}. This potential increment includes the pulses {X[j,t]} of
    each input neuron {j} in the network, modulated by its output modulator
    {H[j,t]} and by the rested synaptic gain {w[j-->i]}. It also
    includes the external input signal {I[i,t]}. */
    
#endif
