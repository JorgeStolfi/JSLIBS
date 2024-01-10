#ifndef nmsim_cohort_H
#define nmsim_cohort_H
 
/* Types and functions for the Galves-Löcherbach neuron model. */
/* Last edited on 2019-01-11 18:06:49 by jstolfi */

#define _GNU_SOURCE

#include <nmsim_neuron_elem.h>

typedef struct nmsim_cohort_state_t
  { double V_avg; /* Mean potential. */
    double V_dev; /* Deviation of potential. */
    double G;     /* Firing gain. */
    double H;     /* Output gain. */
    double eta;   /* Fraction of population that is in this cohort. */
    double rho;   /* Fraction of this cohort that will fire next. */
  } nmsim_cohort_state_t;
  /* Stochastic state of a /cohort/, a set of neurons with the 
    same firing age {k} within some homogeneous population,
    at some discrete time {t}.  
    
    The neurons are assumed to have a normal distribution of potentials,
    with mean {V_avg} and deviation {V_dev}. They are assumed to have
    the same firing gain factor {G}, and same output gain factor {H}.
    
    The cohort comprises a fraction {eta} of the population. A fraction
    {rho} of the cohort will fire during the next time step. */

void nmsim_cohort_state_set
  ( nmsim_cohort_state_t *cs, 
    double V_avg,
    double V_dev,
    double G,
    double H,
    double eta,
    nmsim_neuron_class_t *parms
  );
  /* Sets the state parameters of a cohort {cs} to potential {V_avg ±
    V_dev}, input gain factor {G}, output gain factor {H}, and relative
    fraction {eta}.
    
    The firing fraction {rho} is set based on the firing function
    {parms->Phi}, the potential distribution {V_avg,V_dev},
    modulated by the firing gain {G}. 
    
    Currently, the procedure assumes that {G*V_dev} is small enough for
    {parms->Phi} to be approximated by a linear function at {G*V_avg}.
    
    !!! The procedure should try to handle the case {V_dev} large. !!! */
    
/* MEAN FIELD SIMULATION */

void nmsim_cohort_mf_evolve
  ( nmsim_cohort_state_t *cso, 
    double DV_avg,
    double DV_dev,
    nmsim_neuron_class_t *parms, 
    nmsim_cohort_state_t *csn_fire,
    nmsim_cohort_state_t *csn_fail
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
    each input neuron {j} in the network, modulated by its output gain
    factor {H[j,t]} and by the rested synaptic gain {w[j-->i]}. It also
    includes the external input signal {I[i,t]}. This potential
    increment should not be modulated by the firing gain {G[i,t]} of
    neuron {i}. */

#endif
