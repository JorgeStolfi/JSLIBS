#ifndef nmsim_elem_net_sim_H
#define nmsim_elem_net_sim_H
 
/* Simulation of neuron-level networks of Galves-Löcherbach neurons. */
/* Last edited on 2019-02-13 17:40:48 by jstolfi */

#define _GNU_SOURCE

#include <nmsim_basic.h>
#include <nmsim_elem_net.h>
#include <nmsim_elem_net_trace.h>

void nmsim_elem_net_sim_step
  ( nmsim_elem_net_t *enet,    /* Network description. */
    nmsim_time_t t,            /* Time at start of step. */
    double I[],                /* External neuron inputs (IN). */
    double V[],                /* Neuron potentials (IN/OUT). */
    nmsim_step_count_t age[],  /* Firing ages of neurons (IN/OUT). */
    bool_t X[],                /* Firing indicator of each neuron (OUT). */
    double H[],                /* Output modulator of each neuron (OUT). */
    double M[],                /* Recharge modulator of each neuron (OUT). */
    double J[],                /* Total input of each neuron (OUT). */
    nmsim_elem_net_trace_t *etrace /* Monitored neuron traces. */
  );
  /* Simulates the evolution of a network {enet} of {nne=enet.nne} GL
    neurons during the time step from discrete time {t} to time
    {t+1}.
    
    On input, {V[i]} and {age[i]} must be the membrane potential and firing age
    of each neuron {i} at time {t}, for {i} in {0..nne-1}; which defines its
    state at that time. These
    parameters are updated by the procedure to reflect the state of
    each neuron at time {t+1}.

    In addition to the synaptic inputs, each neuron {i} is also
    assumed to have received, in that time step, an extra input
    current that, by itself, would raise its potential by {I[i]}.
    
    The procedure ets {X[i]} to {TRUE} if it decided that neuron {i} fired during
    that time step.
    
    The procedure also computes the recharge modulator {M[i]} and output
    modulator {H[i]} for each neuron at time {t}. Note that the 
    recharge modulator is irrelevant if the neuron fires in that
    time step. 
    
    The procedure also sets {J[i]} to the total
    input (synaptic inputs plus external input) of neuron {i} between
    {t} and {t+1}. 
    
    The input contents of {X,M,H,J} are ignored.
    
    If {etrace} is not {NULL}, the procedure also stores into it 
    the states of the selected neurons for that time step. Specifically, for each monitored
    neuron trace {trne} listed in {etrace.trne[0..etrace.nne-1]} such 
    that {t} is in {trne.tini .. trne.tfin}, saves in the entry of {trne.st} corresponding to
    time {t} a record with {V,age,M,H} at time {t}, and {X,I,J} between
    times {t} and {t+1}.  Note that the entries of {etrace} corresponding
    to time {t+1} are NOT set. */

#endif
