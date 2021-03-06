#ifndef nmsim_elem_net_sim_H
#define nmsim_elem_net_sim_H
 
/* Simulation of neuron-level networks of Galves-Löcherbach neurons. */
/* Last edited on 2020-12-17 01:16:35 by jstolfi */

#define _GNU_SOURCE

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_elem_net.h>
#include <nmsim_elem_net_trace.h>
#include <nmsim_elem_net_sim_group_stats.h>

void nmsim_elem_net_sim_step
  ( nmsim_elem_net_t *enet,    /* Network description. */
    nmsim_time_t t,            /* Time at start of step. */
    double V[],                /* Neuron potentials (IN/OUT,mV). */
    nmsim_step_count_t age[],  /* Firing ages of neurons (IN/OUT). */
    double M[],                /* Recharge modulator of each neuron (IN/OUT). */
    double H[],                /* Output modulator of each neuron (IN/OUT). */
    bool_t X[],                /* Firing indicator of each neuron (OUT). */
    double I[],                /* External neuron inputs (IN,mV). */
    double J[],                /* Total input of each neuron (OUT,mV). */
    nmsim_elem_net_trace_t *etrace,          /* Traces of monitored neurons. */
    nmsim_elem_net_sim_group_stats_t *gstats /* Statistics of neuron state and activity per group. */
  );
  /* Simulates the evolution of a network {enet} of {nne=enet.nne} GL
    neurons during the time step from discrete time {t} to time
    {t+1}.
    
    On input, {V[i]} and {age[i]} must be the membrane potential and
    firing age of each neuron {i} at time {t}, for {i} in {0..nne-1};
    which defines its state at that time.
    
    On input, {M[i]} and {H[i]} must be the modulators of the recharge
    factor and of the output synaptic strength of neuron {i} at time
    {t}, respectively, as determined from its firing age {age[i]}.
    On output, they will be updated to time {t+1}, taking into 
    account the possible firing of the neuron in that interval.

    On input, {I[i]} must be the extra voltage increment that will be
    imposed on the potential of neuron {i} by external sources, in
    addition to the synaptic inputs from other neurons in the net,
    between times {t} and {t+1}.
    
    The {double} input variables cannot be {NAN}, and the ages must be 
    non-negative.
    
    On output, {X[i]} will be set to {TRUE} if the procedure decided that neuron {i} fired during
    that time step.
    
    On output, {J[i]} will have been set to the total
    input received by of neuron {i} between
    {t} and {t+1}; which is the external input {I[i]} plus the synaptic inputs
    from other neurons in the net. 
    
    On output, {V[i]} and {age[i]} will have been updated by the
    procedure to the simulated state of the neuron at time {t+1},
    according to the GL evolution recurrence. The modulators {M[i]}
    and {H[i]} will be updated too, to match the new {age[i]}.
    
    The input contents of {X} and {J} are ignored.
    
    If {etrace} is not {NULL}, the procedure also stores into it the
    states and evolution data of the selected neurons relevant to that
    time step. Specifically, let {trk = etrace.trne[k]}, for each {k} in
    {0..etrace.nne-1}, be the trace of a monitored neuron {i}.
    
    If {t} is in {trk.tLo .. trk.tHi}, saves in the trace entry 
    {tst=trk.ts[t - trk.tLo]} corresponding to time {t} the parameters
    {V,age,M,H} that describe the state at time {t}, and the parametes
    {X,I,J} that describe what happened between times {t} and {t+1}.
    Note that the updated parameter {V,age,M,H} are NOT stored in the
    entry {t+1} of {trk}.
    
    If {gstats} is not {NULL}, the procedure accumulates into it per-group statistics
    about the state and activity variables of the neurons, between the times
    {gstats.tLo} and {gstats.tHi}. Namely, it accumulates the parameters
    {V,age,M,H} that describe the state at time {t}, the parametes
    {X,I,J} that describe what happened between times {t} and {t+1},
    and the parameters {V,age} at time {t} for neurons that fired 
    between {t} and {t+1}.  Note that the updated parameters {V,age,M,H} 
    at {t+1} are NOT accumulated.  Note also that the {gstats} structure
    must be initialized with {nmsim_elem_net_sim_group_stats_initialize} 
    and finalized with {nmsim_elem_net_sim_group_stats_finalize}. */

void nmsim_elem_net_sim_compute_modulators
  ( nmsim_elem_net_t *enet, 
    nmsim_time_t t,           /* Time at END of step. */
    nmsim_step_count_t age[], /* Firing age of each neuron (IN). */
    double M[],               /* Recharge modulator of each neuron (OUT). */
    double H[]                /* Output modulator of each neuron at time {t} (OUT). */
  );
  /* On input, {age[i]} must be the non-negative firing age of each neuron {i} at 
    some time {t}.  On output, {M[i]} and {H[i]} are the output and recharge modulators of 
    each neuron {i} at that time, as determined from its age.
    This procedure may be useful at the start of the simulation, when
    the firing ages are generated randomly. */

#endif
