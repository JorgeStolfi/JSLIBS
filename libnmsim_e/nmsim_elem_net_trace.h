#ifndef nmsim_elem_net_trace_H
#define nmsim_elem_net_trace_H
 
/* Trace of an elem-level simulation of a Galves-Löcherbach neuron net. */
/* Last edited on 2020-12-17 13:43:02 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_neuron_trace.h>

typedef struct nmsim_elem_net_trace_t
  { nmsim_time_t tLo;                   /* Discrete time of first neuron trace entries. */
    nmsim_time_t tHi;                   /* Discrete time of last neuron trace entries. */
    nmsim_elem_neuron_count_t ntr;      /* Number of monitored neuron sets. */
    nmsim_elem_neuron_trace_t **trne;   /* Traces of monitored neuron sets. */
  } nmsim_elem_net_trace_t;
  /* A data structure that records the trace of selected neurons or sets of neurons for a certain 
    time interval, in an elem-level simulation of a GL network. 
    
    The traces monitored neurons are stored in {trne[0..ntr-1]}. For
    each {k} in {0..ntr-1}, the trace {trnek=trne[k]} records, for each
    time {t} in {trnek.tLo .. trnek.tHi}, the average state and
    evolution variables of the neurons with indices
    {trnek.ineLo..trnek.ineHi}. The interval {trnek.tLo .. trnek.tHi}
    must be contained in the overal time range {tLo..tHi}. */

nmsim_elem_net_trace_t *nmsim_elem_net_trace_new(nmsim_elem_neuron_count_t ntr);
  /* Allocates an {nmsim_elem_net_trace_t} record {etrace} with space 
    for {ntr} neuron traces. The neuron trace pointers 
    {etrace.trne[0..ntr-1]} are initially set to {NULL},
    and the time interval {etrace.tLo..etrace.tHi} is set
    to an empty interval. */
  
void nmsim_elem_net_trace_set
  ( nmsim_elem_net_trace_t *etrace, 
    int32_t ktr,
    nmsim_elem_neuron_trace_t *trne
  );
  /* Stores the trace record {trne} in entry {etrace.trne[ktr]}. Also updates
    the time range {etrace.tLo..etrace.tHi} to include {trne.tLo..trne.tHi}. 
    This range should prefereably have been clipped to the actual simulation 
    time range. */

void nmsim_elem_net_trace_free(nmsim_elem_net_trace_t *etrace);
  /* Releases the storage of record {etrace} including all sub-tables,
    and any neuron traces in {etrace.trne} that are not {NULL}. */ 

/* TRACE SAVING PROCEDURES 
  
  These procedures receive vectors of neuron state variables {V,age,M,H}
  for time {t} or neuron evolution variables {X,I,J} for the time step
  from {t} to {t+1}. Each vector is assumed to be indexed by neuron
  number {ine} in {0..nne-1}.
  
  The procedures store those values in the corresponding fields of each 
  applicable non-null trace record {trnek} in {etrace.trne[0..etrace.ntr-1]}.  
  See {nmsim_elem_neuron_trace_set_V_age_M_H} and 
  {}.  
  
  If trace record {trnek} referrs to only one neuron {ine}, each
  variable {V[ine],age[ine]}, etc is stored in the proper entry of
  {trnek}. If {trnek} refers to two or more neurons, the value stored is
  the average of the variable over those neurons. */

void nmsim_elem_net_trace_set_V_age_M_H
  ( nmsim_elem_net_trace_t *etrace,
    nmsim_time_t t,                /* Discrete time. */
    nmsim_elem_neuron_count_t nne, /* Number of neurons in net. */
    double V[],                    /* Membrane potential of each neuron at {t}. */
    nmsim_step_count_t age[],      /* Firing age of each neuron at time {t}. */
    double M[],                    /* Recharge modulator of each neuron at time {t}. */
    double H[]                     /* Output modulator of each neuron at time {t}. */
  );
  /* For each non-null relevant trace record {trnek} in
    {etrace->trne[...]}, stores in the proper entry of {trnek->ts} the
    values of {V[ine],age[ine],M[ine],H[ine]}, averaged for all neurons
    {ine} in {trnek.ineLo..trnek.ineHi}. See
    {nmsim_elem_neuron_trace_set_V_age_M_H}. */

void nmsim_elem_net_trace_set_X_I_J
  ( nmsim_elem_net_trace_t *etrace,
    nmsim_time_t t,                /* Time at start of step. */
    nmsim_elem_neuron_count_t nne, /* Number of neurons in net. */
    bool_t X[],                    /* Firing indicators between {t} and {t+1}. */
    double I[],                    /* External input of each neuron from {t} to {t+1} (mV). */
    double J[]                     /* Total input of each neuron from {t} to {t+1} (mV). */
  );
  /* For each non-null relevant trace record {trnek} in
    {etrace->trne[...]}, stores in the proper entry of {trnek->ts} the
    values of {X[ine],I[ine],J[ine]}, averaged for all neurons
    {ine} in {trnek.ineLo..trnek.ineHi}. See
    {nmsim_elem_neuron_trace_set_X_I_J}. */

/* INPUT-OUTPUT */

void nmsim_elem_net_trace_write(char *prefix, nmsim_elem_net_trace_t *etrace);
  /* Writes to disk each neuron trace {trnek=etrace.trne[k]}, where
    {k} is in {0..etrace.nne-1}. Each neuron trace is written as a
    separate file whose name is "{prefix}_ne{ILO}--{IHI}_trace.txt", where {ILO} and {IHI}
    are the neuron indices {trnek.ineLo,trnek.ineHi}, formatted as 10 decimal digits, zero
    padded.
    
    Each file will have the format described by
    {nmsim_elem_neuron_trace_read_INFO}. See {nmsim_elem_neuron_trace_write}.
    
    Also writes a statistical summary of each neuron trace 
    {trnek} to a separate file whose name is "{prefix}_ne{ILO}--{IHI}_stats.txt".
    See {nmsim_elem_neuron_trace_write_stats}. */

/* TESTING AND DEBUGGING */

nmsim_elem_net_trace_t *nmsim_elem_net_trace_throw
  ( nmsim_elem_neuron_count_t nne,  /* Number of neurons in network. */
    nmsim_time_t tLo,               /* Discrete time of first neuron trace entries. */
    nmsim_time_t tHi,               /* Discrete time of last neuron trace entries. */
    nmsim_elem_neuron_count_t ntr   /* Number of neuron sets to monitor. */
  );
  /* Creates a {nmsim_elem_net_trace_t} structure that monitors {ntr}
    distinct neurons with indices in {0..nne-1}.  Their time ranges
    will be roughly contained in {tLo..tHi}, but may extrapolate
    somewhat from that range. */

#endif
