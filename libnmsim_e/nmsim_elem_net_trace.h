#ifndef nmsim_elem_net_trace_H
#define nmsim_elem_net_trace_H
 
/* Trace of an elem-level simulation of a Galves-Löcherbach neuron net. */
/* Last edited on 2020-12-07 16:01:37 by jstolfi */

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
  { nmsim_time_t tlo;                   /* Discrete time of first neuron trace entries. */
    nmsim_time_t thi;                   /* Discrete time of last neuron trace entries. */
    nmsim_elem_neuron_count_t tne;      /* Number of monitored neurons. */
    nmsim_elem_neuron_trace_t **trne;   /* Traces of monitored neurons. */
  } nmsim_elem_net_trace_t;
  /* A data structure that records the trace of selected neurons for a certain 
    time interval, in an elem-level simulation of a GL network. 
    
    The traces monitored neurons are stored in {trne[0..tne-1]}. For
    each {k} in {0..tne-1}, the trace {trnek=trne[k]} records the
    state summaries and evolution data of the neuron with index
    {trnek.ine}.
    
    The discrete times to monitor are {tlo..thi}.  However, each neuron trace
    {trnek} has its own time range {trnek.tlo .. trnek.thi}, which must be a
    subset of the overall range {tlo..thi}. */

nmsim_elem_net_trace_t *nmsim_elem_net_trace_new 
  ( nmsim_time_t tlo,             /* Discrete time of first neuron trace entries. */
    nmsim_time_t thi,             /* Discrete time of last neuron trace entries. */
    nmsim_elem_neuron_count_t tne  /* Number of neurons to monitor. */
  );
  /* Allocates an {nmsim_elem_net_trace_t} record {etrace} for the
    time range {tlo..thi} with space for {tne} neuron traces. The
    neuron trace pointers vector {etrace.trne[0..tne-1]} are initially
    set to {NULL}. */
  
void nmsim_elem_net_trace_free(nmsim_elem_net_trace_t *etrace);
  /* Releases the storage of record {etrace} including all sub-tables,
    and any neuron traces in {etrace.trne} that are not {NULL}. */ 

/* TRACE SAVING PROCEDURES 
  
  These procedures receive vectors {V,X,age,M,H,I,J} indexed by
  neuron number {ine} in {0..nne-1}, for time {t}
  or the time step from {t} to {t+1}.  They store those values
  in the corresponding fields of the traces of monitored neurons.
  
  More precisely, for each non-null neuron trace {trnek = etrace.trne[k]}
  with {t} in {trnek.tlo .. trnek.thi}, the procedures store in 
  neuron trace entry {tst=trnek->ts[t - trnek.tlo]} the value of {V[trnek.ine]}
  in the {.V} field, {X[trnek.ine]} in the {.X} field, etc.. */

void nmsim_elem_net_trace_set_V_age_M_H
  ( nmsim_elem_net_trace_t *etrace,
    nmsim_time_t t,                /* Discrete time. */
    nmsim_elem_neuron_count_t nne, /* Number of neurons in net. */
    double V[],                    /* Membrane potential of each neuron at {t}. */
    nmsim_step_count_t age[],      /* Firing age of each neuron at time {t}. */
    double M[],                    /* Recharge modulator of each neuron at time {t}. */
    double H[]                     /* Output modulator of each neuron at time {t}. */
  );

void nmsim_elem_net_trace_set_X_I_J
  ( nmsim_elem_net_trace_t *etrace,
    nmsim_time_t t,                /* Time at start of step. */
    nmsim_elem_neuron_count_t nne, /* Number of neurons in net. */
    bool_t X[],                    /* Firing indicators between {t} and {t+1}. */
    double I[],                    /* External input of each neuron from {t} to {t+1} (mV). */
    double J[]                     /* Total input of each neuron from {t} to {t+1} (mV). */
  );

/* INPUT-OUTPUT */

void nmsim_elem_net_trace_write(char *prefix, nmsim_elem_net_trace_t *etrace);
  /* Writes to disk each neuron trace {trnek=etrace.trne[k]}, where
    {k} is in {0..etrace.nne-1}. Each neuron trace is written as a
    separate file whose name is "{prefix}_n{NN}_trace.txt", where {NN} is
    the neuron index {trnek.ine}, formatted as 10 decimal digits, zero
    padded.
    
    Each file will have the format described by
    {nmsim_elem_neuron_trace_read_INFO}. See {nmsim_elem_neuron_trace_write}.
    
    Also writes a statistical summary of each neuron trace 
    {trnek} to a separate file whose name is "{prefix}_n{NN}_stats.txt".
    See {nmsim_elem_neuron_trace_write_stats}. */

/* TESTING AND DEBUGGING */

nmsim_elem_net_trace_t *nmsim_elem_net_trace_throw
  ( nmsim_elem_neuron_count_t nne,  /* Number of neurons in network. */
    nmsim_time_t tlo,               /* Discrete time of first neuron trace entries. */
    nmsim_time_t thi,               /* Discrete time of last neuron trace entries. */
    nmsim_elem_neuron_count_t tne   /* Number of neurons to monitor. */
  );
  /* Creates a {nmsim_elem_net_trace_t} structure that monitors {tne}
    distinct neurons with indices in {0..nne-1}.  Their time ranges
    will be roughly contained in {tlo..thi}, but may extrapolate
    somewhat from that range. */

#endif
