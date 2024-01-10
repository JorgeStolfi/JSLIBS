#ifndef nmsim_elem_net_trace_H
#define nmsim_elem_net_trace_H
 
/* Trace of an elem-level simulation of a Galves-Löcherbach neuron net. */
/* Last edited on 2019-02-13 18:29:15 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_neuron_trace.h>

typedef struct nmsim_elem_net_trace_t
  { nmsim_time_t tini;                   /* Discrete time of first recorded states. */
    nmsim_time_t tfin;                   /* Discrete time of last recorded states. */
    nmsim_elem_neuron_count_t nne;       /* Number of neurons to monitor. */
    nmsim_elem_neuron_trace_t **trne;      /* Traces of monitored neurons. */
  } nmsim_elem_net_trace_t;
  /* A data structure that records the trace of {nne} selected neurons for a certain 
    time interval, in an elem-level simulation of a GL netowrk. 
    
    The indices of the neurons to monitor are assumed to be in {ine[0..nne-1]}.
    
    The discrete times to monitor are {tini..tfin}.
    
    The states of the neuron with index {ine[i]} are 
    recorded in {trne[i]}.  The time range 
    {trne[i].ini .. trne[i].tfin} must be a subset of 
    the range {tini..tfin}. */

nmsim_elem_net_trace_t *nmsim_elem_net_trace_new 
  ( nmsim_time_t tini,                   /* Discrete time of first recorded states. */
    nmsim_time_t tfin,                   /* Discrete time of last recorded states. */
    nmsim_elem_neuron_count_t nne_tr     /* Number of neurons to monitor. */
  );
  /* Allocates an {nmsim_elem_net_trace_t} record {etrace} for the time range 
    {tini..tfin} with a vector {etrace.trne} with space for {nne_tr} 
    neuron traces.  The elements {etrace.trne[i]} are all set to {NULL};
    they must be allocated and set by the client 
    (see {nmsim_elem_neuron_trace_new}). */
  
void nmsim_elem_net_trace_free(nmsim_elem_net_trace_t *etrace);
  /* Releases the storage of record {etrace} including all sub-tables,
    and any neuron traces in {etrace.trne} that are not {NULL}. */ 

/* TRACE SAVING PROCEDURES 
  
  These procedures receive vectors {V,X,age,H,M,I,J} indexed by
  neuron number {ine} in {0..nne-1}, for time {t}
  or the time step from {t} to {t+1}.  They store those values
  in the corresponding fields of the traces of monitored neurons.
  
  More precisely, for each non-null neuron trace {trne = etrace.trne[i]}
  with {t} in {trne.tini .. trne.tfin}, the procedures store in 
  state record {st = trne->st[t - trne.tini]} the value of {V[trne.ine]}
  in the {.V} field, {X[trne.ine]} in the {.X} field, etc.. */

void nmsim_elem_net_trace_set_V_X
  ( nmsim_elem_net_trace_t *etrace,
    nmsim_time_t t,                /* Time at start of step. */
    nmsim_elem_neuron_count_t nne, /* Number of neurons in net. */
    double V[],                    /* Membrane potential of each neuron at {t}. */
    bool_t X[]                     /* Firing indicators between {t} and {t+1}. */
  );

void nmsim_elem_net_trace_set_age_H_I_J
  ( nmsim_elem_net_trace_t *etrace,
    nmsim_time_t t,                /* Time at start of step. */
    nmsim_elem_neuron_count_t nne, /* Number of neurons in net. */
    nmsim_step_count_t age[],      /* Firing age of each neuron at time {t}. */
    double H[],                    /* Output modulator of each neuron at time {t}. */
    double I[],                    /* External input of each neuron from {t} to {t+1}. */
    double J[]                     /* Total input of each neuron from {t} to {t+1}. */
  );

void nmsim_elem_net_trace_set_M
  ( nmsim_elem_net_trace_t *etrace,
    nmsim_time_t t,                /* Time at start of step. */
    nmsim_elem_neuron_count_t nne, /* Number of neurons in net. */
    double M[]                     /* Recharge modulator of each neuron at time {t}. */
  );

/* INPUT-OUTPUT */

void nmsim_elem_net_trace_write(char *prefix, nmsim_elem_net_trace_t *etrace);
  /* Writes each neuron trace {etrace.trne[i]} to disk, as {etrace.nne} files whose names
    are "{prefix}_n{NN}.txt", where {NN} is the neuron index {etrace.trAXS[i].ine},
    formatted as 10 decimal digits, zero padded.
    
    Each file will have the format described by
    {nmsim_elem_neuron_trace_read_INFO}. */

#endif
