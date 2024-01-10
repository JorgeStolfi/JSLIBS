#ifndef nmsim_group_net_trace_H
#define nmsim_group_net_trace_H
 
/* Trace of an elem-level simulation of a Galves-Löcherbach neuron net. */
/* Last edited on 2019-02-13 17:43:10 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_neuron_trace.h>

typedef struct nmsim_group_net_trace_t
  { nmsim_time_t tini;                    /* Discrete time of first recorded states. */
    nmsim_time_t tfin;                    /* Discrete time of last recorded states. */
    nmsim_group_neuron_count_t nng;       /* Number of neuron groups to monitor. */
    nmsim_group_neuron_trace_t **trng;      /* Traces of monitored neuron groups. */
  } nmsim_group_net_trace_t;
  /* A data structure that records the trace of {nng} selected neuron groups for a certain 
    time interval, in a group-level simulation of a GL netowrk. 
    
    The indices of the neuron groups to monitor are assumed to be in {ing[0..nng-1]}.
    
    The discrete times to monitor are {tini..tfin}.
    
    The states of the neuron group with index {ing[i]} are recorded in
    {trng[i]}. The time range {trng[i].ini .. trng[i].tfin} must be a subset
    of the range {tini..tfin}. */

nmsim_group_net_trace_t *nmsim_group_net_trace_new 
  ( nmsim_time_t tini,                   /* Discrete time of first recorded states. */
    nmsim_time_t tfin,                   /* Discrete time of last recorded states. */
    nmsim_group_neuron_count_t nng_tr    /* Number of neuron groups to monitor. */
  );
  /* Allocates an {nmsim_group_net_trace_t} record {gtrace} for the time range 
    {tini..tfin} with a vector {gtrace.trng} with space for {nng_tr} 
    neuron group traces.  The elements {gtrace.trng[i]} are all set to {NULL};
    they must be allocated and set by the client 
    (see {nmsim_group_neuron_trace_new}). */
  
void nmsim_group_net_trace_free(nmsim_group_net_trace_t *gtrace);
  /* Releases the storage of record {gtrace} including all sub-tables,
    and any neuron group traces in {gtrace.trng} that are not {NULL}. */ 

/* INPUT-OUTPUT */

void nmsim_group_net_trace_write(char *prefix, nmsim_group_net_trace_t *gtrace);
  /* Writes each neuron group trace {gtrace.trng[i]} to disk, as {gtrace.nng} files whose names
    are "{prefix}_g{NN}.txt", where {NN} is the neuron group index {gtrace.trAXS[i].ing},
    formatted as 10 decimal digits, zero padded.
    
    Each file will have the format described by
    {nmsim_group_neuron_trace_read_INFO}. */

#endif
