#ifndef nmsim_elem_neuron_trace_H
#define nmsim_elem_neuron_trace_H
 
/* Trace of a single neuron in an elem-level simulation of a Galves-Löcherbach net. */
/* Last edited on 2019-02-13 17:41:56 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_neuron_state.h>

typedef struct nmsim_elem_neuron_trace_t
  { nmsim_elem_neuron_ix_t ine;       /* Index of the neuron in the network. */
    nmsim_group_neuron_ix_t ing;      /* Index of the neuron's group in the network. */
    nmsim_class_neuron_ix_t inc;      /* Index of the neuron's class in the network. */
    nmsim_time_t tini;                /* Discrete time of first recorded state. */
    nmsim_time_t tfin;                /* Discrete time of last recorded state. */
    nmsim_elem_neuron_state_t *st;    /* States of the neuron. */
  } nmsim_elem_neuron_trace_t;
  /* A data structure that records the trace of a selected neuron for a certain 
    time interval, in an elem-level simulation of a GL netowrk. 
    
    The fields {ine,ing,inc} are provided for documentation, but are not
    otherwise used.
    
    The discrete times considered are {tini..tfin}.
    
    The state of the neuron at discrete time {t} is to be stored in 
    {st[t-tini]}. Note that the array {st} must be allocated with at least
    {tfin-tini+1} entries. */
    
nmsim_elem_neuron_trace_t *nmsim_elem_neuron_trace_new
  ( nmsim_elem_neuron_ix_t ine,       /* Index of the neuron in the network. */
    nmsim_group_neuron_ix_t ing,      /* Index of the neuron's group in the network. */
    nmsim_class_neuron_ix_t inc,      /* Index of the neuron's class in the network. */
    nmsim_time_t tini,                /* Discrete time of first recorded state. */
    nmsim_time_t tfin                 /* Discrete time of last recorded state. */
  );
  /* Allocates a new {nmsim_elem_neuron_trace_t} record on the heap
    with space for states at times {tini..tfin}. */
    
#define nmsim_elem_neuron_trace_states_MAX (((int64_t)1) << 24)
  /* Max number of recorded states, for memory allocation purposes. */

void nmsim_elem_neuron_trace_free(nmsim_elem_neuron_trace_t *trne);
  /* Releases the storage of the record {trne}, including internal tables. */

void nmsim_elem_neuron_trace_write(FILE *wr, nmsim_elem_neuron_trace_t *trne);
  /* Writes the neuron trace {trne} to file {wr}, in the format described by 
    {nmsim_elem_neuron_trace_read_INFO} below. */
    
nmsim_elem_neuron_trace_t *nmsim_elem_neuron_trace_read(FILE *rd);
  /* Reads a neuron trace {trne} from file {rd}, in the format described by 
    {nmsim_elem_neuron_trace_read_INFO} below. */
    
#define nmsim_elem_neuron_trace_read_INFO \
  "  The file begins with a line \n" \
  "      \"begin " nmsim_elem_neuron_trace_FILE_TYPE " format of " nmsim_elem_neuron_trace_VERSION "\"\n" \
  "  Then follow a set of lines in the format {NAME} = {VALUE}, where {NAME} is \"neuron_elem\", " \
  " \"neuron_group\", \"neuron_class\", \"initial_time\", and \"final time\". Then" \
  " follow one line for each discrete time in the monitored interval.\n" \
  "\n" \
  "  " nmsim_elem_neuron_state_read_INFO ""

#define nmsim_elem_neuron_trace_FILE_TYPE "nmsim_elem_neuron_trace"
    
#define nmsim_elem_neuron_trace_VERSION "2019-01-31"

#endif
