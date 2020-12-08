#ifndef nmsim_elem_neuron_trace_stats_H
#define nmsim_elem_neuron_trace_stats_H
 
/* Statistical summary of a single neuron trace in a Galves-Löecherbach net. */
/* Last edited on 2020-12-07 00:43:50 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>
#include <nmsim_stats.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_neuron_trace.h>

typedef struct nmsim_elem_neuron_trace_stats_t
  { nmsim_elem_neuron_ix_t ine;           /* Index of the neuron in the network. */
    nmsim_time_t tlo;                     /* Discrete time of first trace entry. */
    nmsim_time_t thi;                     /* Discrete time of last trace entry. */
    
    nmsim_stats_t V;                      /* Stats of potential. */
    nmsim_stats_t V_fire;                 /* Stats of potential just before firing. */
    nmsim_stats_t age;                    /* Stats of age. */
    nmsim_stats_t age_fire;               /* Stats of age just before firing. */
    nmsim_stats_t M;                      /* Stats of potential decay modulator. */
    nmsim_stats_t H;                      /* Stats of output synapse strength modulator. */
    nmsim_stats_t X;                      /* Stats of firing indicator. */
    nmsim_stats_t I;                      /* Stats of external input. */
    nmsim_stats_t J;                      /* Stats of total input. */
  } nmsim_elem_neuron_trace_stats_t;
  /* A data structure that records the trace of a selected neuron with
    index {ine} for the discrete times {tlo..thi} inclusive, in
    an elem-level simulation of a GL network.
    
    The state of the neuron at discrete time {t} and the evolution
    data from {t} to {t+1} are to be stored in {ts[t-tlo]}. Note that
    the array {ts} must be allocated with at least {thi-tlo+1}
    entries. */
   
void nmsim_elem_neuron_trace_stats_compute
  ( nmsim_elem_neuron_trace_t *trne,
    nmsim_elem_neuron_trace_stats_t *trS
  );
  /* Stores into {*trS} a statistical summary of the neuron trace {trne}.  Ignores
    any values that are {±INF} or {NAN}. */

void nmsim_elem_neuron_trace_stats_write(FILE *wr, nmsim_elem_neuron_trace_stats_t *trS);
  /* Writes the statistical summary {trS} of a neuron trace to file {wr}. */
   
#define nmsim_elem_neuron_trace_stats_FILE_TYPE "nmsim_elem_neuron_trace_stats"
    
#define nmsim_elem_neuron_trace_stats_VERSION "2020-12-07"

#endif
