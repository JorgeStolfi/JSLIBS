#ifndef nmsim_elem_neuron_trace_stats_H
#define nmsim_elem_neuron_trace_stats_H
 
/* Statistical summary of a single neuron trace in a Galves-Löecherbach net. */
/* Last edited on 2020-12-15 22:23:05 by jstolfi */

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
#include <nmsim_elem_net_sim_stats.h>
#include <nmsim_elem_neuron_trace.h>

void nmsim_elem_neuron_trace_stats_compute
  ( nmsim_elem_neuron_trace_t *trne,
    nmsim_elem_net_sim_stats_t *ntS
  );
  /* Stores into {*ntS} a statistical summary of the neuron trace {trne}.  Ignores
    any values that are {±INF} or {NAN}. */

#endif
