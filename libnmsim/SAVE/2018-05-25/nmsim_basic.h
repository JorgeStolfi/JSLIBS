#ifndef nmsim_basic_H
#define nmsim_basic_H
 
/* Basic types for neuromat network simulation. */
/* Last edited on 2018-03-24 05:02:50 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>

/* COUNTS AND INDICES FOR NEURON-LEVEL MODELING */

typedef int64_t nmsim_neuron_count_t;
  /* Count of neurons. */
  
typedef int64_t nmsim_neuron_ix_t;
  /* Index of a neuron in a network. */
  
typedef int64_t nmsim_synapse_count_t;
  /* Count of synapses in a neuron or network. */
  
typedef int64_t nmsim_synapse_ix_t;
  /* Index of a synapse in a neuron or network. */

/* COUNTS AND INDICES FOR POPULATION-LEVEL MODELING */

typedef int32_t nmsim_pop_count_t;
  /* Count of neuron populations. */

typedef int32_t nmsim_pop_ix_t;
  /* Index of a neuron population. */
  
typedef int32_t nmsim_bundle_count_t;
  /* Count of synaptic bundles populations. */

typedef int32_t nmsim_bundle_ix_t;
  /* Index of a synaptic bundle. */

#endif
