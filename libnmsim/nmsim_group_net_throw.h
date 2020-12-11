#ifndef nmsim_group_net_throw_H
#define nmsim_group_net_throw_H
 
/* Generates random population-level model of a network of Galves-Löcherbach neurons. */
/* Last edited on 2020-12-11 06:26:36 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <nmsim_basic.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_synapse.h>
#include <nmsim_class_net.h>
#include <nmsim_group_net.h>

nmsim_group_net_t *nmsim_group_net_throw
  ( nmsim_class_net_t *cnet,
    nmsim_group_neuron_count_t nng, 
    nmsim_group_synapse_count_t nsg,
    nmsim_elem_neuron_count_t nne, 
    nmsim_elem_synapse_count_t nse
  );
  /* Creates a random group-level network description with the 
    class-level description {cnet}, with {nng} neuron neuron populations
    and {nsg} synaptic bundles between them, assuming a total of {nne}
    neurons and {nse} synapses. */

#endif

