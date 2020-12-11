#ifndef nmsim_elem_net_throw_H
#define nmsim_elem_net_throw_H
 
/* Generates a random neuron-level model of a network of Galves-Löcherbach neurons. */
/* Last edited on 2020-12-10 18:14:07 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_class_net.h>
#include <nmsim_group_net.h>
#include <nmsim_elem_net.h>
  
nmsim_elem_net_t *nmsim_elem_net_throw(nmsim_group_net_t *gnet);
  /* Allocates an {nmsim_elem_net_t} structure with the group-level
    network {gnet}, and generates individual neurons and synapses
    as described in the neuron and synapse group attributes. */

#endif

