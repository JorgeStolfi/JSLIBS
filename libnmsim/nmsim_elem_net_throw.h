#ifndef nmsim_elem_net_throw_H
#define nmsim_elem_net_throw_H
 
/* Generates a random neuron-level model of a network of Galves-Löcherbach neurons. */
/* Last edited on 2021-01-07 21:12:36 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_basic.h>

#include <nmsim_class_net.h>
#include <nmsim_group_net.h>
#include <nmsim_elem_net.h>
  
nmsim_elem_net_t *nmsim_elem_net_throw(nmsim_group_net_t *gnet);
  /* Allocates an {nmsim_elem_net_t} structure with the group-level
    network {gnet}, and generates individual neurons and synapses
    as described in the neuron and synapse group attributes. */

#endif

