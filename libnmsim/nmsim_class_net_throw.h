#ifndef nmsim_class_net_throw_H
#define nmsim_class_net_throw_H
 
/* Generates a random class-level description of a network of Galves-Löcherbach neurons. */
/* Last edited on 2020-12-11 13:36:46 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_class_net.h>

nmsim_class_net_t *nmsim_class_net_throw
  ( nmsim_class_neuron_count_t nnc, 
    nmsim_class_synapse_count_t nsc
  );
  /* Creates a random class-level network description with {nnc} neuron classes
    (at least 1) and {nsc} synapse classes. */

#endif

