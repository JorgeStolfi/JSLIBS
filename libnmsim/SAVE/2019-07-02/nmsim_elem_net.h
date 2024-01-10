#ifndef nmsim_elem_net_H
#define nmsim_elem_net_H
 
/* Neuron-level model of a network of Galves-Löcherbach neurons. */
/* Last edited on 2019-06-18 11:18:49 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_synapse.h>
#include <nmsim_elem_synapse.h>

#include <nmsim_class_net.h>
#include <nmsim_group_net.h>

#include <nmsim_elem_net.h>
  
typedef struct nmsim_elem_net_t
  { nmsim_group_net_t *gnet;  /* Neuron and synapse groups. */
    /* Neuron and synapse counts: */
    nmsim_elem_neuron_count_t nne;       /* Total number of neurons. */
    nmsim_elem_synapse_count_t nse;      /* Total number of synapses. */
    /* Indexed with neuron index {0 .. nne-1}: */
    nmsim_elem_neuron_vec_t neu;       /* List of neurons. */
    /* Indexed with synapse index {0 .. nse-1}: */
    nmsim_elem_synapse_vec_t syn;      /* List of synapses. */
  } nmsim_elem_net_t;
  /* Element-level description of a network of Galves-Loecherbach neurons.
    See {nmsim_elem_net_INFO}. */
  
#define nmsim_elem_net_INFO \
  "The element-level decription comprises a group-level" \
  " description {gnet} of the network, a list of individual" \
  " neuron descriptions {neu.e[0..nne-1]}, and" \
  " a list of individual synapses {syn.e[0..nse-1]}."
    
/* NETWORK CREATION */

nmsim_elem_net_t *nmsim_elem_net_throw(nmsim_group_net_t *gnet);
  /* Allocates an {nmsim_elem_net_t} structure with the group-level
    network {gnet}, and generates individual neurons and synapses
    as described in the neuron and synapse group attributes. */
    
void nmsim_elem_net_free(nmsim_elem_net_t *enet);
  /* Releases all storage of the element-level network {enet},
    including internal tables and the top-level
    descriptor itself. */

/* INPUT/OUTPUT */

void nmsim_elem_net_write(FILE *wr, nmsim_elem_net_t *enet, double timeStep);
  /* Writes a description of the network net {enet} to file {wr},
    in the format described by {nmsim_elem_net_read_INFO} below. 
    
    The {timeStep} parameter (in milliseconds) is used to convert
    neuron decay factors to step-invariant characteristic decay
    times. */

nmsim_elem_net_t *nmsim_elem_net_read(FILE *rd, double timeStep);
  /* Reads a group_level network description from {rd}, in the format described 
    by {nmsim_elem_net_read_INFO} below. 
    
    The {timeStep} parameter (in milliseconds) is used to convert
    step-invariant characteristic decay times of the neurons into the
    step-dependent decay factors. */
    
#define nmsim_elem_net_read_INFO \
  "The element-level description of the network begins with a line \n" \
  "    begin " nmsim_elem_net_FILE_TYPE " (format of " nmsim_elem_net_VERSION ")\n" \
  " followed by the group-level description of the network (see below), and two blocks of lines that describe the individual neurons and synapses, and" \
  " a closing line \"end " nmsim_elem_net_FILE_TYPE "\".\n" \
  "\n" \
  "  The neuron list is preceded by a line \"neuron_elems = {nne}\" where {nne} is the number of neurons.  Then follow {nne} lines with the attributes of neurons (see below).\n" \
  "\n" \
  "  The synapse list is likewise preceded by a line \"synapse_elems = {nse}\", where {nse} is the number of synapses.  There follow {nse} lines with the attributes of the synapses (see below).\n" \
  "\n" \
  "GROUP-LEVEL DESCRIPTION\n" \
  "  " nmsim_group_net_read_INFO "\n" \
  "\n" \
  "NEURON DESCRIPTION\n" \
  "  " nmsim_elem_neuron_read_INFO "\n" \
  "\n" \
  "SYNAPSE DESCRIPTION\n" \
  "  " nmsim_elem_synapse_read_INFO ""
    
#define nmsim_elem_net_FILE_TYPE "nmsim_elem_net"
    
#define nmsim_elem_net_VERSION "2019-01-11"

/* TESTING AND DEBUGGING */

void nmsim_elem_net_compare(nmsim_elem_net_t *enet_read, nmsim_elem_net_t *enet_orig);
  /* Compares the elem-level network description {enet_read} read from a file
    with the expected description {enet_orig}.  Also compares the group- and class-level
    network descriptions.  Aborts if there are discrepancies (beyond
    roundoff tolerance). */

#endif

