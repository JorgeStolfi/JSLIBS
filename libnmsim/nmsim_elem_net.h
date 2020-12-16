#ifndef nmsim_elem_net_H
#define nmsim_elem_net_H
 
/* Neuron-level model of a network of Galves-Löcherbach neurons. */
/* Last edited on 2020-12-13 16:03:25 by jstolfi */

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
    nmsim_elem_neuron_t *neu;  /* List of neurons. */
    /* Indexed with synapse index {0 .. nse-1}: */
    nmsim_elem_synapse_t *syn; /* List of synapses. */
  } nmsim_elem_net_t;
  /* An object {enet} of this type is an element-level description 
    of a network of Galves-Loecherbach neurons. See {nmsim_elem_net_INFO}.
    
    For efficient processing, the neurons in memory are sorted by increasing.
    group index. The irst neuron of group {ing} has index {gnet->ngrp[ing].ine_start}.
    
    For efficient simulation, the synapses in memory are sorted by the
    index of the pre-synaptic neuron. Thus all output synapses of the
    same neuron {ine} have consecutive indices, between those of neurons
    {ine-1} and {ine+1}. The field {.ise_out_start} of each neuron is the sum
    of the {.nse_out} fields of all previous neurons (even if the neuron 
    has zero outgoing synapses).
    
    Thus the output synapses of neuron {ine} are
    {enet->syn[ini..fin]} where 
    
      {ini = enet->neu[ine].ise_out_start} 
      {fin = ini + enet->neu[ine].nse_out - 1}
      
    */
  
#define nmsim_elem_net_INFO \
  "The element-level decription comprises a group-level" \
  " description {gnet} of the network, a list of individual" \
  " neuron descriptions {neu[0..nne-1]}, and" \
  " a list of individual synapses {syn[0..nse-1]}."
    
/* NETWORK CREATION */

nmsim_elem_net_t *nmsim_elem_net_new
  ( nmsim_group_net_t *gnet,
    nmsim_elem_neuron_count_t nne,
    nmsim_elem_synapse_count_t nse
  );
  /* Allocates an {nmsim_elem_net_t} structure {enet} with the group-level
    network {gnet}, with space for {nne} neurons and {nse} synapses.
    
    The tables {enet.neu} and {enet.syn} are fully allocated,
    their contents is undefined.  Thus the description is invalid, 
    until those tables has been properly filled. */
    
void nmsim_elem_net_free(nmsim_elem_net_t *enet);
  /* Releases all storage of the element-level network {enet},
    including internal tables and the top-level
    descriptor itself. */

/* INPUT/OUTPUT */

void nmsim_elem_net_write(FILE *wr, nmsim_elem_net_t *enet, double timeStep);
  /* Writes a description of the element-level network net {enet} to file {wr},
    in the format described by {nmsim_elem_net_read_INFO} below.
    
    The {timeStep} parameter (in milliseconds) is used to convert
    neuron decay factors to step-invariant characteristic decay
    times. */

nmsim_elem_net_t *nmsim_elem_net_read(FILE *rd, double timeStep);
  /* Reads an element-level network description from {rd}, in the format
    described by {nmsim_elem_net_read_INFO} below. 
    
    The {timeStep} parameter (in milliseconds) is used to convert
    step-invariant characteristic decay times of the neurons into the
    step-dependent decay factors. 
    
    The synapses are reordered and re-indexed if necessary 
    so that the pre-synaptic neuron indices are not decreasing.
    The field {.nse_out} of each neuron is checked, and its 
    {.ise_out_start} field is set appropriately. */
    
#define nmsim_elem_net_read_INFO \
  "The element-level description of the network begins with a line \n" \
  "    begin " nmsim_elem_net_FILE_TYPE " (format of " nmsim_elem_net_VERSION ")\n" \
  " followed by:\n" \
  "      * a line \"neuron_elems = {nne}\" where {nne} is the number of neurons\n" \
  "      * a line \"synapse_elems = {nse}\" where {nse} is the number of synapses\n" \
  "      * the group-level description of the network (see below)\n" \
  "      * a set of {nne} lines that describe the neurons (see below)\n" \
  "      * a set of {nse} lines that describe synapse bundles (ditto)\n" \
  "  The description ends with the line \"end " nmsim_elem_net_FILE_TYPE "\".\n" \
  "\n" \
  "  The neurons must be sorted by increasing group indices. The synapses must be sorted" \
  " by increasing index of the pre-synaptic neuron.\n" \
  "\n" \
  "GROUP-LEVEL NETWORK DESCRIPTION\n" \
  "  " nmsim_group_net_read_INFO "\n" \
  "\n" \
  "DESCRIPTION OF A NEURON\n" \
  "  " nmsim_elem_neuron_read_INFO "\n" \
  "\n" \
  "DESCRIPTION OF A SYNAPSE\n" \
  "  " nmsim_elem_synapse_read_INFO ""
    
#define nmsim_elem_net_FILE_TYPE "nmsim_elem_net"
    
#define nmsim_elem_net_VERSION "2020-12-10"

/* TESTING AND DEBUGGING */

void nmsim_elem_net_compare(nmsim_elem_net_t *enet_read, nmsim_elem_net_t *enet_orig);
  /* Compares the elem-level network description {enet_read} read from a file
    with the expected description {enet_orig}.  Also compares the group- and class-level
    network descriptions.  Aborts if there are discrepancies (beyond
    roundoff tolerance). */

#endif

