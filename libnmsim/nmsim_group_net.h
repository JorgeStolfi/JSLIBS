#ifndef nmsim_group_net_H
#define nmsim_group_net_H
 
/* Population-level model of a network of Galves-Löcherbach neurons. */
/* Last edited on 2020-12-11 14:53:25 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <nmsim_basic.h>
#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_synapse.h>
#include <nmsim_class_net.h>

typedef struct nmsim_group_net_t
   { nmsim_class_net_t *cnet;            /* Neuron and synapse classes.*/
    /* Neuron groups (populations): */
    nmsim_group_neuron_count_t nng;      /* Count of neuron groups (at least 1). */
    nmsim_group_neuron_t *ngrp;          /* Attributes of neuron groups, indexed {0..nng-1}. */
    /* Synapse groups (bundles): */
    nmsim_group_synapse_count_t nsg;     /* Count of synapse groups (may be zero). */
    nmsim_group_synapse_t *sgrp;         /* Attributes of synapse groups, indexed {0..nsg-1}. */
  } nmsim_group_net_t;
  /* A record {gnet}of this type is a group-level description of a 
    network of Galves-Loecherbach neurons,
    
    The neurons are grouped in {.nng} groups (neuron populations),
    belonging {cnet.nnc} distinct neuron classes. Each neuron group has
    an index {ing} in {0..nng-1}; its attributes record is {P =
    gnet.ngrp[ing]}, and additional parameters are specified 
    in {cnet.nclass[P.inc]}.
    
    The synapses are divided in {gnet.nsg} synapse groups (synaptic
    bundles), belonging to {cnet.nsc} distinct synapse classes. Each
    synapse group has an index {isg} in {0..nsg-1}; its attributes recordis {B =
    gnet.sgrp[isg]}, and additional parameters are specified in {cnet.sclass[B.inc]}. Each
    bundle connects neurons in the population with index {B.ing_pre}
    to those in the population {B.ing.pos}.  If {cnet} has no synapse
    classes ({cnet-.sc = 0}) then {gnet} cannot have any synaptic 
    bundles ({gnet.nsg = 0}).
    
    There must be at least one neuron group, and therefore at least one
    neuron class in {cnet}.  If there is at least one synapse group,
    there must be at least one synapse class in {cnet}. */
    
/* NETWORK CREATION */

nmsim_group_net_t *nmsim_group_net_new
  ( nmsim_class_net_t *cnet, 
    nmsim_group_neuron_count_t nng,
    nmsim_group_synapse_count_t nsg
  );
  /* Allocates an {nmsim_group_net_t} structure {gnet}, with the neuron and
    synapse classes defined in {cnet}.
    
    The tables {gnet.ngrp} and {gnet.sgrp} are allocated for {nng}
    neuron groups and {nsg} synapse groups, respectively. The contents
    of those tables will be undefined, so {gnet} will be invalid until
    they are properly filled. */

void nmsim_group_net_free(nmsim_group_net_t *gnet);
  /* Releases all storage of the population network {gnet},
    including the class network {gnet.cnet}, all internal tables, 
    and the top-level  descriptor itself. */

/* INPUT/OUTPUT */

void nmsim_group_net_write(FILE *wr, nmsim_group_net_t *gnet, double timeStep);
  /* Writes a description of the network net {gnet} to file {wr},
    in the format described by {nmsim_group_net_read_INFO} below. 
    
    The {timeStep} parameter (in milliseconds) is used to convert
    neuron decay factors to step-invariant characteristic decay
    times. */

nmsim_group_net_t *nmsim_group_net_read
  ( FILE *rd, 
    nmsim_elem_neuron_count_t nne_g_max,
    nmsim_elem_neuron_count_t nse_g_max,
    double timeStep
  );
  /* Reads a group_level network description from {rd}, in the format described 
    by {nmsim_group_net_read_INFO} below. 
    
    Fails if the number of neurons in any neuron group exceeds {nne_g_max},
    or the number of synapses in any synaptic bundle exceeds {nse_g_max}.
    
    The {timeStep} parameter (in milliseconds) is used to convert
    step-invariant characteristic decay times of the neurons into the
    step-dependent decay factors.  See {nmsim_class_net_read}.
    
    The procedure ensures that the synaptic bundles are sorted
    by the pre-synaptic group index.  It sets the field {.ine_start} of each
    neuron group to zero. */
    
#define nmsim_group_net_read_INFO \
  "The group-level description of the network begins with a line \n" \
  "    begin " nmsim_group_net_FILE_TYPE " (format of " nmsim_group_net_VERSION ")\n" \
  " followed by:\n" \
  "      * a line \"neuron_groups = {nng}\" where {nng} is the number of neuron groups (populations)\n" \
  "      * a line \"synapse_groups = {nsg}\" where {nsg} is the number of synaptic bundles\n" \
  "      * the class-level description of the network (see below)\n" \
  "      * a set of {nng} lines that describe the neuron poulations (see below)\n" \
  "      * a set of {nsg} lines that describe synapse bundles (ditto)\n" \
  "  The description ends with the line \"end " nmsim_group_net_FILE_TYPE "\".\n" \
  "\n" \
  "  The synaptic bundles must be sorted by increasing index of the pre-synaptic neuron group.\n" \
  "\n" \
  "CLASS-LEVEL NETWORK DESCRIPTION\n" \
  "  " nmsim_class_net_read_INFO "\n" \
  "\n" \
  "DESCRIPTION OF A NEURON POPULATION\n" \
  "  " nmsim_group_neuron_read_INFO "\n" \
  "\n" \
  "DESCRIPTION OF A SYNAPSE BUNDLE\n" \
  "  " nmsim_group_synapse_read_INFO ""
    
#define nmsim_group_net_FILE_TYPE "nmsim_group_net"
    
#define nmsim_group_net_VERSION "2020-12-10"

/* TESTING AND DEBUGGING */

void nmsim_group_net_compare(nmsim_group_net_t *gnet_read, nmsim_group_net_t *gnet_orig);
  /* Compares the group-level network description {gnet_read} read from a file
    with the expected description {gnet_orig}.  Also compares the class-level
    network descriptions.  Aborts if any parameter of any class or
    group does not match (within roundoff tolerance). */

#endif

