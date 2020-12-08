#ifndef nmsim_group_net_H
#define nmsim_group_net_H
 
/* Population-level model of a network of Galves-Löcherbach neurons. */
/* Last edited on 2020-12-07 23:37:53 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_synapse.h>
#include <nmsim_class_net.h>

typedef struct nmsim_group_net_t
   { nmsim_class_net_t *cnet;            /* Neuron and synapse classes.*/
    /* Neuron groups (populations): */
    nmsim_group_neuron_count_t nng;      /* Count of neuron groups). */
    nmsim_group_neuron_vec_t ngrp;       /* Attributes of neuron groups. */
    /* Synapse groups (bundles): */
    nmsim_group_synapse_count_t nsg;     /* Count of synapse groups. */
    nmsim_group_synapse_vec_t sgrp;      /* Attributes of synapse groups. */
  } nmsim_group_net_t;
  /* Group-level description of a network of Galves-Loecherbach neurons,
    
    The neurons are groupled in {nng} groups (neuron populations),
    belonging {cnet->nnc} distinct neuron classes. Each neuron group has
    an index {ing} in {0..nng-1}; its attributes are {P =
    ngrp.e[ing]}, and its parameters are {cnet->nclass.e[P.inc]}.
    
    The synapses are divided in {nsg} synapse groups (synaptic
    bundles), belonging to {cnet->nsc} distinct synapse classes. Each
    synapse group has an index {isg} in {0..nsg-1}; its attributes are {B =
    sgrp.e[isg]}, and its parameters are {cnet->sclass.e[B.inc]}. Each
    bundle connects neurons in the population with index {B.ing_pre}
    to those in the population {B.ing.pos}.
    
    The allocated size of the {ngrp} vector is {ngrp.ne} which is
    usually larger than then number {nsg} of used entries. The same
    applies to the {sgrp} vector.
    
    If there is at least one neuron group, there must be at least one
    neuron class in {cnet}.  If there is at least one synapse group,
    there must be at least one synapse class in {cnet}, and there 
    must be at least one neuron group for the synapses to connect to/from. */
    
/* NETWORK CREATION */

nmsim_group_net_t *nmsim_group_net_new(nmsim_class_net_t *cnet);
  /* Allocates an {nmsim_group_net_t} structure, with the neuron and
    synapse classes defined in {cnet}, initially with zero neuron
    populations and zero synaptic bundles. */

void nmsim_group_net_free(nmsim_group_net_t *gnet);
  /* Releases all storage of the population network {gnet},
    including the class network {gnet.cnet}, all internal tables, and the top-level
    descriptor itself. */

nmsim_group_neuron_ix_t nmsim_group_net_add_neuron_group
  ( nmsim_group_net_t *gnet, 
    nmsim_group_neuron_t *ngrp
  );
  /* Adds to the network {gnet} a new neuron group with attributes
    {*ngrp}. The population is stored in the {gnet->ngrp} vector, which
    is expanded if necessary, and incrementng {gnet->nng}. Returns the
    assigned population index {ing}. */
    
nmsim_group_synapse_ix_t nmsim_group_net_add_synapse_group
  ( nmsim_group_net_t *gnet,
    nmsim_group_synapse_t *sgrp
  );
  /* Adds to the network {gnet} a new synaptic bundle with attributes
    {*sgrp}. The bundle is stored in the {gnet->sgrp} vector, which is
    expanded if necessary, and incrementng {gnet->nsg}. Returns the
    assigned bundle index {isg}. */

/* INPUT/OUTPUT */

void nmsim_group_net_write(FILE *wr, nmsim_group_net_t *gnet, double timeStep);
  /* Writes a description of the network net {gnet} to file {wr},
    in the format described by {nmsim_group_net_read_INFO} below. 
    
    The {timeStep} parameter (in milliseconds) is used to convert
    neuron decay factors to step-invariant characteristic decay
    times. */

nmsim_group_net_t *nmsim_group_net_read(FILE *rd, double timeStep);
  /* Reads a group_level network description from {rd}, in the format described 
    by {nmsim_group_net_read_INFO} below. 
    
    The {timeStep} parameter (in milliseconds) is used to convert
    step-invariant characteristic decay times of the neurons into the
    step-dependent decay factors. */
    
#define nmsim_group_net_read_INFO \
  "The group-level description of the network begins with a line \n" \
  "    begin " nmsim_group_net_FILE_TYPE " (format of " nmsim_group_net_VERSION ")\n" \
  " followed by the class-level description of the network (see below), and" \
  " two blocks of lines that describe the neuron poulations and synapse bundles, and" \
  " a closing line \"end " nmsim_group_net_FILE_TYPE "\".\n" \
  "\n" \
  "  The neuron population list is preceded by a line \"neuron_groups = {nng}\" where" \
  " {nng} is the number of populations.  Then follow {nng} lines with the attributes of" \
  " neuron groups (see below).\n" \
  "\n" \
  "  The synapse bundle list is likewise preceded by a line \"synapse_groups = {nsg}\", where" \
  " {nsg} is the number of bundles.  There follow {nsg} lines with the attributes of the" \
  " synapse bundles (see below).\n" \
  "\n" \
  "CLASS-LEVEL DESCRIPTION\n" \
  "  " nmsim_class_net_read_INFO "\n" \
  "\n" \
  "NEURON POPULATION DESCRIPTION\n" \
  "  " nmsim_group_neuron_read_INFO "\n" \
  "\n" \
  "SYNAPSE BUNDLE DESCRIPTION\n" \
  "  " nmsim_group_synapse_read_INFO ""
    
#define nmsim_group_net_FILE_TYPE "nmsim_group_net"
    
#define nmsim_group_net_VERSION "2019-01-12"

/* TESTING AND DEBUGGING */

nmsim_group_net_t *nmsim_group_net_throw
  ( nmsim_class_net_t *cnet,
    nmsim_group_neuron_count_t nng, 
    nmsim_group_synapse_count_t nsg,
    nmsim_elem_neuron_count_t nne_g_min,
    nmsim_elem_neuron_count_t nne_g_max,
    double K_min,
    double K_max
  );
  /* Creates a random group-level network description with the 
    class-level description {cnet}, with {nng} neuron neuron populations
    and {nsg} synaptic bundles between them.  Each neuron group will have
    a random number of neurons in {nne_g_min .. nne_g_max}, and the average in-degree 
    {K} of the post-synaptic neurons in each synapse group will be a random
    number in {[K_min_ K_max]}. */

void nmsim_group_net_compare(nmsim_group_net_t *gnet_read, nmsim_group_net_t *gnet_orig);
  /* Compares the group-level network description {gnet_read} read from a file
    with the expected description {gnet_orig}.  Also compares the class-level
    network descriptions.  Aborts if any parameter of any class or
    group does not match (within roundoff tolerance). */

#endif

