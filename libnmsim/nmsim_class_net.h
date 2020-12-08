#ifndef nmsim_class_net_H
#define nmsim_class_net_H
 
/* Class-level description of a network of Galves-Löcherbach neurons. */
/* Last edited on 2019-06-18 11:17:49 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>

typedef struct nmsim_class_net_t
  { /* Neuron classes: */
    nmsim_class_neuron_count_t nnc;        /* Count of distinct neuron classes. */
    nmsim_class_neuron_ref_vec_t nclass;   /* Parameters for each neuron class. */
    /* Synapse classes: */
    nmsim_class_synapse_count_t nsc;       /* Count of distinct synapse classes. */
    nmsim_class_synapse_ref_vec_t sclass;  /* Parameters for each synapse class. */
  } nmsim_class_net_t;
  /* Class-level description of a network of Galves-Loecherbach neurons,
    
    The neurons are divided into {nnc} distinct classes. Each neuron
    class has an index {inc} in {0.. nnc-1}, and its parameters are
    {nclass.e[inc]}.
    
    The synapses are likewise divided into {nsc} distinct classes.
    Each synapse class has an index {isc} in {0..nsc-1}, and its
    parameters are {sclass.e[isc]}.
    
    The allocated size of the {nclass} vector is {nclass.ne} which is
    usually larger than then number of used entries {nnc}. The same
    applies to the vector {sclass}, which has only {nsc} used entries. */
    
/* NETWORK CREATION */

nmsim_class_net_t *nmsim_class_net_new(void);
  /* Allocates an {nmsim_class_net_t} structure, initially with zero neuron
    and synapse classes. */

void nmsim_class_net_free(nmsim_class_net_t *cnet);
  /* Releases all storage of the class-level network description {cnet},
    including internal tables and parameter records, and the top-level
    descriptor itself. */

nmsim_class_neuron_ix_t nmsim_class_net_add_neuron_class
  ( nmsim_class_net_t *cnet, 
    nmsim_class_neuron_t *nclass
  );
  /* Adds to the network {cnet} a new neuron class, described by the
    given neuron parameter record {*nclass}. Assigns to those
    parameters the next neuron class index {inc}and stores {parm} in
    the {cnet->nclass} vector, expanding it if needed, and incrementng
    {cnet->nnc}. Returns the class index {inc}. */

nmsim_class_synapse_ix_t nmsim_class_net_add_synapse_class
  ( nmsim_class_net_t *cnet, 
    nmsim_class_synapse_t *sclass
  );
  /* Adds to the network {cnet} a new synapse class, described by the
    given synapse parameter record {*sclass}. Assigns to those
    parameters the next synapse class index {isc} and stores {parm} in
    the {cnet->sclass} vector, expanding it if needed, and incrementng
    {cnet->nsc}. Returns the index {isc}. */

/* INPUT/OUTPUT */

void nmsim_class_net_write(FILE *wr, nmsim_class_net_t *cnet, double timeStep);
  /* Writes the class-level description {cnet} of a network to file {wr},
    in the format described by {nmsim_class_net_read_INFO} below.  
    
    The {timeStep} parameter (in milliseconds) is used to convert
    the decay factors to step-invariant characteristic decay
    times. */

nmsim_class_net_t *nmsim_class_net_read(FILE *rd, double timeStep);
  /* Reads the class-level description of a networkfrom file {rd}, in
    the format described by {nmsim_class_net_read_INFO} below. 
    
    The {timeStep} parameter (in milliseconds) is used to convert
    step-invariant characteristic decay times of the neurons into the
    step-dependent decay factors. */
    
#define nmsim_class_net_read_INFO \
  "The description of the class-level network begins with a line \n" \
  "    begin " nmsim_class_net_FILE_TYPE " (format of " nmsim_class_net_VERSION ")\n" \
  " followed by two blocks of lines that describe the neuron and synapse classes, and" \
  " a closing line \"end " nmsim_class_net_FILE_TYPE "\".\n" \
  "\n" \
  "  The first block is preceded by a line \"neuron_classes = {nnc}\", where" \
  " {nnc} is the number of distinct neuron classes.  Then there" \
  " follows {nnc} neuron class descriptions (see below), each preceded by a" \
  " line \"class = {inc}\" where {inc} the class index in {0 .. nnc}.\n" \
  "\n" \
  "  The second block of lines is similar, except that it is preceded" \
  " by a line \"synapse_classes = {nsc}\" with the number {nsc} of" \
  " synapse classes; and each entry (see below) is a synapse class description.\n" \
  "\n" \
  "NEURON CLASS DESCRIPTION\n" \
  "  " nmsim_class_neuron_read_INFO "\n" \
  "\n" \
  "SYNAPSE CLASS DESCRIPTION\n" \
  "  " nmsim_class_synapse_read_INFO ""
    
#define nmsim_class_net_FILE_TYPE "nmsim_class_net"
    
#define nmsim_class_net_VERSION "2019-01-10"

/* TESTING AND DEBUGGING */

nmsim_class_net_t *nmsim_class_net_throw
  ( nmsim_class_neuron_count_t nnc, 
    nmsim_class_synapse_count_t nsc
  );
  /* Creates a random class-level network description with {nnc} neuron classes
    and {nsc} synapse classes.*/

void nmsim_class_net_compare(nmsim_class_net_t *cnet_read, nmsim_class_net_t *cnet_orig);
  /* Compares the class-level network description {cnet_read} read from a file
    with the expected description {cnet_orig}.  Aborts if any 
    parameter of any class does not match (within roundoff tolerance). */

#endif

