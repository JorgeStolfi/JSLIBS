#ifndef nmsim_class_net_H
#define nmsim_class_net_H
 
/* Class-level description of a network of Galves-Löcherbach neurons. */
/* Last edited on 2020-12-12 08:19:59 by jstolfi */

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
    nmsim_class_neuron_count_t nnc;     /* Count of distinct neuron classes (at least 1). */
    nmsim_class_neuron_t **nclass;      /* Parameters for each neuron class. */
    /* Synapse classes: */
    nmsim_class_synapse_count_t nsc;    /* Count of distinct synapse classes (may be zero). */
    nmsim_class_synapse_t **sclass;     /* Parameters for each synapse class. */
  } nmsim_class_net_t;
  /* Class-level description of a network of Galves-Loecherbach neurons,
    
    The neurons are divided into {nnc} distinct classes. Each neuron
    class has an index {inc} in {0.. nnc-1}, and its parameters are
    {*(nclass[inc])}. There must be at least one neuron class.
    
    The synapses are likewise divided into {nsc} distinct classes.
    Each synapse class has an index {isc} in {0..nsc-1}, and its
    parameters are {*(sclass[isc])}. There may be zero synapse classes. */
    
/* NETWORK CREATION */

nmsim_class_net_t *nmsim_class_net_new
  ( nmsim_class_neuron_count_t nnc, 
    nmsim_class_synapse_count_t nsc
  );
  /* Allocates an {nmsim_class_net_t} structure, with tables for {nnc} neuron classes
    and {nsc} synapse classes. */

void nmsim_class_net_free(nmsim_class_net_t *cnet);
  /* Releases all storage of the class-level network description {cnet},
    including internal tables and parameter records, and the top-level
    descriptor itself. */

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
  " followed by:\n" \
  "      * a line \"neuron_classes = {nnc}\" where {nnc} is the number of neuron classes\n" \
  "      * a line \"synapse_classes = {nsc}\" where {nsc} is the number of synapse classes\n" \
  "      * a set of {nsc} lines that describe the neuron classes (see below)\n" \
  "      * a set of {nsc} lines that describe the synapse classes (ditto)\n" \
  "  The description ends with the line \"end " nmsim_class_net_FILE_TYPE "\".\n" \
  "\n" \
  "DESCRIPTION OF A NEURON CLASS\n" \
  "  Each neuron class is preceded by a line \"neuron_class = {inc}\" where" \
  " {inc} is the index of the class.\n" \
  "\n" \
  "  " nmsim_class_neuron_read_INFO "\n" \
  "\n" \
  "DESCRIPTION OF A SYNAPSE CLASS\n" \
  "  Each synapse class is preceded by a line \"synpse_class = {isc}\" where" \
  " {isc} is the index of the class.\n" \
  "\n" \
  "  " nmsim_class_synapse_read_INFO ""
    
#define nmsim_class_net_FILE_TYPE "nmsim_class_net"
    
#define nmsim_class_net_VERSION "2020-12-10"

/* TESTING AND DEBUGGING */

void nmsim_class_net_compare(nmsim_class_net_t *cnet_read, nmsim_class_net_t *cnet_orig);
  /* Compares the class-level network description {cnet_read} read from a file
    with the expected description {cnet_orig}.  Aborts if any 
    parameter of any class does not match (within roundoff tolerance). */

#endif

