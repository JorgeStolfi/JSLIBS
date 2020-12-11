#ifndef nmsim_elem_neuron_H
#define nmsim_elem_neuron_H
 
/* Neurons in an element-level network. */
/* Last edited on 2020-12-11 01:42:20 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

typedef struct nmsim_elem_neuron_t 
  { 
    nmsim_group_neuron_ix_t ing;            /* Index of neuron group. */
    nmsim_elem_synapse_ix_t nse_out;        /* Number of output synapses (may be zero). */
    nmsim_elem_synapse_ix_t ise_out_start;  /* Index of first output synapse, or 0. */
  } nmsim_elem_neuron_t;
  /* Specifies a neuron.  See {nmsim_elem_neuron_INFO}.
  
    The field {.nse_out} is the number of synapses out of this neuron.
    The field {.ise_out_start} (that is not present in network description
    files) is the index of the first of those synapses. 
    
    More precisely, for a neuron with index {ine}, the field
    {.ise_out_start} will be the be the sum of {.nse_out} for all previous
    neurons {0..ine-1}. This will be true even if neuron {ine} has no
    output synapses. This rule assumes that the synapses are sorted by
    the index of the pre-synaptic neuron. */
  
#define nmsim_elem_neuron_INFO \
  "The attributes of a neuron in memory are the index {ing} of the neuron" \
  " group to which it belongs, and the number {nse_out} of synapses" \
  " that go from it to other neurons."

/* INPUT/OUTPUT */

void nmsim_elem_neuron_write(FILE *wr, nmsim_elem_neuron_ix_t ine, nmsim_elem_neuron_t *neu);
  /* Writes the attributes of neuron {*neu} with assumed index {ine}
    to file {wr}, as described by {nmsim_elem_neuron_read_INFO} below.
    The field {neu.ise_out_start} is not written. */

nmsim_elem_neuron_t nmsim_elem_neuron_read
  ( FILE *rd, 
    nmsim_elem_neuron_ix_t ine, 
    nmsim_group_neuron_ix_t ing_max, 
    nmsim_elem_synapse_count_t nse_out_max
  );
  /* Reads the index and attributes of a neuron from file {rd}, in the
    format described by {nmsim_elem_neuron_read_INFO} below. Checks that
    the index read is indeed {ine}, that the neuron group index {.ing}
    is in {0..ing_max}, that the number of output synapses {.nse_out}
    is in {0..nse_out_max}.  The field {.ise_out_start} is set to 0, and
    should be defined by the caller. */
    
#define nmsim_elem_neuron_read_INFO \
  "The neuron description is three numbers \"{ine} {ing} {nse_out}\" on" \
  " the same line of the file, separated by one or more blanks.\n" \
  "\n" \
  "  " nmsim_elem_neuron_INFO

/* DEBUGGING AND TESTING */

void nmsim_elem_neuron_show(FILE *wr, char *pref, nmsim_elem_neuron_t *neu, char *suff);
  /* Writes the neuron attributes {*neu} as a compact line, with prefix {pref} 
    and suffix {suff}. */
  
nmsim_elem_neuron_t nmsim_elem_neuron_throw(nmsim_group_neuron_ix_t ing_max);
  /* Generates a random neuron description.  The neuron group index {ing}
    will be random in {0 .. ing_max}.  Both the output synapse count
    {nse_out} and first output synapse index {ise_out_start} will be 
    set to zero, and must be fixed by the client. */

void nmsim_elem_neuron_compare
  ( nmsim_elem_neuron_t *neu_read, 
    nmsim_elem_neuron_t *neu_orig
  );
  /* Compares the neuron attributes {*neu_read} read from a file
    with the expected attributes {*neu_orig}.  
    Ignores the auxiliary {ise_out_start} fields.  Aborts if any attribute
    does not match (within roundoff tolerance). */

#endif
