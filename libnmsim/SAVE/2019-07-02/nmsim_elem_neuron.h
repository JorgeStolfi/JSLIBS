#ifndef nmsim_elem_neuron_H
#define nmsim_elem_neuron_H
 
/* Neurons in an element-level network. */
/* Last edited on 2019-03-28 18:17:07 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <vec.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

typedef struct nmsim_elem_neuron_t 
  { 
    nmsim_group_neuron_ix_t ing; /* Index of neuron group. */
  } nmsim_elem_neuron_t;
  /* Specifies a neuron.  See {nmsim_elem_neuron_INFO}. */
  
#define nmsim_elem_neuron_INFO \
  "The only attribute of a neuron is its group index {ing}."

/* INPUT/OUTPUT */

void nmsim_elem_neuron_write(FILE *wr, nmsim_elem_neuron_ix_t ine, nmsim_elem_neuron_t *neu);
  /* Writes the attributes of neuron {*neu} with assumed index {ine}
    to file {wr}, as described by {nmsim_elem_neuron_read_INFO} below. */

nmsim_elem_neuron_t nmsim_elem_neuron_read
  ( FILE *rd, 
    nmsim_elem_neuron_ix_t ine, 
    nmsim_group_neuron_ix_t ing_max
  );
  /* Reads the index and attributes of a neuron from file {rd}, in the format described 
    by {nmsim_elem_neuron_read_INFO} below.  Checks that the index read is indeed {ine} and
    that the neuron group index {inc} is in {0..ing_max}.  */
    
#define nmsim_elem_neuron_read_INFO \
  "The neuron description is two numbers \"{ine} {ing}\" on the same line of the" \
  " file, separated by one or more blanks.\n" \
  "\n" \
  "  " nmsim_elem_neuron_INFO

/* DEBUGGING AND TESTING */

void nmsim_elem_neuron_show(FILE *wr, char *pref, nmsim_elem_neuron_t *neu, char *suff);
  /* Writes the neuron attributes {*neu} as a compact line, with prefix {pref} 
    and suffix {suff}. */
  
nmsim_elem_neuron_t nmsim_elem_neuron_throw(nmsim_group_neuron_ix_t ing_max);
  /* Generates a random neuron description.  The neuron group index {ing}
    will be random in {0 .. ing_max}. */

void nmsim_elem_neuron_compare(nmsim_elem_neuron_t *neu_read, nmsim_elem_neuron_t *neu_orig);
  /* Compares the neuron attributes {*neu_read} read from a file
    with the expected attributes {*neu_orig}.  Aborts if any attribute
    does not match (within roundoff tolerance). */

vec_typedef(nmsim_elem_neuron_vec_t, nmsim_elem_neuron_vec, nmsim_elem_neuron_t);
  /* Type of an extensible vector of neurons. */

#endif
