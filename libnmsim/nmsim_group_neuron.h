#ifndef nmsim_group_neuron_H
#define nmsim_group_neuron_H
 
/* Attributes of a population of Galves-Löcherbach neurons. */
/* Last edited on 2019-03-28 18:22:10 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <vec.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

typedef struct nmsim_group_neuron_t
  { nmsim_class_neuron_ix_t inc;    /* Class of neurons in population. */
    nmsim_elem_neuron_count_t nne;  /* Number of neurons in group. */
  } nmsim_group_neuron_t;
  /* Attributes of a GL neuron population. See {nmsim_group_neuron_INFO} for details. */
  
#define nmsim_group_neuron_INFO \
  "The field {inc} specifies the parameters of the neurons in this" \
  " population.  It an index into the table of neuron classes" \
  " ({nmsim_class_neuron_t} records) of the network.  The field {nne} is" \
  " the number of neurons in the group (at least 1)."

/* INPUT/OUTPUT */

void nmsim_group_neuron_write(FILE *wr, nmsim_group_neuron_ix_t ing, nmsim_group_neuron_t *ngrp);
  /* Writes the attributes of neuron population {*ngrp} with assumed index {ing}
    to file {wr}, as described by {nmsim_group_neuron_read_INFO} below. */

nmsim_group_neuron_t nmsim_group_neuron_read
  ( FILE *rd, 
    nmsim_group_neuron_ix_t ing, 
    nmsim_class_neuron_ix_t inc_max
  );
  /* Reads the index and attributes of a neuron population from file {rd}, in the format described 
    by {nmsim_group_neuron_read_INFO} below.  Checks that the index read from the file is indeed {ing}
    and that the neuron class index {inc} is in {0..inc_max}. */
    
#define nmsim_group_neuron_read_INFO \
  "The neuron group description consists of three integers -- the neuron group" \
  " index {ing}, the neuron class index {inc}, and" \
  " the number {nne} of neurons in the group -- separated by one" \
  " or more blanks, on the same line."

/* DEBUGGING AND TESTING */

void nmsim_group_neuron_show(FILE *wr, char *pref, nmsim_group_neuron_t *ngrp, char *suff);
  /* Writes the attributes {ngrp} of a neuron population as a compact line, with prefix {pref} 
    and suffix {suff}. */
 
nmsim_group_neuron_t nmsim_group_neuron_throw
  ( nmsim_class_neuron_ix_t inc_max,
    nmsim_elem_neuron_count_t nne_g_min,
    nmsim_elem_neuron_count_t nne_g_max
  );
  /* Generates a random neuron group description.  The 
    neuron class index {inc} will be a random integer in {0 .. inc_max},
    and the number of neurons {nne} in the group will be a random integer
    in {nne_g_min .. nne_g_max}. */

void nmsim_group_neuron_compare(nmsim_group_neuron_t *ngrp_read, nmsim_group_neuron_t *ngrp_orig);
  /* Compares the neuron group attributes {*ngrp_read} read from a file
    with the expected group attributes {*ngrp_orig}.  Aborts if any attribute
    does not match (within roundoff tolerance). */

vec_typedef(nmsim_group_neuron_vec_t,nmsim_group_neuron_vec,nmsim_group_neuron_t);

#endif
   
