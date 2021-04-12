#ifndef nmsim_group_neuron_H
#define nmsim_group_neuron_H
 
/* Attributes of a population of Galves-Löcherbach neurons. */
/* Last edited on 2021-01-06 11:09:02 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <vec.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

typedef struct nmsim_group_neuron_t
  { nmsim_class_neuron_ix_t inc;             /* Class of neurons in population. */
    nmsim_elem_neuron_count_t nne;           /* Number of neurons in group. */
    nmsim_elem_neuron_ix_t ine_start;        /* Index of first neuron in group, or 0. */
    nmsim_group_synapse_count_t nsg_out;     /* Number of synaptic bundles out of this group. */
  } nmsim_group_neuron_t;
  /* A record {ngrp} of this type contains the main attributes
    of a neuron population in a GL network.  See {nmsim_group_neuron_INFO} for details. 
  
    The field {.ine_start}(that is not present in network description
    files) is used only when this record is part of an element-level
    network description.  Otherwise it is zero.
    
    The field {.nsg_out} is the number of synaptic bundles out of this neuron group. */
  
#define nmsim_group_neuron_INFO \
  "The attributes of a neuron population in memory are the index {inc} of the" \
  " neuron class to which it belongs, the number of neurons {nne} in" \
  " the group (must be at least 1), and the number {.nsg_out} of" \
  " synaptic bundles that go from that group to other groups."

/* INPUT/OUTPUT */

void nmsim_group_neuron_write(FILE *wr, nmsim_group_neuron_ix_t ing, nmsim_group_neuron_t *ngrp);
  /* Writes the attributes of neuron population {*ngrp} with assumed index {ing}
    to file {wr}, as described by {nmsim_group_neuron_read_INFO} below.  
    The fields {.ine_start} are not written. */

nmsim_group_neuron_t nmsim_group_neuron_read
  ( FILE *rd, 
    nmsim_group_neuron_ix_t ing, 
    nmsim_class_neuron_ix_t inc_max,
    nmsim_elem_neuron_count_t nne_g_max,
    nmsim_group_synapse_count_t nsg_out_max
  );
  /* Reads the index and attributes of a neuron population from file {rd}, in the format described 
    by {nmsim_group_neuron_read_INFO} below.  Checks that the index read from the file is indeed {ing},
    that the neuron class index {inc} is in {0..inc_max}, that the number of neurons
    in the group is in {1..nne_g_max}, and that the number of output synapic bundles {.nsg_out}
    is in 0..nsg_out_max}. The field {.ine_start}, not in the file, is set to zero */
    
#define nmsim_group_neuron_read_INFO \
  "The neuron group description consists of four integers -- the neuron group" \
  " index {ing}, the neuron class index {inc}, the" \
  " number {nne} of neurons in the group, and the number {nsg_out} of" \
  " synaptic bundles out of this neuron group -- separated by one" \
  " or more blanks, on the same line."

/* DEBUGGING AND TESTING */

void nmsim_group_neuron_show(FILE *wr, char *pref, nmsim_group_neuron_t *ngrp, char *suff);
  /* Writes the attributes {ngrp} of a neuron population as a compact line, with prefix {pref} 
    and suffix {suff}. */
 
nmsim_group_neuron_t nmsim_group_neuron_throw
  ( nmsim_class_neuron_ix_t inc_max,
    nmsim_elem_neuron_count_t nne_min, 
    nmsim_elem_neuron_count_t nne_max,
    nmsim_group_synapse_count_t nsg_out_min, 
    nmsim_group_synapse_count_t nsg_out_max
  );
  /* Generates a random neuron group description {ngrp}.  

    The neuron class index {ngrp.inc} will be a random integer in {0 .. inc_max},
    and the number of neurons {ngrp.nne} in the group will be a random integer
    in {nne_min .. nne_max}. Requires {1 <= nne_min <= nne_max} so that there
    will always be at least one neuron in the group.  
    
    The number of outut symaptic bundles {ngrp.nsg_out} from this neuron
    group will be a random number between {nsg_out_min .. nsg_out_max}
    (but not exceeding the number of neurons {ngrp.nne}, so that there
    can be at least one synapse per bundle).  Requires {0 <= nsg_out_min
    <= nsg_out_max}. Note that the neuron group may (will) have no
    output synapses If {nsg_out_min} ({nsg_out_max}) is zero.
    
    The fields {.ine_start} are set to zero.  */

void nmsim_group_neuron_compare(nmsim_group_neuron_t *ngrp_read, nmsim_group_neuron_t *ngrp_orig);
  /* Compares the neuron group attributes {*ngrp_read} read from a file
    with the expected group attributes {*ngrp_orig}.  Ignores the fields
    {.ine_start}. Aborts if any attribute does not match (within roundoff tolerance). */

vec_typedef(nmsim_group_neuron_vec_t,nmsim_group_neuron_vec,nmsim_group_neuron_t);

#endif
   
