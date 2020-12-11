#ifndef nmsim_group_synapse_H
#define nmsim_group_synapse_H
 
/* Attributes of synapse groups (synaptic bundles) in a Galves-Löcherbach network. */
/* Last edited on 2020-12-11 14:39:16 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <vec.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

typedef struct nmsim_group_synapse_t
  { nmsim_class_synapse_ix_t isc;      /* Class of synapses in synapse group. */
    nmsim_group_neuron_ix_t ing_pre;   /* Index of pre-synaptic neuron group. */
    nmsim_group_neuron_ix_t ing_pos;   /* Index of post-synaptic neuron group. */
    nmsim_elem_synapse_count_t nse;    /* Number of synapses in group. */
  } nmsim_group_synapse_t;
  /* Attributes of a group of synapses (synaptic bundle) between
    two neuron populations. See {nmsim_group_synapse_INFO} for details. */
  
#define nmsim_group_synapse_INFO \
  "The field {isc}  specifies the class of the synapses in this" \
  " bundle.  It is an index into the table" \
  " of synapse casses ({nmsim_class_synapse_t} records) of the network.\n" \
  "\n" \
  "  The fields {ing_pre} and {ing_pos} are the indices of pre- and" \
  " post-synaptic neuron populations {A} and {B}, respectively.  Every" \
  " synapse in this group connects a neuron of group {A} to a neuron of group {B}.\n" \
  "\n" \
  "  The field {nse} is the number of synapses in the group (may be 0).  It follows" \
  " that if the neuron groups {A} and {B} have {NA} and {NB} neurons," \
  " respectively, the average number of synapses of this group that" \
  " exit a neuron of {A} is {nse/NA}, the averag number that" \
  " enter a neuron of {B} is {nse/NB}, and the average number" \
  " between two such neurons is {nse/(NA*NB)}.\n" \
  "\n" \
  "  Note that as the latter ratio approaches or exceeds 1 there is significant" \
  " probability of multiple synapses between two given neurons.."

/* INPUT/OUTPUT */

void nmsim_group_synapse_write(FILE *wr, nmsim_group_synapse_ix_t isg, nmsim_group_synapse_t *sgrp);
  /* Writes the synapse group attributes {*sgrp} with assumed index {isg}
    to file {wr}, as described by {nmsim_group_synapse_read_INFO} below. */

nmsim_group_synapse_t nmsim_group_synapse_read
  ( FILE *rd, 
    nmsim_group_synapse_ix_t isg, 
    nmsim_class_synapse_ix_t isc_max, 
    nmsim_group_neuron_ix_t ing_max, 
    nmsim_elem_synapse_count_t nse_max
  );
  /* Reads the index and attributes of a synapse group from file {rd}, in the format described 
    by {nmsim_group_synapse_read_INFO} below.  Checks that the index read is indeed {isg},
    that synapse class index {isc} is in {0..isc_max}, that the neuron population indices
    {ing_pre} and {ing_pos} are in {0..ing_max}, and that the number of synapses
    in the group is in {0..nse_max}.  */
    
#define nmsim_group_synapse_read_INFO \
  "A synapse group description is five numbers\n" \
  "     \"{isg} {isc} {ing_pre} {ing_pos} {nse}\"\n" \
  " The numbers must be all on the same line of the" \
  " file, separated by one or more blanks.\n" \
  "\n" \
  "  " nmsim_group_synapse_INFO

/* DEBUGGING AND TESTING */

void nmsim_group_synapse_show(FILE *wr, char *pref, nmsim_group_synapse_t *sgrp, char *suff);
  /* Writes the attributes of {*sgrp} as a compact line, with prefix {pref} 
    and suffix {suff}. */
  
nmsim_group_synapse_t nmsim_group_synapse_throw
  ( nmsim_class_synapse_ix_t isc_max, 
    nmsim_group_neuron_ix_t ing_pre, 
    nmsim_group_neuron_ix_t ing_pos,
    nmsim_elem_synapse_count_t nse_min,
    nmsim_elem_synapse_count_t nse_max
  );
  /* Generates a random synapse group description.  The synapse
    class index {isc} will be a random integer in {0 .. isc_max},
    the pre- and post-synaptic neuron group indices will be 
    {ing_pre} and {ing_pos}, and the number of synapses {nse} 
    will be a random integer in {nse_min .. nse_max}. */

void nmsim_group_synapse_compare(nmsim_group_synapse_t *sgrp_read, nmsim_group_synapse_t *sgrp_orig);
  /* Compares the synapse group attributes {*sgrp_read} read from a file
    with the expected group attributes {*sgrp_orig}.  Aborts if any attribute
    does not match. */

#endif
   
