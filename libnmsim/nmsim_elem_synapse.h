#ifndef nmsim_elem_synapse_H
#define nmsim_elem_synapse_H
 
/* Synapses in an element-level network. */
/* Last edited on 2020-12-06 18:57:35 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <vec.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

typedef struct nmsim_elem_synapse_t 
  { nmsim_group_synapse_ix_t isg;      /* Index of synaptic group. */
    nmsim_elem_neuron_ix_t ine_pre;    /* Index of pre-synaptic neuron. */
    nmsim_elem_neuron_ix_t ine_pos;    /* Index of post-synaptic neuron. */
    double W;                          /* Resting weight of synapse (mV). */
  } nmsim_elem_synapse_t;
  /* Specifies a chemical synapse.  See {nmsim_elem_synapse_INFO}. */
  
#define nmsim_elem_synapse_INFO \
  "A synapse is defined by a synapse group index {isg}, the" \
  " indices of the presynaptic neuron {ine_pre} and post-synaptic" \
  " neuron {ine_pos}, and a resting weight (or strength) {W}.\n" \
  "\n" \
  "  When the pre-synaptic neuron fires, after a long interval without firing, the" \
  " potential of the post-synaptic neuron increases by {W} millivots.\n" \
  "\n" \
  "  The strength {W}  can be negative; this option can be used to" \
  " simulate inhibitory synapses.  If {W} is zero, the synapse" \
  " is inactive (as if it did not exist)."

/* INPUT/OUTPUT */

void nmsim_elem_synapse_write(FILE *wr, nmsim_elem_synapse_ix_t ise, nmsim_elem_synapse_t *syn);
  /* Writes the attributes of synapse {*syn} with assumed index {ise}
    to file {wr}, as described by {nmsim_elem_synapse_read_INFO} below. */

nmsim_elem_synapse_t nmsim_elem_synapse_read
  ( FILE *rd, 
    nmsim_elem_synapse_ix_t ise, 
    nmsim_group_synapse_ix_t isg_max, 
    nmsim_elem_neuron_ix_t ine_max
  );
  /* Reads the index and attributes of a synapse from file {rd}, in the format described 
    by {nmsim_elem_synapse_read_INFO} below.  Checks that the index read is indeed {ise},
    that synapse group index {isc} is in {0..isg_max}, and that the neuron indices
    {ine_pre} and {ine_pos} are in {0..ine_max}.  */
    
#define nmsim_elem_synapse_read_INFO \
  "The synapse description is five numbers\n" \
  "     \"{ise} {isg} {ine_pre}{inw_pos} {W}\"\n" \
  " The numbers must be all on the same line of the" \
  " file, separated by one or more blanks.\n" \
  "\n" \
  "  " nmsim_elem_synapse_INFO

/* DEBUGGING AND TESTING */

void nmsim_elem_synapse_show(FILE *wr, char *pref, nmsim_elem_synapse_t *syn, char *suff);
  /* Writes the synapse attributes {*syn} as a compact line, with prefix {pref} 
    and suffix {suff}. */
  
nmsim_elem_synapse_t nmsim_elem_synapse_throw
  ( nmsim_group_synapse_ix_t isg_min, 
    nmsim_group_synapse_ix_t isg_max, 
    nmsim_elem_neuron_ix_t ine_pre_min,
    nmsim_elem_neuron_ix_t ine_pre_max, 
    nmsim_elem_neuron_ix_t ine_pos_min,
    nmsim_elem_neuron_ix_t ine_pos_max,
    double W_avg,
    double W_dev
  );
  /* Generates a random synapse description.  The synapse group {isg}
    will be random in {isg_min .. isg_max}, and the pre-synaptic
    neuron index will be random in {ine_pre_min .. ine_pre_max},
    and the  and post-synaptic one will be in {ine_pos_min .. ine_pos_max}. 
    
    If {W_dev} is zero, the synaptic weight {W} will be equal to {W_avg}.
    Otherwise {W_avg} must not be zero, and {W} will be a variable with log-normal distribution of 
    mean {|W_avg|} and deviation {W_dev}, with the same
    sign as {W_avg}. */

void nmsim_elem_synapse_compare(nmsim_elem_synapse_t *syn_read, nmsim_elem_synapse_t *sin_orig);
  /* Compares the synapse attributes {*syn_read} read from a file
    with the expected attributes {*sin_orig}.  Aborts if any attribute
    does not match (within roundoff tolerance). */

vec_typedef(nmsim_elem_synapse_vec_t, nmsim_elem_synapse_vec, nmsim_elem_synapse_t);
  /* Type of an extensible vector of synapses. */

#endif
