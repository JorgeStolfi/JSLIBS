#ifndef nmsim_group_synapse_stats_H
#define nmsim_group_synapse_stats_H
 
/* Representing and printing summary data of a synapse group. */
/* Last edited on 2020-12-15 10:01:32 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <sign.h>

#include <nmsim_basic.h>

typedef struct nmsim_group_synapse_stats_t {
    nmsim_group_synapse_ix_t isg;        /* Index of the synapse group. */
    nmsim_elem_synapse_count_t nse;      /* Total count of synapses in group. */
    nmsim_group_neuron_ix_t ing_pre;     /* Index of pre-synaptic neuron group. */
    nmsim_elem_neuron_count_t nne_pre;   /* Number of neurons in pre-synaptic group. */
    nmsim_group_neuron_ix_t ing_pos;     /* Index of post-synaptic neuron group. */
    nmsim_elem_neuron_count_t nne_pos;   /* Number of neurons in post-synaptic group. */
    /* Statistics about the synapses considered in statistics: */
    nmsim_elem_synapse_count_t nse_cons; /* Count of synapses onsidered. */
    double Wavg; /* Average of synapse weights in group. */
    double Wdev; /* Deviation of synapse weights in group. */
    /* Statistics about total inputs of post-synaptic neurons: */
    double navg; /* Average number of input synapses among post-synaptic neurons. */
    double Wina; /* Average total input synaptic weight among post-synaptic neurons. */
    double Wind; /* Deviation of total input synaptic weight among post-synaptic neurons. */
  } nmsim_group_synapse_stats_t;
  /* A record {st} ofthis type stores summary data for a synapse group,
    or for a subset thereof (such as only the excitatory synapses or only
    the inhibitory ones). 
    
    The {nse} field is the total count of synapses in the set. Synapses
    with zero weight are counted in {nse_zero} and excluded from the averages 
    and deviations. 
    
    The {Wina} and {Wdev} fields are obtained by computing the sum of
    weights of all synapses into each of the {nne_pos} post-synaptic
    neurons, and then taking the average and deviation of those total
    input weights.
    
    The record may be used to store data on an arbitrary set of synapses,
    which may not be contained in a single synaptic group and need not connect two
    specific neuron groups.  In that case the fields {isg,ing_pre,ing_pos} should be set
    to {-1}, and the counts {nne_pre,nne_pos} must include all pre- and post-synaptic
    neurons of the synapses in the set. */
  
void nmsim_group_synapse_stats_print_gen(FILE *wr, nmsim_group_synapse_stats_t *st, int32_t hdr);
void nmsim_group_synapse_stats_print_uns(FILE *wr, nmsim_group_synapse_stats_t *st, int32_t hdr);
void nmsim_group_synapse_stats_print_sgn(FILE *wr, nmsim_group_synapse_stats_t *st, sign_t sgn, int32_t hdr);
  /* If {hdr} is zero, these procedures write to {wr} some of the fields in the record {st}. If
    {hdr} is 1 to 4, they ignore {st} (which may be {NULL}) and instead write
    line {hdr} of the header for those same fields.
    
    {nmsim_group_synapse_stats_print_gen} prints the generic fields {isg}, {nse},
    {ing_pre}, {ing_pos}, {nne_pre} and {nne_pos}.
    
    {nmsim_group_synapse_stats_print_uns} prints only the fields
    {nse_cons}, {Wavg} and {Wina}. It is suitable when {st} refers to 
    a mix of excitatory and inhibitory synapses.
    
    {nmsim_group_synapse_stats_print_sgn} prints the fields {nse_cons},
    {navg}, {Wavg}, {Wdev}, {Wina}, and {Wind}; where {navg} is
    {nse_cons/nne_pos}. It is suitable when {st} refers to a set of
    synapses that are all excitatory or all inhibitory. */

#define nmsim_group_synapse_stats_print_INFO \
  "The table fields in the columns at left have the following meanings:\n" \
  "\n" \
  "  \"syn/group\" : index of the synaptic group.\n" \
  "\n" \
  "  \"total synapses/total\" : total count of synapses in group.\n" \
  "\n" \
  "  \"neur groups/pre\" : index of pre-synaptic neuron group.\n" \
  "\n" \
  "  \"neur groups/post\" : index of post-synaptic neuron group.\n" \
  "\n" \
  "  \"num neurons/pre\" : count of neurons in the pre-synaptic group.\n" \
  "\n" \
  "  \"num neurons/post\" : count of neurons in the post-synaptic group.\n" \
  "\n" \
  "The following fields under the top headings \"nonzero synapses\",\n" \
  "\"excitatory synapses\", and \"inhibitory synapses\" consider only\n" \
  "the subset of synapses whose weight is non-zero, strictly positive,\n" \
  "and strictly negative, respectively:\n" \
  "\n" \
  "  \"count\" : count of synapses in that subset.\n" \
  "\n" \
  "  \"Wavg\" : Average of the weights of the synapses in that subset.\n" \
  "\n" \
  "  \"Wdev\" : Deviation of the weights of the synapses in that subset.\n" \
  "\n" \
  "  \"navg\" : Average, over the neurons in post-synaptig group, of\n" \
  "    the count of synapses from that subset that enter each neuron.\n" \
  "\n" \
  "  \"Wina\" : Average, over the neurons in post-synaptic group, of\n" \
  "    the total weight of input synapses from that subset that enter each neuron.\n" \
  "\n" \
  "  \"Wind\" : Deviation of those weight totals."

#endif

