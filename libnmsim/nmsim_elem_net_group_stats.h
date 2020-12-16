#ifndef nmsim_elem_net_group_stats_H
#define nmsim_elem_net_group_stats_H
 
/* Summaries of neuron synapse group data in an element-level GL network model. */
/* Last edited on 2020-12-13 21:17:27 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <nmsim_group_synapse_stats.h>
#include <nmsim_elem_net.h>

nmsim_group_synapse_stats_t *nmsim_elem_net_group_stats_get
  ( nmsim_elem_net_t *enet, 
    nmsim_group_synapse_ix_t isg,
    sign_t sgn
  );
  /* Returns a record {st} with startistcs of the synapses in synapse group {isg} of the
    element-level network [enet}.  
    
    If {isg} is zero, considers all synapses; but uses only those with
    nonzero weight in averages {st->Wavg}, {st->Wdev}, and {st->navg},
    counting those with zero weight in st->nse_zero).
    
    If {isg} is {+1} or {-1}, considers only those synapses whose weight
    has the same sign as {isg}.  In these cases the field {st->nse_zero} will
    be zero. */
  
void nmsim_elem_net_group_stats_print_one
  ( FILE *wr, 
    nmsim_elem_net_t *enet,
    nmsim_group_synapse_ix_t isg,
    char *pref,
    int32_t hdr,
    char *suff
  );
  /* If {hdr} is zero, writes to {wr} one line of a table with
    statistics of the synaptic group {isg} in network [enet},
    considering all non-zero synapses, only the exitatory ones, axsand only
    the inhibitory ones.
    
    If {hdr} is 1 to 4, ignores {enet} and {isg}, and prints instead
    header line number {hdr} of the table.
    
    In any case, the line is prefixed by {pref} (if not {NULL}) and
    followed by {suff} (ditto). */
  
void nmsim_elem_net_group_stats_print_all
  ( FILE *wr, 
    nmsim_elem_net_t *enet,
    char *pref,
    char *suff
  );
  /* Writes to {wr} a table with statistics of all synaptic groups in
    network {enet}, including the table header and an explanatory
    legend below.
    
    Each line, including in the header, is prefixed by {pref} (if not
    {NULL}) and followed by {suff} (ditto). The legend is printed flush
    left. */

#endif

