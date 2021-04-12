#ifndef nmsim_elem_net_sim_group_stats_H
#define nmsim_elem_net_sim_group_stats_H
 
/* Per-group statistical summaries of neuron states through the simulation of a GL net. */
/* Last edited on 2020-12-17 01:18:38 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>
#include <nmsim_stats.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_net.h>
#include <nmsim_elem_net_sim_stats.h>

typedef struct nmsim_elem_net_sim_group_stats_t
  { nmsim_group_neuron_count_t nng;       /* Number of neuron groups. */
    nmsim_time_t tLo;                    /* Discrete time of first state considered. */
    nmsim_time_t tHi;                    /* Discrete time of last state considered. */
    nmsim_elem_net_sim_stats_t **S;       /* Ststistics records for each group. */
  } nmsim_elem_net_sim_group_stats_t;
  /* A data structure that records per-group statistics of the states of neurons
    for the discrete times {tLo..tHi} inclusive, in an elem-level simulation
    of a GL network.
    
    The statistics for each neuron group {ing} are gathered in {gstats.S[ing]},
    which specifies the range of neuron indices in that group.
    The time interval {tLo..tHi} is the union of the time intervals in all
    those records. */
  
nmsim_elem_net_sim_group_stats_t *nmsim_elem_net_sim_group_stats_new
  ( nmsim_group_net_t *gnet,  /* Data on neuron groups, including neuron index ranges. */
    nmsim_time_t tLo,         /* Discrete time of first state to consider. */
    nmsim_time_t tHi          /* Discrete time of last state to consider. */
  );
  /* Allocates an {nmsim_elem_net_sim_group_stats_t} record for the neuron groups
    specified in {gnet} and sets the time ranges to {tLo..tHi}.  Does not initialize
    the individual statistics fields; client must call
    {nmsim_elem_net_sim_group_stats_initialize} for that. 
    
    If {tLo > tHi}, the result is {NULL}. */

void nmsim_elem_net_sim_group_stats_write(char *outPrefix, nmsim_elem_net_sim_group_stats_t *gstats);
  /* Writes the statistical summary {S} for each group in {gstats} in a separate file.
    See {nmsim_elem_net_sim_group_stats_write_INFO}. for file names. If {gstats} is {NULL},
    no file is written. */
    
#define nmsim_elem_net_sim_group_stats_write_INFO \
  "The data for each group is written to a file called \"{outPrefix}_ng{ING}_stats.txt\" where {ING} is" \
  " the neuron group formatted with 10 digits, zero-padded."
    
/* INCREMENTAL STATS GATHERING

  An {nmsim_elem_net_sim_group_stats_t} record set {gstats} can be computed incrementally by 
  
    calling {nmsim_elem_net_sim_group_stats_clear} just once, then
    
    calling {nmsim_elem_net_sim_group_stats_accumulate_V_age_M_H} once for each time {t}, then
    
    calling {nmsim_elem_net_sim_group_stats_accumulate_VF_AF_X_I_J} once for each time step
    from {t} to {t+1} (after computing {X,I,J}, but before updating the {V,age}),
    and then
    
    calling {nmsim_elem_net_sim_group_stats_finalize} just once.
    
  Between clearing and finalizing, the subfields {.avg} and {.dev} of all
  {gstats} fields are temporarily used to hold the sum of samples and the sum of their
  squares, respetively. */ 

void nmsim_elem_net_sim_group_stats_initialize(nmsim_elem_net_sim_group_stats_t *gstats);
  /* Clears all the fields of {S} in preparation for 
    {nmsim_elem_net_sim_group_stats_accumulate}. A no-op if {gstats} is {NULL}. */
   
void nmsim_elem_net_sim_group_stats_accumulate_V_age_M_H
  ( nmsim_elem_net_sim_group_stats_t *S,
    nmsim_time_t t,
    nmsim_elem_neuron_count_t nne,
    double V[],                /* Neuron potentials at time {t}. */
    nmsim_step_count_t age[],  /* Ages of neurons at time {t}. */
    double M[],                /* Recharge modulators at time {t}. */
    double H[]                 /* Output modulators at time {t}. */
  );
  /* Accumulates statistics of neuron state variables {V,age,M,H}
    for one time value {t} and all neuron groups in {gstats}.  See 
    {nmsim_elem_net_sim_stats_accumulate_VF_AF_X_I_J}.
    
    The procedure does nothing if{gstats} is {NULL} or {t} is not in {S.tLo..S.tHi}. */
    
void nmsim_elem_net_sim_group_stats_accumulate_VF_AF_X_I_J
  ( nmsim_elem_net_sim_group_stats_t *S,
    nmsim_time_t t,
    nmsim_elem_neuron_count_t nne,
    double V[],                /* Neuron potentials at time {t}. */
    nmsim_step_count_t age[],  /* Ages of neurons at time {t}. */
    bool_t X[],                /* Firing indicators for step {t} to {t+1}. */
    double I[],                /* External neuron inputs for step {t} to {t+1}. */
    double J[]                 /* Total inputs for step {t} to {t+1}. */
  );
  /* Accumulates statistics of neuron state variables variables that refer to one time step from
    {t} to {t+1}, namely {X,I,J}, for all neuron groups in {gstats}.
    Also accumulates statistcs on {V,age} only for neurons that fired in that time step.
    
    This procedure should be called after computing {X,I,J} but before 
    updating {V,age}, so that they still refer to time {t}.
    
    The procedure does nothing if{gstats} is {NULL} or {t} is not in {S.tLo..S.tHi}. */

void nmsim_elem_net_sim_group_stats_finalize(nmsim_elem_net_sim_group_stats_t *S);
  /* Converts all {.avg} subfields in {gstats} from 
    sum of sample values to their average, and the {.dev} subfields
    from sum of squared samples to their standard  deviation.
    See {nmsim_stats_finalize} for more details. */

#endif
