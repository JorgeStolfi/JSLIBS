#ifndef nmsim_elem_net_sim_stats_H
#define nmsim_elem_net_sim_stats_H
 
/* Statistical summary of a range of neurons in a Galves-Löecherbach net. */
/* Last edited on 2020-12-25 11:49:06 by jstolfi */

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

typedef struct nmsim_elem_net_sim_stats_t
  { nmsim_elem_neuron_ix_t ineLo;           /* Index of the first neuron of the set. */
    nmsim_elem_neuron_ix_t ineHi;           /* Index of the last neuron of the set. */
    nmsim_time_t tLo;                     /* Discrete time of first state considered. */
    nmsim_time_t tHi;                     /* Discrete time of last state considered. */
    
    nmsim_stats_t V;                      /* Stats of potential. */
    nmsim_stats_t age;                    /* Stats of neuron age (steps since last firing). */
    nmsim_stats_t M;                      /* Stats of recharge factor modulator. */
    nmsim_stats_t H;                      /* Stats of output synapse strength modulator. */
    nmsim_stats_t X;                      /* Stats of firing indicator. */
    nmsim_stats_t I;                      /* Stats of external input. */
    nmsim_stats_t S;                      /* Stats of synaptic input. */
    nmsim_stats_t J;                      /* Stats of total input. */

    nmsim_stats_t V_fire;                 /* Stats of potential just before firing. */
    nmsim_stats_t age_fire;               /* Stats of age just before firing. */
  } nmsim_elem_net_sim_stats_t;
  /* A data structure that records statistics of the states of neurons with indices
    {ineLo..ineHi} over the discrete times {tLo..tHi} inclusive, in an 
    elem-level simulation of a GL network.
    
    The statistics for {V,age,M,H} reflect their values at the times {tLo..tHi}.
    The statistics for {X,I,J,V_fire,age_fire} reflects only the values of those
    variables for time {tLo..tHi-1}, since they are defined for time steps, not 
    for single times. */
  
nmsim_elem_net_sim_stats_t *nmsim_elem_net_sim_stats_new
  ( nmsim_elem_neuron_ix_t ineLo,    /* Index of the first neuron of the set. */
    nmsim_elem_neuron_ix_t ineHi,    /* Index of the last neuron of the set. */
    nmsim_time_t tLo,                /* Discrete time of first state to consider. */
    nmsim_time_t tHi                 /* Discrete time of last state to consider. */
  );
  /* Allocates an {nmsim_elem_net_sim_stats_t} record and sets its neuron and time 
    ranges as specified.  Does not initialize the fields; client must call
    {nmsim_elem_net_sim_stats_initialize} for that.  */

#define nmsim_elem_net_sim_stats_FILE_TYPE "nmsim_elem_net_sim_stats"
    
#define nmsim_elem_net_sim_stats_VERSION "2020-12-16"

void nmsim_elem_net_sim_stats_write(FILE *wr, nmsim_elem_net_sim_stats_t *S);
  /* Writes the statistical summary {S} of a set of neurons to file {wr}. */
    
/* INCREMENTAL STATS GATHERING

  An {nmsim_elem_net_sim_stats_t} record {S} can be computed incrementally by 
  
    calling {nmsim_elem_net_sim_stats_clear} just once, then
    
    calling {nmsim_elem_net_sim_stats_accumulate_V_age_M_H} once for each time {t}, then
    
    calling {nmsim_elem_net_sim_stats_accumulate_VF_AF_X_I_J} once for each time step
    from {t} to {t+1} (after computing {X,I,J}, but before updating the {V,age}), 
    and then
    
    calling {nmsim_elem_net_sim_stats_finalize} just once.
    
  Between clearing and finalizing, the subfields {.avg} and {.dev} of all
  {S} fields are temporarily used to hold the sum of samples and the sum of their
  squares, respetively. */ 

void nmsim_elem_net_sim_stats_initialize(nmsim_elem_net_sim_stats_t *S);
  /* Clears all the fields of {S} in preparation for {nmsim_elem_net_sim_stats_accumulate}. */
   
void nmsim_elem_net_sim_stats_accumulate_V_age_M_H
  ( nmsim_elem_net_sim_stats_t *S,
    nmsim_time_t t,
    nmsim_elem_neuron_count_t nne,
    double V[],                /* Neuron potentials at time {t}. */
    nmsim_step_count_t age[],  /* Ages of neurons at time {t}. */
    double M[],                /* Recharge modulators at time {t}. */
    double H[]                 /* Output modulators at time {t}. */
  );
  /* Accumulates statistics of neuron state variables {V,age,M,H}
    for one time value {t}.
    
    Specifically, accumulates {V[ine],age[ine],M[ine],H[ine]}
    into {S.V,S.age,S.M,S.H}, respectively.
    
    Considers only neurons in the range {S.ineLo..S.ineHi} and only if
    {t} is in {S.tLo..S.tHi}. */
    
void nmsim_elem_net_sim_stats_accumulate_VF_AF_X_I_J
  ( nmsim_elem_net_sim_stats_t *S,
    nmsim_time_t t,
    nmsim_elem_neuron_count_t nne,
    double V[],                /* Neuron potentials at time {t}. */
    nmsim_step_count_t age[],  /* Ages of neurons at time {t}. */
    bool_t X[],                /* Firing indicators for step {t} to {t+1}. */
    double I[],                /* External neuron inputs for step {t} to {t+1}. */
    double J[]                 /* Total inputs for step {t} to {t+1}. */
  );
  /* Accumulates statistics variables that refer to one time step from
    {t} to {t+1}, namely {X,I,J}. Also accumulates statistcs on {V,age} only 
    for neurons that fired in that time step. 
    
    This procedure should be called after computing {X,I,J} but before 
    updating {V,age}, so that they still refer to time {t}.
    
    Considers only neurons in the range {S.ineLo..S.ineHi} and only if
    {t} is in {S.tLo..S.tHi-1}. */

void nmsim_elem_net_sim_stats_finalize(nmsim_elem_net_sim_stats_t *S);
  /* For each field of of {S}, converts the {.avg} subfield from 
    sum of sample values to their average, and the {.dev} subfield
    from sum of squared samples to their standard  deviation.
    
    See {nmsim_stats_finalize} for more details. */

#endif
