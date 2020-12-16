#ifndef nmsim_elem_net_stats_H
#define nmsim_elem_net_stats_H
 
/* Statistical summary of a range of neurons in a Galves-Löecherbach net. */
/* Last edited on 2020-12-15 21:52:31 by jstolfi */

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

typedef struct nmsim_elem_net_stats_t
  { nmsim_elem_neuron_ix_t ilo;           /* Index of the first neuron of the set. */
    nmsim_elem_neuron_ix_t ihi;           /* Index of the last neuron of the set. */
    nmsim_time_t tlo;                     /* Discrete time of first state considered. */
    nmsim_time_t thi;                     /* Discrete time of last state considered. */
    
    nmsim_stats_t V;                      /* Stats of potential. */
    nmsim_stats_t age;                    /* Stats of neuron age (steps since last firing). */
    nmsim_stats_t M;                      /* Stats of recharge factor modulator. */
    nmsim_stats_t H;                      /* Stats of output synapse strength modulator. */
    nmsim_stats_t X;                      /* Stats of firing indicator. */
    nmsim_stats_t I;                      /* Stats of external input. */
    nmsim_stats_t J;                      /* Stats of total input. */

    nmsim_stats_t V_fire;                 /* Stats of potential just before firing. */
    nmsim_stats_t age_fire;               /* Stats of age just before firing. */
  } nmsim_elem_net_stats_t;
  /* A data structure that records statistics of the neurons with indices
    {ilo..ihi} for the discrete times {tlo..thi} inclusive, in an elem-level simulation
    of a GL network.
    
    The statistics for {V,age,M,H} reflect their values at the times {tlo..thi}.
    The statistics for {X,I,J,V_fire,age_fire} reflects only the values of those
    variables for time {tlo..thi-1}, since they are defined for time steps, not 
    for single times. */
  
#define nmsim_elem_net_stats_FILE_TYPE "nmsim_elem_net_stats"
    
#define nmsim_elem_net_stats_VERSION "2020-12-15"

void nmsim_elem_net_stats_write(FILE *wr, nmsim_elem_net_stats_t *ntS);
  /* Writes the statistical summary {ntS} of a set of neurons to file {wr}. */
    
/* INCREMENTAL STATS GATHERING

  An {nmsim__elem_netstats_t} record {ntS} can be computed incrementally by 
  
    calling {nmsim_elem_net_stats_clear} just once, then
    
    calling {nmsim_elem_net_stats_accumulate_V_age_M_H} once for each time {t}, then
    
    calling {nmsim_elem_net_stats_accumulate_VF_AF_X_I_J} once for each time step
    from {t} to {t+1}, and then
    
    calling {nmsim_elem_net_stats_finalize} just once.
    
  Between clearing and finalizing, the subfields {.avg} and {.dev} of all
  {ntS} fields are temporarily used to hold the sum of samples and the sum of their
  squares, respetively. */
   
void nmsim_elem_net_stats_initialize(nmsim_elem_net_stats_t *ntS);
  /* Clears all the fields of {ntS} in preparation for {nmsim_elem_net_stats_accumulate}. */
   
void nmsim_elem_net_stats_accumulate_V_age_M_H
  ( nmsim_elem_net_stats_t *ntS,
    nmsim_time_t t,
    nmsim_elem_neuron_count_t nne,
    double V[],                /* Neuron potentials at time {t}. */
    nmsim_step_count_t age[],  /* Ages of neurons at time {t}. */
    double M[],                /* Recharge modulators at time {t}. */
    double H[],                /* Output modulators at time {t}. */
  );
  /* Accumulates statistics of neuron state variables {V,age,M,H}
    for one time value {t}.  See 
    
    Specifically, accumulates {V[ine],age[ine],M[ine],H[ine]}
    into {ntS.V,ntS.age,ntS.M,ntS.H}, respectively.
    
    Considers only neurons in the range {ntS.ilo..ntS.ihi} and only if
    {t} is in {ntS.tlo..ntS.thi}. */
    
void nmsim_elem_net_stats_accumulate_VF_AF_X_I_J
  ( nmsim_elem_net_stats_t *ntS,
    nmsim_time_t t,
    nmsim_elem_neuron_count_t nne,
    double V[],                /* Neuron potentials at time {t}. */
    nmsim_step_count_t age[],  /* Ages of neurons at time {t}. */
    bool_t X[],                /* Firing indicators for step {t} to {t+1}. */
    double I[],                /* External neuron inputs for step {t} to {t+1}. */
    double J[]                 /* Total inputs for step {t} to {t+1}. */
  );
  /* Accumulates statistics variables that refer to one time step from
    {t} to {t+1}.  Specifically, accumulates {V[ine],age[ine],X[ine],I[ine],J[ine]}
    into {ntS.V_fire,ntS.age_fire,ntS.X,ntS.I,ntS.J}. The first two are accumulated
    only for neurons that fired in that time step.
    
    Considers only neurons in the range {ntS.ilo..ntS.ihi} and only if
    {t} is in {ntS.tlo..ntS.thi-1}. */

void nmsim_elem_net_stats_finalize(nmsim_elem_net_stats_t *ntS);
  /* For each field of of {ntS}, converts the {.avg} subfield from 
    sum of sample values to their average, and the {.dev} subfield
    from sum of squared samples to their standard  deviation.
    
    See {nmsim_stats_finalize} for more details. */

#endif
