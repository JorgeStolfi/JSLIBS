#ifndef nmsim_basic_H
#define nmsim_basic_H
 
/* Basic types for neuromat network simulation. */
/* Last edited on 2020-12-07 15:31:57 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

/* COUNTS AND INDICES FOR INDIVIDUAL NEURONS */

typedef int32_t nmsim_elem_neuron_count_t;
  /* A count of individual neurons in a network (or in a list of neurons).
    We cannot define it as {int64_t} because we cannot allocate arrays
    with more than {2^30} elements or so. */
  
typedef int32_t nmsim_elem_neuron_ix_t;
  /* Index of a neuron in a network (or in a list of neurons). */
  
#define nmsim_elem_neuron_count_MAX (((uint32_t)1) << 30)
#define nmsim_elem_neuron_ix_MAX (nmsim_elem_neuron_count_MAX - 1)
  /* Maximum neuron count and index. */

/* COUNTS AND INDICES FOR INDIVIDUAL SYNAPSES */
  
typedef int32_t nmsim_elem_synapse_count_t;
  /* A count of individual synapses in a network (or in a list of synapses). */
  
typedef int32_t nmsim_elem_synapse_ix_t;
  /* Index of a synapse in a network (or in a list of synapses). */

#define nmsim_elem_synapse_count_MAX (((uint32_t)1) << 30)
#define nmsim_elem_synapse_ix_MAX (nmsim_elem_synapse_count_MAX - 1)
  /* Maximum synapses count and index in a network. */

/* COUNTS AND INDICES FOR NEURON GROUPS */
  
/* A /neuron group/ (or /neuron population/, or /population/ for short) is a 
   subset of the neurons, with parameters and connections
   that are the homogeneous (or drawn from the same distribution). */

typedef int32_t nmsim_group_neuron_count_t;
  /* A count of neuron groups. */

typedef int32_t nmsim_group_neuron_ix_t;
  /* Index of a neuron group. */
  
#define nmsim_group_neuron_count_MAX (((uint32_t)1) << 24)
#define nmsim_group_neuron_ix_MAX (nmsim_group_neuron_count_MAX - 1)
  /* Maximum neuron group count and index. */
  
/* COUNTS AND INDICES FOR NEURON COHORTS */
  
/* A /neuron cohort/ is subset of the neurons in a neuron group
  which have the same firing age (or firing age exceeding some cutoff). */

typedef int32_t nmsim_cohort_neuron_count_t;
  /* A count of neuron cohorts. */

typedef int32_t nmsim_cohort_neuron_ix_t;
  /* Index of a neuron cohort. */
  
#define nmsim_cohort_neuron_count_MAX (((uint32_t)1) << 14)
#define nmsim_cohort_neuron_ix_MAX (nmsim_cohort_neuron_count_MAX - 1)
  /* Maximum neuron cohort count and index. */

/* COUNTS AND INDICES FOR SYNAPTIC GROUPS */
   
/* A /synapse group/ (or /synaptic bundle/, or /bundle/ for short) is a 
   subset of the synapses with parameters that are homogeneous
   (or drawn from the same distribution), connecting neurons in
   one specific neuron group to another specific neuron group. */
 
typedef int32_t nmsim_group_synapse_count_t;
  /* A count of synapse groups. */

typedef int32_t nmsim_group_synapse_ix_t;
  /* Index of a synapse group. */
  
#define nmsim_group_synapse_count_MAX (((int32_t)1) << 24)
#define nmsim_group_synapse_ix_MAX (nmsim_group_synapse_count_MAX)
  /* Maximum synapse group count and index. */

/* COUNTS AND INDICES FOR NEURON CLASSES */

/* A /neuron class/ is a subset of neurons with the same parameters (or the
  same statistical distribution of parameters).  It comprises one or more
  neuron populations, all with the same parameters. */

typedef nmsim_group_neuron_count_t nmsim_class_neuron_count_t;
  /* A count of neuron classes. */

typedef nmsim_group_neuron_ix_t nmsim_class_neuron_ix_t;
  /* Index of a neuron class. */

#define nmsim_class_neuron_count_MAX (((int32_t)1) << 20)
#define nmsim_class_neuron_ix_MAX (nmsim_class_neuron_count_MAX - 1)
  /* Maximum neuron class count and index. */
  
/* COUNTS AND INDICES FOR SYNAPTIC CLASSES */

/* A /synapse class/ is a set of synapses with the same parameters (or the 
  same statistical distribution of parameters).  It comprises one or
  more synaptic groups, all with the same parameters.  The neurons connected
  by a given synapse class are NOT restricted to praticular neuron classes. */

typedef nmsim_group_synapse_count_t nmsim_class_synapse_count_t;
  /* A count of synapse classes. */

typedef nmsim_group_synapse_ix_t nmsim_class_synapse_ix_t;
  /* Index of a synapse class. */

#define nmsim_class_synapse_count_MAX (((int32_t)1) << 20)
#define nmsim_class_synapse_ix_MAX (nmsim_class_synapse_count_MAX - 1)
  /* Maximum synapse class count and index. */

/* DISCRETE TIMES AND COUNTS AND INDICES FOR TIME STEPS */
    
typedef int64_t nmsim_step_count_t;
  /* A count of time steps, e.g. the firing age of a neuron. */
  
typedef int64_t nmsim_step_ix_t;
  /* Index of a time step in a simulation.  The step with index {i}
    goes from discrete time {i} to discrete time {i+1}. */

typedef int64_t nmsim_time_t;
  /* A discrete sampling time, expressed as the number of time steps
    since some arbitrary origin time. */

#define nmsim_step_count_MAX (((int64_t)1) << 40)
#define nmsim_step_ix_MAX (nmsim_step_count_MAX-1)
#define nmsim_time_MAX (nmsim_step_count_MAX)
#define nmsim_time_MIN (-(nmsim_step_count_MAX))
  /* Maximum count of time steps and minimum maximum sampling time
    index. Note that in a simulation with {N} time steps the discrete
    times go from 0 to {N}. */

/* RANGE CLIPPING */

void nmsim_int64_range_clip
  ( int64_t alo,
    int64_t ahi, 
    int64_t blo, 
    int64_t bhi, 
    int64_t *cloP, 
    int64_t *chiP,
    char *itname
  );
void nmsim_int32_range_clip
  ( int32_t alo,
    int32_t ahi, 
    int32_t blo, 
    int32_t bhi, 
    int32_t *cloP, 
    int32_t *chiP,
    char *itname
  );
  /* Computes the intersection of two integer ranges {alo..ahi} and {blo..bhi}, 
    and returns the low and high ends of the result in {*cloP} and {*chiP}.
    If the intersection is empty, returns {*cloP=0} and {*chiP=-1}.
    
    If {itname} is not {NULL}, assumes that it is the name of an item in those
    ranges (e.g. "index", "time", etc.) and prints warnings if {alo..ahi} is not fully inside
    the range {blo..bhi}, or is completely outside it. */

void nmsim_int64_range_unite
  ( int64_t alo,
    int64_t ahi, 
    int64_t blo, 
    int64_t bhi, 
    int64_t *uloP, 
    int64_t *uhiP
  );
void nmsim_int32_range_unite
  ( int32_t alo,
    int32_t ahi, 
    int32_t blo, 
    int32_t bhi, 
    int32_t *uloP, 
    int32_t *uhiP
  );
  /* Computes the union of two integer ranges {alo..ahi} and {blo..bhi}, 
    and returns the low and high ends of the result in {*uloP} and {*uhiP}. */

/* OTHER UTILITIES */
  
double nmsim_throw_double(double vlo, double vhi);
  /* Generates a random double-precision value in the range {[vlo _ vhi]},
    rounded to a proper number of decimal digits. */

double nmsim_select_rounding_mod(double vlo, double vhi);
  /* Selects a suitable rounding modulus for rounding a random number between
    {vlo} and {vhi}. */
    
void nmsim_throw_time_range(nmsim_time_t tlo, nmsim_time_t thi, nmsim_time_t *tloP, nmsim_time_t *thiP);
  /* Generates a random time interval {*tloP ..* thiP} roughly contained in {tlo .. thi}, but 
    may extrapolate somewhat from that range. */

/* MISCELLANEOUS */

#define INF INFINITY

#endif
