#ifndef nmsim_basic_H
#define nmsim_basic_H
 
/* Basic types for neuromat network simulation. */
/* Last edited on 2019-02-13 14:43:05 by jstolfi */

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

/* NEURON CLASSES */

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
  
/* SYNAPTIC CLASSES */

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

/* TIME STEPS */
    
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
  /* Maximum count of time steps and maximum sampling time index.
    Note that in a simulation with {N} time steps
    the discrete times go from 0 to {N}, not {N-1}. */

/* UTILITIES */

void check_int64_value(char *name, int64_t v, int64_t vmin, int64_t vmax);
  /* Aborts if the value {v} is not in the range {vmin..vmax}. */

void check_double_value(char *name, double v, double vmin, double vmax);
  /* Aborts if the value {v} is not in the range {[vmin _ vmax]}. */

int64_t read_int64_value(FILE *rd, char *name, int64_t vmin, int64_t vmax);
  /* Reads an integer from {rd}, which must be on the current line. 
    Aborts if the value read is not in the range {vmin..vmax}.  */

double read_double_value(FILE *rd, char *name, double vmin, double vmax);
  /* Reads a {double} value from {rd}, which must be on the current line. 
    Aborts if the value read is not in the range {[vmin _ vmax]}.  */

void write_double_value(FILE *wr, double v, char *fmt);
  /* Writes {v} using the given format (of the 'f' or 'e' kind),
    but without leading blanks and trailing fractional zeros. */

void write_int64_parm(FILE *wr, char *pref, char *name, int64_t v, char *fmt);
  /* Writes a line "{pref}{name} = {v}" using the given format (of the 'd' kind). */

void write_double_parm(FILE *wr, char *pref, char *name, double v, char *fmt);
  /* Writes a line "{pref}{name} = {v}" using the given format (of the 'f' or 'e' kind).
    Removes leading blanks and trailing fractional zeros. */

int64_t read_int64_parm(FILE *rd, char *name, int64_t vmin, int64_t vmax);
  /* Reads a line \"{name} = {value}\", including the end-of-line.
    Aborts if the {value} is not an integer in the range {vmin..vmax}..  */

double read_double_parm(FILE *rd, char *name, double vmin, double vmax);
  /* Reads a line \"{name} = {value}\", including the end-of-line.
    Aborts if the {value} is in the range {[vmin _ vmax]}. */

void compare_int64_parm(char *name, int64_t vr, int64_t ve);
  /* Compares a parameter value {vr} read from a file
    with the expected value {ve}. Aborts if different. */
    
void compare_double_parm(char *name, double vr, double ve, double tol);
  /* Compares a parameter value {vr} read from a file
    with the expected value {ve}. Aborts if abs difference
    is more than {tol}.  Requires identity if either is 0.0
    or 1.0. */
    
double throw_double(double vlo, double vhi);
  /* Generates a random double-precision value in the range {[vlo _ vhi]},
    rounded to a proper number of decimal digits. */

double select_rounding_mod(double vlo, double vhi);
  /* Selects a suitable rounding modulus for rounding a random number between
    {vlo} and {vhi}. */
    
#endif
