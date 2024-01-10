#ifndef nmsim_group_neuron_state_H
#define nmsim_group_neuron_state_H

/* Mean-field state and tools for group-level network simulation. */
/* Last edited on 2019-07-02 15:10:52 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_group_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_class_neuron.h>
#include <nmsim_firing_func.h>
#include <nmsim_cohort_neuron_state.h>

typedef struct nmsim_group_neuron_state_t
  { nmsim_cohort_neuron_count_t nc;     /* Number of cohorts individually represented. */
    nmsim_cohort_neuron_state_vec_t cs; /* State of neurons in each cohort. */
  } nmsim_group_neuron_state_t;
  /* Stochastic state of a homogeneous population of neurons in
    group-level simulation. 
    
    The neurons with firing age {k} have state {cs.e[k]}, for {k} in
    {0..nc-1}. All neurons with age {nc} or greater are lumped
    together and assumed to have state {cs.e[nc].}
    The allocated size of the vector, {cs.ne} will be at least {nc+1} elements.
    
    Specifically, the parameters {cs.e[k].V_avg} and {cs.e[k].V_dev}
    specify the potential distribution of the neurons (assumed to be
    Gaussian) of the cohort that has firing age {k} at some time {t}.
    the field {cs[k].eta} is the probability that a neuron of the group
    is in in that cohort. */
    
nmsim_group_neuron_state_t nmsim_group_neuron_state_make(void);
  /* Returns a neuron group state record, initially with zero cohorts 
    and all value fields set to {NAN}. */
    
void nmsim_group_neuron_state_free(nmsim_group_neuron_state_t *st);
  /* Reclaims the space of the internal tables of {st}, such as 
    {st->cs}. Does NOT reclaim the record {st} itself. */
    
void nmsim_group_neuron_state_throw
  ( nmsim_group_neuron_state_t *st,
    nmsim_class_neuron_t *nclass,
    nmsim_cohort_neuron_count_t nc
  );
  /* Generates a "random" group-level state {*st} for a neuron group
    belonging to neuron class {nclass}.
    Allocates the cohort vector of {st}, which should be {NULL} on entry,
    for {nc+1} cohorts.
    
    In the current version, the cohorts will be all empty except the last one, 
    which will have all neurons and a potential precisely equal to the
    resting potential {V_B} of the neurons.  That is equivalent to assuming
    that the neurons have not fired nor received any input for a long time. */
    
void nmsim_group_neuron_state_copy(nmsim_group_neuron_state_t *src, nmsim_group_neuron_state_t *dst);
  /* Copies the state {src} into the state {dst}. The cohort vector of {dst}
    is allocated or re-allocated as needed. */
    
void nmsim_group_neuron_state_write(FILE *wr, nmsim_group_neuron_state_t *st);
  /* Writes the state {st} of a neuron group at some time {t} to file {wr},
    in the format described by {nmsim_group_neuron_state_read_INFO} below.
    The written data is NOT preceded or terminated by a newline. */

void nmsim_group_neuron_state_read
  ( FILE *rd,
    nmsim_group_neuron_state_t *st
  );
  /* Reads from file {rd} the state of a neuron group at some time {t}
    and saves it into {*st}. Assumes the format described in
    {nmsim_group_neuron_state_read_INFO} below. */
    
#define nmsim_group_neuron_state_read_INFO \
  "Data for the state of a neuron group at some time {t} consists of the following numbers, all on the same line:\n" \
  "\n" \
  "    \"{nc} {COHORTS}\"\n" \
  "\n" \
  "  where {nc} is the number of cohorts stored.\n" \
  "\n" \
  "  Then there follows the data for the {nc+1} represented" \
  " coords, each in the format described below.\n" \
  "\n" \
  "COHORT STATE\n" \
  "  " nmsim_cohort_neuron_state_read_INFO ""

/* DEBUGGING */

void nmsim_group_neuron_state_compare
  ( nmsim_group_neuron_state_t *st_read, 
    nmsim_group_neuron_state_t *st_orig
  );
  /* Compares the neuron group state {st_read} read from a file
    with the expected class parameters {st_orig}.  Aborts if any parameter
    does not match (within roundoff tolerance). */

#endif
