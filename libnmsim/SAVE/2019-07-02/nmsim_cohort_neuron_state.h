#ifndef nmsim_cohort_neuron_state_H
#define nmsim_cohort_neuron_state_H
 
/* Mean-field state and evolution of a neuron cohort. */
/* Last edited on 2019-07-02 14:40:26 by jstolfi */

#define _GNU_SOURCE

#include <nmsim_class_neuron.h>
#include <nmsim_firing_func.h>

typedef struct nmsim_cohort_neuron_state_t
  { double V_avg; /* Mean potential. */
    double V_dev; /* Deviation of potential. */
    double M;     /* Recharge modulator. */
    double H;     /* Output synapse strength modulator. */
    double eta;   /* Probability that a neuron of the group is in this cohort. */
  } nmsim_cohort_neuron_state_t;
  /* Stochastic state of a /cohort/, a set of neurons with the 
    same non-negative firing age {k} within some homogeneous population.  
    
    The neurons are assumed to have a normal distribution of potentials,
    with mean {.V_avg} and deviation {.V_dev}, at some time {t}. Their 
    recharge and output strength modulators are all equalto {.M} and {.H}.  The field
    {.eta} is the probability that a neuron of the population is in this cohort.  */

void nmsim_cohort_neuron_state_set
  ( nmsim_cohort_neuron_state_t *cs, 
    double V_avg,
    double V_dev,
    double M,
    double H,
    double eta
  );
  /* Sets the state parameters of a cohort {cs} to potential {V_avg ±
    V_dev}, modulators {M,H}, and cohort pertinence probability {eta}. */
    
void nmsim_cohort_neuron_state_write
  ( FILE *wr,
    nmsim_cohort_neuron_ix_t k,
    nmsim_cohort_neuron_state_t *cs
  );
  /* Writes the state {cs} of a cohort to file {wr}, in the format described by
    {nmsim_cohort_neuron_state_read_INFO} below. The cohort is assumed to have
    the given non-negative age {k}, in the format described by
    {nmsim_cohort_neuron_state_read_INFO} below.
    The data is NOT preceded or terminated by a newline. */

void nmsim_cohort_neuron_state_read
  ( FILE *rd,
    nmsim_cohort_neuron_ix_t k,
    nmsim_cohort_neuron_state_t *cs
  );
  /* Reads from file {rd} the state of a cohort, assumed to be of the
    given non-negative age {k}, and saves it into {*cs}. Assumes the format described
    in {nmsim_cohort_neuron_state_read_INFO} below. Checks whether the
    cohort index in the file is {k}. */
    
#define nmsim_cohort_neuron_state_read_INFO \
  "Data for the state of a cohort at some time {t} consists" \
  " of seven numbers, all on the same line:\n" \
  "\n" \
  "    \"{k} {V_avg} {V_dev} {M} {H} {eta}\"\n" \
  "\n" \
  "  where {k} is the non-negative cohort's index (neuron" \
  " age, i.e. count of full steps since last firing)," \
  " {V_avg} and {V_dev} are the average and deviation of the membrane" \
  " potential for neurons in that cohort, and {eta} is the" \
  " probability that a neuron of the group is in this cohort, all at time {t}."

/* AUXILIARY PROCEDURES */

void nmsim_cohort_neuron_state_clear(nmsim_cohort_neuron_state_t *cs);
  /* Resets the state {cs} to zero probability {cs.eta} and {NAN}
    for the other parameters. */

bool_t nmsim_cohort_neuron_state_merge
  ( nmsim_cohort_neuron_state_t *csa, 
    nmsim_cohort_neuron_state_t *cst
  );
  /* Tries to merge the cohort state {csa} into the cohort state {cst}.
    The total cohort inclusion probabilities {.eta} are added, the mean and deviation
    of the potentials {.V_avg,.V_dev} are combined as appropriate.
    The modulators are averaged.
    
    In particular, if {csa.eta} is zero then {cst} is not changed;
    if {cst.eta} is zero, then {cst} is set to {csa}. 
    
    If the merge succeeds, the procedure returns {TRUE}.
    If the two distributions are too dissimilar to combine
    into a single gaussian, the merge fails: the {*cst} record
    is not changed, and the procedure returns {FALSE}. */

double nmsim_cohort_neuron_state_compute_firing_rate
  ( nmsim_cohort_neuron_state_t *cs, 
    nmsim_firing_func_t *Phi
  );
  /* Computes the firing probability {cs.rho} for cohorts {cs} whose 
    neurons are assumed to have firing funtion {Phi} and a Gaussian potential
    distribution with parameters {cs.V_avg,cs.V_dev}. */

/* DEBUGGING */

void nmsim_cohort_neuron_state_compare
  ( nmsim_cohort_neuron_state_t *cs_read, 
    nmsim_cohort_neuron_state_t *cs_orig
  );
  /* Compares the neuron cohort state {cs_read} read from a file
    with the expected class parameters {cs_orig}.  Aborts if any parameter
    does not match (within roundoff tolerance). */

vec_typedef(nmsim_cohort_neuron_state_vec_t,nmsim_cohort_neuron_state_vec,nmsim_cohort_neuron_state_t);  

#endif
