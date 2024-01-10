#ifndef nmsim_group_neuron_sim_H
#define nmsim_group_neuron_sim_H

/* Mean-field state and tools for group-level network simulation. */
/* Last edited on 2019-02-13 16:40:14 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_cohort_neuron_state.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_neuron_state.h>
#include <nmsim_class_neuron.h>

void nmsim_group_neuron_sim_step
  ( nmsim_group_neuron_state_t *pso,
    double DV_avg,
    double DV_dev,
    nmsim_class_neuron_t *nclass,
    nmsim_group_neuron_state_t *psn
  );
  /* ??? */
    
#endif
