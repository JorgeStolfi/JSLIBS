/* See {nmsim_group_neuron_state.h} */
/* Last edited on 2019-02-13 16:48:21 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


#include <bool.h>
#include <affirm.h>

#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron_state.h>
#include <nmsim_cohort_neuron_state.h>

nmsim_group_neuron_state_t nmsim_group_neuron_state_make(void)
  { nmsim_group_neuron_state_t st = (nmsim_group_neuron_state_t)
      { .nc = 0, .cs = NULL,
        .I_avg = NAN, .I_dev = NAN,
        .X_avg = NAN,
        .J_avg = NAN, .J_dev = NAN
      };
    return st;
  }

void nmsim_group_neuron_state_free(nmsim_group_neuron_state_t *st)
  { 
    if (st->cs != NULL) { free(st->cs); }
  }

void nmsim_group_neuron_state_throw
  ( nmsim_group_neuron_state_t *st, 
    nmsim_class_neuron_t *nclass,
    nmsim_cohort_neuron_count_t nc
  )
  {
    /* Get the neuron's reset potential: */
    double V_R = nclass->V_R;
    
    /* Allocate and initialize the cohort vector: */
    st->nc = nc;
    assert(st->cs == NULL);
    st->cs = notnull(malloc((nc+1)*sizeof(nmsim_cohort_neuron_state_t)), "no mem");
    for (nmsim_cohort_neuron_ix_t k = 0; k <= nc; k++)
      { /* These cohort parameters are irrelevant for {k < nc}. */
        double V_avg = V_R;
        double V_dev = 0.0;
        double M = 1.0;
        double H = 1.0;
        double eta = (k == nc ? 1.0 : 0.0);
        nmsim_cohort_neuron_state_set(&(st->cs[k]), V_avg, V_dev, M, H, eta, nclass);
      }
    st->I_avg = NAN; st->I_dev = NAN;
    st->X_avg = NAN;
    st->J_avg = NAN; st->J_dev = NAN;
    
  }

