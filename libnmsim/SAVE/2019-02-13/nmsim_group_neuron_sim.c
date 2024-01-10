/* See {nmsim_group_neuron_sim.h} */
/* Last edited on 2019-02-13 17:36:16 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron_state.h>
#include <nmsim_cohort_neuron_state.h>
#include <nmsim_cohort_neuron_sim.h>

#include <nmsim_group_neuron_sim.h>

void nmsim_group_neuron_sim_step
  ( nmsim_group_neuron_state_t *pso,
    double DV_avg,
    double DV_dev,
    nmsim_class_neuron_t *nclass,
    nmsim_group_neuron_state_t *psn
  )
  {
    /* The state of the new cohort of age zero: */
    nmsim_cohort_neuron_state_t csn_fire;
    nmsim_cohort_neuron_state_clear(&csn_fire);
    /* The state of the new cohort of age {nc}: */
    nmsim_cohort_neuron_state_t csn_lump;

    int32_t nc = pso->nc;
    int32_t k;
    for (k = nc; k > 0; k--)
      { /* Evolve cohort {S[k,t]} to time {t+1}: */
        nmsim_cohort_neuron_state_t *cso = &(pso->cs[k]);
        nmsim_cohort_neuron_state_t csn_k_fire;
        nmsim_cohort_neuron_state_t csn_k_fail;
        nmsim_cohort_neuron_sim_step
	  ( cso, DV_avg, DV_dev, 
	    nclass, 
	    &csn_k_fire, &csn_k_fail
	  );

        /* Merge the part of the cohort that fired into {csn_fire}: */
        nmsim_cohort_neuron_state_merge(&csn_k_fire, &csn_fire, nclass);

        /* Store or merge the part of the cohort that did not fire: */
        if (k == nc)
          { /* Save the new cohort {nc+1} to merge later: */
            csn_lump = csn_k_fail; 
	  }
        else 
          { nmsim_cohort_neuron_state_t *csn = &(psn->cs[k+1]); 
            (*csn) = csn_k_fail;
            if (k == nc-1) { nmsim_cohort_neuron_state_merge(&csn_lump, csn, nclass); }
	  }
      }
   
   }
