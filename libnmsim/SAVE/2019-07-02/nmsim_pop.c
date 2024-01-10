/* See {nmsim_pop.h} */
/* Last edited on 2016-12-08 11:54:30 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_neuron.h>
#include <nmsim_cohort.h>
#include <nmsim_pop.h>

void nmsim_pop_mf_evolve
  ( nmsim_pop_mf_state_t *pso,
    double DV_avg,
    double DV_dev,
    nmsim_neuron_parms_t *parms,
    nmsim_pop_mf_state_t *psn
  )
  {
    /* The state of the new cohort of age zero: */
    nmsim_cohort_mf_state_t csn_fire;
    nmsim_cohort_mf_state_clear(&csn_fire);
    /* The state of the new cohort of age {nc}: */
    nmsim_cohort_mf_state_t csn_lump;

    int32_t nc = pso->nc;
    int32_t k;
    for (k = nc; k > 0; k--)
      { /* Evolve cohort {S[k,t]} to time {t+1}: */
        nmsim_cohort_mf_state_t *cso = &(pso->cs[k]);
        nmsim_cohort_mf_state_t csn_k_fire;
        nmsim_cohort_mf_state_t csn_k_fail;
        nmsim_cohort_mf_evolve
	  ( cso, DV_avg, DV_dev, 
	    parms, 
	    &csn_k_fire, &csn_k_fail
	  );

        /* Merge the part of the cohort that fired into {csn_fire}: */
        nmsim_cohort_mf_state_merge(&csn_k_fire, &csn_fire, parms);

        /* Store or merge the part of the cohort that did not fire: */
        if (k == nc)
          { /* Save the new cohort {nc+1} to merge later: */
            csn_lump = csn_k_fail; 
	  }
        else 
          { nmsim_cohort_mf_state_t *csn = &(psn->cs[k+1]); 
            (*csn) = csn_k_fail;
            if (k == nc-1) { nmsim_cohort_mf_state_merge(&csn_lump, csn, parms); }
	  }
      }
   
   }

