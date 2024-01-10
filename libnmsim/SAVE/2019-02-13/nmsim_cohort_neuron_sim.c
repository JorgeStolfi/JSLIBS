/* See {nmsim_cohort_neuron_sim.h} */
/* Last edited on 2019-02-13 17:33:54 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <nmsim_class_neuron.h>
#include <nmsim_cohort_neuron_state.h>
#include <nmsim_cohort_neuron_state.h>

#include <nmsim_cohort_neuron_sim.h>

void nmsim_cohort_neuron_sim_step
  ( nmsim_cohort_neuron_state_t *cso, 
    double DV_avg, 
    double DV_dev, 
    nmsim_class_neuron_t *nclass, 
    nmsim_cohort_neuron_state_t *csn_fire, 
    nmsim_cohort_neuron_state_t *csn_fail
  )
  {
    /* Grab some variables for convenience: */
    double V_avg = cso->V_avg; /* Mean potential of cohort. */
    double V_dev = cso->V_dev; /* Potential deviation in cohort. */
    double M = cso->M;         /* Recharge modulator of cohort. */
    double H = cso->H;         /* Output gain of cohort. */
    double eta = cso->eta;     /* Fraction of population in cohort. */
    double rho = cso->rho;     /* Fraction of cohort that is about to fire. */

    double c_B = nclass->c_B;  /* Recharge factor at rest. */
    
    double M_mu = nclass->M_mu;  /* Recovery factor of recharge modulator. */
    double H_mu = nclass->H_mu;  /* Recovery factor of output modulator. */

    /* Store the new state of the firing neuron elems in {csn_fire}: */
    double V_avg_fire = nclass->V_R;
    double V_dev_fire = 0.0;
    nmsim_cohort_neuron_state_set
      ( csn_fire,
        V_avg_fire, V_dev_fire, nclass->M_R, nclass->H_R, eta*rho,
        nclass
      );
   
    /* Compute the average potential of the non-firing neurons: */
    double c = c_B * M; /* Present recharge factor. */
    double V_avg_new = c*V_avg + DV_avg;
    /* Compute the deviation of the potential of those neurons: */
    /* ??? Can we ignore correlations here? ??? */
    double V_dev_new = hypot(c*V_dev, DV_dev);
    /* Compute the new recharge and output modulators of those neurons: */
    double M_new = 1.0 - M_mu*(1.0 - M);
    double H_new = 1.0 - H_mu*(1.0 - H);
    /* Store the new state of the non-firing neurons in {csn_fail}: */
    nmsim_cohort_neuron_state_set
      ( csn_fail,
        V_avg_new, V_dev_new, M_new, H_new, eta*(1.0 - rho),
        nclass
      );
  }
