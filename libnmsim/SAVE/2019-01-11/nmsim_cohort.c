/* See {nmsim_cohort.h} */
/* Last edited on 2019-01-11 18:05:34 by jstolfi */

#define _GNU_SOURCE

#include <nmsim_firing_func.h>
#include <nmsim_neuron_elem.h>
#include <nmsim_cohort.h>

void nmsim_cohort_set
  ( nmsim_cohort_state_t *cs, 
    double V_avg,
    double V_dev,
    double G,
    double H,
    double eta,
    nmsim_parms_t *parm
  )
  {
    double rho;
    nmsim_firing_func_eval(parm->Phi, G*V_avg, &pr, NULL);
    (*cs) = (nmsim_cohort_state_t)
      { .V_avg = V_avg, .V_dev = V_dev,
        .G = G, .H = H,
        .eta = eta, .rho = rho
      };
  }

void nmsim_cohort_evolve
  ( nmsim_cohort_state_t *cso, 
    double DV_avg,
    double DV_dev,
    nmsim_neuron_class_t *parms, 
    nmsim_cohort_state_t *csn_fire,
    nmsim_cohort_state_t *csn_fail
  )
  {
    /* Grab some variables for convenience: */
    double V_avg = cso->V_avg; /* Mean potential of cohort. */
    double V_dev = cso->V_dev; /* Potential deviation in cohort. */
    double G = cso->G;         /* Input gain factor of cohort. */
    double H - cso->H;         /* Output gain factor of cohort. */
    double eta = cso->eta;     /* Fraction of population in cohort. */
    double rho = cso->rho;     /* Fraction of cohort that is about to fire. */

    double mu_V = parm->mu_V;  /* Potential recovery factor. */
    double mu_G = parm->mu_G;  /* Input gain recovery factor. */
    double mu_H = parm->mu_H;  /* Output gain recovery factor. */

    /* Store the new state of the firing neurons in {csn_fire}: */
    nmsim_cohort_set
      ( csn_fire,
        parms->V_R, 0.0, parms->G_R, parms->H_R, eta*rho,
        parm
      );
   
    /* Compute the average potential of the non-firing neurons: */
    double V_avg_new = mu_V*V_avg + DV_avg;
    /* Compute the deviation of the potential of those neurons: */
    /* ??? Can we ignore correlations here? ??? */
    double V_dev_new = hypot(mu_V*V_dev, DV_dev);
    /* Compute the new input and  output gain factors of those neurons: */
    double G_new = 1.0 - mu_G*(1.0 - G);
    double H_new = 1.0 - mu_H*(1.0 - H);
    /* Store the new state of the non-firing neurons in {csn_fail}: */
    nmsim_cohort_set
      ( csn_fail,
        V_new_avg, V_new_dev, G_new, H_new, eta*(1.0 - rho),
        parm
      );
  }

