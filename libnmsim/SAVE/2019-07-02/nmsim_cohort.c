/* See {nmsim_cohort.h} */
/* Last edited on 2016-12-08 11:51:13 by jstolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>

#include <nmsim_firing_func.h>
#include <nmsim_neuron.h>
#include <nmsim_cohort.h>

void nmsim_cohort_mf_state_set
  ( nmsim_cohort_mf_state_t *cs, 
    double V_avg,
    double V_dev,
    double G,
    double H,
    double eta,
    nmsim_neuron_parms_t *parms
  )
  {
    double rho;
    nmsim_firing_func_eval(parms->Phi, G*V_avg, &rho, NULL);
    (*cs) = (nmsim_cohort_mf_state_t)
      { .V_avg = V_avg, .V_dev = V_dev,
        .G = G, .H = H,
        .eta = eta, .rho = rho
      };
  }

void nmsim_cohort_mf_state_clear(nmsim_cohort_mf_state_t *cs)
  {
    (*cs) = (nmsim_cohort_mf_state_t)
      { .V_avg = 0.0, .V_dev = 0.0,
        .G = 1.0, .H = 1.0,
        .eta = 0.0, .rho = 0.0
      };
  }

void nmsim_cohort_mf_state_merge
  ( nmsim_cohort_mf_state_t *csa, 
    nmsim_cohort_mf_state_t *cst,
    nmsim_neuron_parms_t *parms
  )
  {
    double ea = csa->eta;
    double et = cst->eta;
    assert(ea >= 0);
    if (ea == 0) { return; }

    assert(et >= 0);
    if (et == 0)
      { /* Just replace: */
        (*cst) = (*csa); 
      }
    else
      { /* Get the main parameters: */
        double avga = csa->V_avg;
        double deva = csa->V_dev;
        double avgt = cst->V_avg;
        double devt = cst->V_dev;
        /* Total fraction of the two cohorts in population: */
        double eta_new = ea + et;
        /* Compute average potential of merged populations: */
        double avg_new = (ea*avga + et*avgt)/eta_new; 
        /* Compute the variance of the new population: */
        double da = avga - avg_new;
        double dt = avgt - avg_new;
        double dev_new = sqrt((ea*(da*da + deva*deva) + et*(dt*dt + devt*devt))/eta_new);
        /* Gains (for now, to be deleted): */
        double G_new = (ea*(csa->G) + et*(cst->G))/eta_new;
        double H_new = (ea*(csa->H) + et*(cst->H))/eta_new;
        /* Store results: */
        nmsim_cohort_mf_state_set(cst, avg_new, dev_new, G_new, H_new, eta_new, parms);
      }        
  }

void nmsim_cohort_mf_evolve
  ( nmsim_cohort_mf_state_t *cso, 
    double DV_avg,
    double DV_dev,
    nmsim_neuron_parms_t *parms, 
    nmsim_cohort_mf_state_t *csn_fire,
    nmsim_cohort_mf_state_t *csn_fail
  )
  {
    /* Grab some variables for convenience: */
    double V_avg = cso->V_avg; /* Mean potential of cohort. */
    double V_dev = cso->V_dev; /* Potential deviation in cohort. */
    double G = cso->G;         /* Input gain factor of cohort. */
    double H = cso->H;         /* Output gain factor of cohort. */
    double eta = cso->eta;     /* Fraction of population in cohort. */
    double rho = cso->rho;     /* Fraction of cohort that is about to fire. */

    double mu_V = parms->mu_V;  /* Potential recovery factor. */
    double mu_G = parms->mu_G;  /* Input gain recovery factor. */
    double mu_H = parms->mu_H;  /* Output gain recovery factor. */

    /* Store the new state of the firing neurons in {csn_fire}: */
    nmsim_cohort_mf_state_set
      ( csn_fire,
        parms->V_R, 0.0, parms->G_R, parms->H_R, eta*rho,
        parms
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
    nmsim_cohort_mf_state_set
      ( csn_fail,
        V_avg_new, V_dev_new, G_new, H_new, eta*(1.0 - rho),
        parms
      );
  }

