/* See {nmsim_elem_net_stats.h} */
/* Last edited on 2020-12-15 21:53:54 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <nget.h>
#include <fget.h>
#include <filefmt.h>

#include <nmsim_basic.h>
#include <nmsim_stats.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_net.h>
#include <nmsim_elem_net_entry.h>

#include <nmsim_elem_net_stats.h>

void nmsim_elem_net_stats_write(FILE *wr, nmsim_elem_net_stats_t *trS)
  {
    char *ind1 = "";    /* Indentation for whole summary. */
    char *ind2 = "  ";  /* Indentation for each parameter and state. */
    
    /* Write the file header: */
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_elem_net_stats_FILE_TYPE, nmsim_elem_net_stats_VERSION);
    
    /* Write the neuron index and time range: */
    fprintf(wr, "%sneuron_elem = %d times %ld..%ld\n", ind2, trS->ine, trS->tlo, trS->thi);
    fputs(ind2, wr); nmsim_stats_print(wr, "V",        &(trS->V),        nmsim_write_VIJ_PREC, TRUE,  FALSE, FALSE);
    fputs(ind2, wr); nmsim_stats_print(wr, "V_fire",   &(trS->V_fire),   nmsim_write_VIJ_PREC, TRUE,  FALSE, FALSE);
    fputs(ind2, wr); nmsim_stats_print(wr, "age",      &(trS->age),      0.06, FALSE,  FALSE, FALSE);
    fputs(ind2, wr); nmsim_stats_print(wr, "age_fire", &(trS->age_fire), 0.06, FALSE,  FALSE, FALSE);
    fputs(ind2, wr); nmsim_stats_print(wr, "M",        &(trS->M),        nmsim_write_MH_PREC, FALSE,  TRUE, TRUE);
    fputs(ind2, wr); nmsim_stats_print(wr, "H",        &(trS->H),        nmsim_write_MH_PREC, FALSE,  TRUE, TRUE);
    fputs(ind2, wr); nmsim_stats_print(wr, "X",        &(trS->X),        nmsim_write_rho_PREC, FALSE,  TRUE, TRUE);
    fputs(ind2, wr); nmsim_stats_print(wr, "I",        &(trS->I),        nmsim_write_VIJ_PREC, TRUE,  FALSE, FALSE);
    fputs(ind2, wr); nmsim_stats_print(wr, "J",        &(trS->J),        nmsim_write_VIJ_PREC, TRUE,  FALSE, FALSE);

    /* Write the file footer: */
    fputs(ind1, wr);
    filefmt_write_footer(wr, nmsim_elem_net_stats_FILE_TYPE);

    fflush(wr);
  }

void nmsim_elem_net_stats_initialize(nmsim_elem_net_stats_t *ntS)
  { 
    nmsim_stats_initialize(&(nts->V));
    nmsim_stats_initialize(&(nts->age));
    nmsim_stats_initialize(&(nts->M));
    nmsim_stats_initialize(&(nts->H));
    nmsim_stats_initialize(&(nts->X));
    nmsim_stats_initialize(&(nts->I));
    nmsim_stats_initialize(&(nts->J));

    nmsim_stats_initialize(&(nts->V_fire));
    nmsim_stats_initialize(&(nts->age_fire));
  }
   
void nmsim_elem_net_stats_accumulate_V_age_M_H
  ( nmsim_elem_net_stats_t *ntS,
    nmsim_time_t t,
    nmsim_elem_neuron_count_t nne,
    double V[],                /* Neuron potentials at time {t}. */
    nmsim_step_count_t age[],  /* Ages of neurons at time {t}. */
    double M[],                /* Recharge modulators at time {t}. */
    double H[],                /* Output modulators at time {t}. */
  )
  {
    if ((t >= ntS->tlo) && (t <= ntS->thi))
      { nmsim_elem_neuron_ix_t ilo = imax(0, ntS->ilo);
        nmsim_elem_neuron_ix_t ihi = imin(ntS->ihi, nne-1);
        for (nmsim_elem_neuron_ix_t ine = ilo; ine <= ihi; ine++)
          { nmsim_stats_accumulate(&(nts->V), V[ine]);
            nmsim_stats_accumulate(&(nts->age), age[ine]);
            nmsim_stats_accumulate(&(nts->M), M[ine]);
            nmsim_stats_accumulate(&(nts->H), H[ine]);
          }
      }
  }
    
void nmsim_elem_net_stats_accumulate_VF_AF_X_I_J
  ( nmsim_elem_net_stats_t *ntS,
    nmsim_time_t t,
    nmsim_elem_neuron_count_t nne,
    double V[],                /* Neuron potentials at time {t}. */
    nmsim_step_count_t age[],  /* Ages of neurons at time {t}. */
    bool_t X[],                /* Firing indicators for step {t} to {t+1}. */
    double I[],                /* External neuron inputs for step {t} to {t+1}. */
    double J[]                 /* Total inputs for step {t} to {t+1}. */
  )
  {
    if ((t >= ntS->tlo) && (t <= ntS->thi))
      { nmsim_elem_neuron_ix_t ilo = imax(0, ntS->ilo);
        nmsim_elem_neuron_ix_t ihi = imin(ntS->ihi, nne-1);
        for (nmsim_elem_neuron_ix_t ine = ilo; ine <= ihi; ine++)
          { nmsim_stats_accumulate(&(nts->X), X[ine]);
            nmsim_stats_accumulate(&(nts->I), I[ine]);
            nmsim_stats_accumulate(&(nts->J), J[ine]);
            if (X[ine])
              { nmsim_stats_accumulate(&(nts->V_fire));
                nmsim_stats_accumulate(&(nts->age_fire));
              }
          }
      }
  }

void nmsim_elem_net_stats_finalize(nmsim_elem_net_stats_t *ntS)
  { 
    nmsim_stats_finalize(&(nts->V));
    nmsim_stats_finalize(&(nts->age));
    nmsim_stats_finalize(&(nts->M));
    nmsim_stats_finalize(&(nts->H));
    nmsim_stats_finalize(&(nts->X));
    nmsim_stats_finalize(&(nts->I));
    nmsim_stats_finalize(&(nts->J));

    nmsim_stats_finalize(&(nts->V_fire));
    nmsim_stats_finalize(&(nts->age_fire));
  }
