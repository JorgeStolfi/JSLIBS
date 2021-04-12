/* See {nmsim_elem_net_sim_stats.h} */
/* Last edited on 2020-12-25 11:50:43 by jstolfi */

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
#include <jsmath.h>

#include <nmsim_basic.h>
#include <nmsim_stats.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_net.h>

#include <nmsim_elem_net_sim_stats.h>

nmsim_elem_net_sim_stats_t *nmsim_elem_net_sim_stats_new
  ( nmsim_elem_neuron_ix_t ineLo,    /* Index of the first neuron of the set. */
    nmsim_elem_neuron_ix_t ineHi,    /* Index of the last neuron of the set. */
    nmsim_time_t tLo,                /* Discrete time of first state to consider. */
    nmsim_time_t tHi                 /* Discrete time of last state to consider. */
  )
  {
    nmsim_elem_net_sim_stats_t *S = notnull(malloc(sizeof(nmsim_elem_net_sim_stats_t)), "no mem");
    S->ineLo = ineLo;
    S->ineHi = ineHi;
    S->tLo = tLo;
    S->tHi = tHi;
    return S;
  }

void nmsim_elem_net_sim_stats_write(FILE *wr, nmsim_elem_net_sim_stats_t *S)
  {
    char *ind1 = "";    /* Indentation for whole summary. */
    char *ind2 = "  ";  /* Indentation for each parameter and state. */
    
    /* Write the file header: */
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_elem_net_sim_stats_FILE_TYPE, nmsim_elem_net_sim_stats_VERSION);
    
    /* Write the neuron index and time range: */
    fprintf(wr, "%sfirst_neuron = %d\n", ind2, S->ineLo);
    fprintf(wr, "%slast_neuron = %d\n", ind2, S->ineHi);
    fprintf(wr, "%sinitial_time = %ld\n", ind2, S->tLo);
    fprintf(wr, "%sfinal_time = %ld\n", ind2, S->tHi);
    fputs(ind2, wr); nmsim_stats_print(wr, "V",        &(S->V),        nmsim_write_VIJ_PREC, TRUE,  FALSE, FALSE);
    fputs(ind2, wr); nmsim_stats_print(wr, "V_fire",   &(S->V_fire),   nmsim_write_VIJ_PREC, TRUE,  FALSE, FALSE);
    fputs(ind2, wr); nmsim_stats_print(wr, "age",      &(S->age),      nmsim_write_age_PREC, FALSE, FALSE, FALSE);
    fputs(ind2, wr); nmsim_stats_print(wr, "age_fire", &(S->age_fire), nmsim_write_age_PREC, FALSE, FALSE, FALSE);
    fputs(ind2, wr); nmsim_stats_print(wr, "M",        &(S->M),        nmsim_write_MH_PREC,  FALSE, TRUE, TRUE);
    fputs(ind2, wr); nmsim_stats_print(wr, "H",        &(S->H),        nmsim_write_MH_PREC,  FALSE, TRUE, TRUE);
    fputs(ind2, wr); nmsim_stats_print(wr, "X",        &(S->X),        nmsim_write_rho_PREC, FALSE, TRUE, TRUE);
    fputs(ind2, wr); nmsim_stats_print(wr, "I",        &(S->I),        nmsim_write_VIJ_PREC, TRUE,  FALSE, FALSE);
    fputs(ind2, wr); nmsim_stats_print(wr, "S",        &(S->S),        nmsim_write_VIJ_PREC, TRUE,  FALSE, FALSE);
    fputs(ind2, wr); nmsim_stats_print(wr, "J",        &(S->J),        nmsim_write_VIJ_PREC, TRUE,  FALSE, FALSE);

    /* Write the file footer: */
    fputs(ind1, wr);
    filefmt_write_footer(wr, nmsim_elem_net_sim_stats_FILE_TYPE);

    fflush(wr);
  }

void nmsim_elem_net_sim_stats_initialize(nmsim_elem_net_sim_stats_t *S)
  { 
    nmsim_stats_initialize(&(S->V));
    nmsim_stats_initialize(&(S->age));
    nmsim_stats_initialize(&(S->M));
    nmsim_stats_initialize(&(S->H));
    nmsim_stats_initialize(&(S->X));
    nmsim_stats_initialize(&(S->I));
    nmsim_stats_initialize(&(S->S));
    nmsim_stats_initialize(&(S->J));

    nmsim_stats_initialize(&(S->V_fire));
    nmsim_stats_initialize(&(S->age_fire));
  }
   
void nmsim_elem_net_sim_stats_accumulate_V_age_M_H
  ( nmsim_elem_net_sim_stats_t *S,
    nmsim_time_t t,
    nmsim_elem_neuron_count_t nne,
    double V[],                /* Neuron potentials at time {t}. */
    nmsim_step_count_t age[],  /* Ages of neurons at time {t}. */
    double M[],                /* Recharge modulators at time {t}. */
    double H[]                 /* Output modulators at time {t}. */
  )
  {
    /* Accept data for {t} in {tLo..tHi}: */
    if ((t >= S->tLo) && (t <= S->tHi))
      { nmsim_elem_neuron_ix_t ineLo = (nmsim_elem_neuron_ix_t)imax(0, S->ineLo);
        nmsim_elem_neuron_ix_t ineHi = (nmsim_elem_neuron_ix_t)imin(S->ineHi, nne-1);
        for (nmsim_elem_neuron_ix_t ine = ineLo; ine <= ineHi; ine++)
          { nmsim_stats_accumulate(&(S->V), V[ine]);
            nmsim_stats_accumulate(&(S->age), (double)age[ine]);
            nmsim_stats_accumulate(&(S->M), M[ine]);
            nmsim_stats_accumulate(&(S->H), H[ine]);
          }
      }
  }
    
void nmsim_elem_net_sim_stats_accumulate_VF_AF_X_I_J
  ( nmsim_elem_net_sim_stats_t *S,
    nmsim_time_t t,
    nmsim_elem_neuron_count_t nne,
    double V[],                /* Neuron potentials at time {t}. */
    nmsim_step_count_t age[],  /* Ages of neurons at time {t}. */
    bool_t X[],                /* Firing indicators for step {t} to {t+1}. */
    double I[],                /* External neuron inputs for step {t} to {t+1}. */
    double J[]                 /* Total inputs for step {t} to {t+1}. */
  )
  {
    /* Accept data for {t} in {tLo..tHi-1} (excluding {tHi}): */
    if ((t >= S->tLo) && (t < S->tHi))
      { nmsim_elem_neuron_ix_t ineLo = (nmsim_elem_neuron_ix_t)S->ineLo;
        nmsim_elem_neuron_ix_t ineHi = (nmsim_elem_neuron_ix_t)S->ineHi;
        assert((ineLo >= 0) && (ineLo <= ineHi) && (ineHi < nne));
        for (nmsim_elem_neuron_ix_t ine = ineLo; ine <= ineHi; ine++)
          { nmsim_stats_accumulate(&(S->X), X[ine]);
            nmsim_stats_accumulate(&(S->I), I[ine]);
            nmsim_stats_accumulate(&(S->S), J[ine]-I[ine]);
            nmsim_stats_accumulate(&(S->J), J[ine]);
            if (X[ine])
              { nmsim_stats_accumulate(&(S->V_fire), V[ine]);
                nmsim_stats_accumulate(&(S->age_fire), (double)age[ine]);
              }
          }
      }
  }

void nmsim_elem_net_sim_stats_finalize(nmsim_elem_net_sim_stats_t *S)
  { 
    nmsim_stats_finalize(&(S->V));
    nmsim_stats_finalize(&(S->age));
    nmsim_stats_finalize(&(S->M));
    nmsim_stats_finalize(&(S->H));
    nmsim_stats_finalize(&(S->X));
    nmsim_stats_finalize(&(S->I));
    nmsim_stats_finalize(&(S->S));
    nmsim_stats_finalize(&(S->J));

    nmsim_stats_finalize(&(S->V_fire));
    nmsim_stats_finalize(&(S->age_fire));
  }
