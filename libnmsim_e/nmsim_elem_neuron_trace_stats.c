/* See {nmsim_elem_neuron_trace_stats.h} */
/* Last edited on 2020-12-16 00:23:59 by jstolfi */

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
#include <nmsim_elem_neuron_trace.h>
#include <nmsim_elem_neuron_trace_entry.h>
#include <nmsim_elem_net_sim_stats.h>

#include <nmsim_elem_neuron_trace_stats.h>

void nmsim_elem_neuron_trace_stats_compute
  ( nmsim_elem_neuron_trace_t *trne,
    nmsim_elem_net_sim_stats_t *S
  )
  {
    /* Copy data: */
    S->ine_lo = trne->ine;
    S->ine_hi = trne->ine;
    S->t_lo = trne->tlo;
    S->t_hi = trne->thi;
    
    /* Clear statistics: */
    nmsim_elem_net_sim_stats_initialize(S);

    /* Scan trace entries and collect {nvs,min,max} for each parameter.  Accumulate sum in {avg}. */
    nmsim_time_t tlo = trne->tlo;
    nmsim_time_t thi = trne->thi;
    for (nmsim_time_t t = tlo; t <= thi; t++)
      { nmsim_elem_neuron_trace_entry_t *tsk = &(trne->ts[t - tlo]);
        nmsim_stats_accumulate(&(S->V), tsk->V);
        nmsim_stats_accumulate(&(S->age), (double)tsk->age);
        nmsim_stats_accumulate(&(S->M), tsk->M);
        nmsim_stats_accumulate(&(S->H), tsk->H);
        /* The variables {X,I,J} may not be defined for {t=thi}: */
        if (t < thi)
          { nmsim_stats_accumulate(&(S->X), (double)tsk->X);
            nmsim_stats_accumulate(&(S->I), tsk->I);
            nmsim_stats_accumulate(&(S->J), tsk->J);
            if (tsk->X)
              { /* Neuron fired in the next time step. Collect stats: */
                nmsim_stats_accumulate(&(S->V_fire), tsk->V);
                nmsim_stats_accumulate(&(S->age_fire), (double)tsk->age);
              }
          }
      }

    nmsim_elem_net_sim_stats_finalize(S);
    return;
  }
