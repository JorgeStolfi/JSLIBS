/* See {nmsim_elem_neuron_trace_stats.h} */
/* Last edited on 2020-12-09 00:21:58 by jstolfi */

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

#include <nmsim_elem_neuron_trace_stats.h>

void nmsim_elem_neuron_trace_stats_write(FILE *wr, nmsim_elem_neuron_trace_stats_t *trS)
  {
    char *ind1 = "";    /* Indentation for whole summary. */
    char *ind2 = "  ";  /* Indentation for each parameter and state. */
    
    /* Write the file header: */
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_elem_neuron_trace_stats_FILE_TYPE, nmsim_elem_neuron_trace_stats_VERSION);
    
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
    filefmt_write_footer(wr, nmsim_elem_neuron_trace_stats_FILE_TYPE);

    fflush(wr);
  }

void nmsim_elem_neuron_trace_stats_compute
  ( nmsim_elem_neuron_trace_t *trne,
    nmsim_elem_neuron_trace_stats_t *trS
  )
  {
    auto void clear(nmsim_stats_t *S);
      /* Initializes {S.nvs,S.min,S.max} to proper values, sets {S.avg,S.dev} to zero. */

    auto void accum1(double v, nmsim_stats_t *S);
      /* If sample {v} is neither {NAN} nor {±INF}, updates {S.nvs,S.min,S.max} 
        to account for it, and adds {v} {S.avg}.  Does not change {S.dev}. */
        
    auto void compavg(nmsim_stats_t *S);
      /* Converts {S.avg} from sum to average by dividing {S.nvs} into it.
        Sets it to {NAN} if {S.nvs} is zero. */

    auto void accum2(double v, nmsim_stats_t *S);
      /* If sample {v} is neither {NAN} nor {±INF}, adds {(v - S.avg)^2} to {S.dev}. */
        
    auto void compdev(nmsim_stats_t *S);
      /* Converts {S.dev} from sum of sqare differences from {S.avg}
        to deviation, by dividing {S.nvs-1} into it.
        Sets it to {NAN} if {S.nvs} is less than 2. */

    /* Copy data: */
    trS->ine = trne->ine;
    trS->tlo = trne->tlo;
    trS->thi = trne->thi;
    
    /* Clear statistics: */
    clear(&(trS->V));
    clear(&(trS->age));
    clear(&(trS->M));
    clear(&(trS->H));
    clear(&(trS->X));
    clear(&(trS->I));
    clear(&(trS->J));
    clear(&(trS->V_fire));
    clear(&(trS->age_fire));

    /* Scan trace entries and collect {nvs,min,max} for each parameter.  Accumulate sum in {avg}. */
    nmsim_time_t tlo = trne->tlo;
    nmsim_time_t thi = trne->thi;
    for (nmsim_time_t t = tlo; t <= thi; t++)
      { nmsim_elem_neuron_trace_entry_t *tsk = &(trne->ts[t - tlo]);
        accum1(tsk->V, &(trS->V));
        accum1((double)tsk->age, &(trS->age));
        accum1(tsk->M, &(trS->M));
        accum1(tsk->H, &(trS->H));
        accum1((double)tsk->X, &(trS->X));
        accum1(tsk->I, &(trS->I));
        accum1(tsk->J, &(trS->J));
        if (tsk->X)
          { /* Neuron fired in the next time step. Collect stats: */
            accum1(tsk->V, &(trS->V_fire));
            accum1((double)tsk->age, &(trS->age_fire));
          }
      }
    
    /* Compute averages: */
    compavg(&(trS->V));
    compavg(&(trS->age));
    compavg(&(trS->M));
    compavg(&(trS->H));
    compavg(&(trS->X));
    compavg(&(trS->I));
    compavg(&(trS->J));
    compavg(&(trS->V_fire));
    compavg(&(trS->age_fire));
    
    /* Scan trace entries again and accumulate sum of squared diferences from average in {dev}. */
    for (nmsim_time_t t = tlo; t <= thi; t++)
      { nmsim_elem_neuron_trace_entry_t *tsk = &(trne->ts[t - tlo]);
        accum2(tsk->V, &(trS->V));
        accum2((double)tsk->age, &(trS->age));
        accum2(tsk->M, &(trS->M));
        accum2(tsk->H, &(trS->H));
        accum2((double)tsk->X, &(trS->X));
        accum2(tsk->I, &(trS->I));
        accum2(tsk->J, &(trS->J));
        if (tsk->X)
          { /* Neuron fired in the next time step. Collect stats: */
            accum2(tsk->V, &(trS->V_fire));
            accum2((double)tsk->age, &(trS->age_fire));
          }
      }

    /* Compute deviations: */
    compdev(&(trS->V));
    compdev(&(trS->age));
    compdev(&(trS->M));
    compdev(&(trS->H));
    compdev(&(trS->X));
    compdev(&(trS->I));
    compdev(&(trS->J));
    compdev(&(trS->V_fire));
    compdev(&(trS->age_fire));
    
    return;
    
    void clear(nmsim_stats_t *S)
      { S->nvs = 0;
        S->min = +INF;
        S->max = -INF;
        S->avg = 0.0;
        S->dev = 0.0;
      }
      
    void accum1(double v, nmsim_stats_t *S)
      {
        if ((! isnan(v)) && (! isinf(v)))
          { S->nvs++;
            S->min = fmin(S->min, v);
            S->max = fmax(S->max, v);
            S->avg += v;
          }
      }
      
    void compavg(nmsim_stats_t *S)
      {
        S->avg /= ((double)S->nvs); 
      }
      
    void accum2(double v, nmsim_stats_t *S)
      {
        if ((! isnan(v)) && (! isinf(v)))
          { double d = v - S->avg;
            S->dev += d*d;
          }
      }
      
    void compdev(nmsim_stats_t *S)
      {
        if (S->nvs < 2)
          { S->dev = NAN; }
        else
          { S->dev = sqrt(S->dev/((double)(S->nvs - 1))); }
      }
      
  }
