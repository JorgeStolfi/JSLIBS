/* See {nmsim_elem_neuron_trace.h} */
/* Last edited on 2020-12-17 09:24:32 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <nget.h>
#include <fget.h>
#include <filefmt.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_neuron_trace_entry.h>

#include <nmsim_elem_neuron_trace.h>

nmsim_elem_neuron_trace_t *nmsim_elem_neuron_trace_new
  ( nmsim_elem_neuron_ix_t ineLo,   /* Index of the first neuron of the set to trace. */
    nmsim_elem_neuron_ix_t ineHi,   /* Index of the last neuron of the set to trace. */
    nmsim_time_t tLo,                /* Discrete time of first recorded state. */
    nmsim_time_t tHi                 /* Discrete time of last recorded state. */
  )
  { 
    demand((ineLo >= 0) && (ineLo <= ineHi) && (ineHi <= nmsim_elem_neuron_ix_MAX), "invalid neuron index");
    demand((tLo >= nmsim_time_MIN) && (tLo <= tHi) && (tHi <= nmsim_time_MAX), "invalid time range");
    
    /* Allocate and initialize the array of full states: */
    nmsim_step_count_t nns = tHi-tLo+1; /* Number of discrete times monitored. */
    demand(nns <= nmsim_elem_neuron_trace_entries_MAX, "too many states");
    nmsim_elem_neuron_trace_entry_t *ts = notnull(malloc(nns*sizeof(nmsim_elem_neuron_trace_entry_t)), "no mem");
    for (nmsim_time_t t = tLo; t <= tHi; t++) { nmsim_elem_neuron_trace_entry_clear(&(ts[t - tLo])); }
       
    /* Allocate and fill the head record: */
    nmsim_elem_neuron_trace_t *trne = notnull(malloc(sizeof(nmsim_elem_neuron_trace_t)), "no mem");
    (*trne) = (nmsim_elem_neuron_trace_t)
      { .ineLo = ineLo, .ineHi = ineHi, .tLo = tLo, .tHi = tHi, .ts = ts };
    return trne;
  }

void nmsim_elem_neuron_trace_free(nmsim_elem_neuron_trace_t *trne)    
  {
    if (trne != NULL)
      { if (trne->ts != NULL) { free(trne->ts); }
        free(trne);
      }
  }
  
void nmsim_elem_neuron_trace_set_V_age_M_H
  ( nmsim_elem_neuron_trace_t *trne,
    nmsim_time_t t,  /* Discrete time. */
    double V,        /* Membrane potential of each neuron at {t}. */
    double age,      /* Firing age of each neuron at time {t}. */
    double M,        /* Recharge modulator of each neuron at time {t}. */
    double H         /* Output modulator of each neuron at time {t}. */
  )
  {
    if (trne == NULL) { return; }
    if ((t >= trne->tLo) && (t <= trne->tHi))
      { /* Get the trace entry for time {t}. */
        nmsim_elem_neuron_trace_entry_t *tst = &(trne->ts[t - trne->tLo]); 
        tst->V = V; tst->age = (double)age;
        tst->M = M; tst->H = H;
      }
  }
  
void nmsim_elem_neuron_trace_set_X_I_J
  ( nmsim_elem_neuron_trace_t *trne,
    nmsim_time_t t,   /* Time at start of step. */
    double X,         /* Firing indicators between {t} and {t+1}. */
    double I,         /* External input of each neuron from {t} to {t+1} (mV). */
    double J          /* Total input of each neuron from {t} to {t+1} (mV). */
  )  
  {
    if (trne == NULL) { return; }
    if ((t >= trne->tLo) && (t < trne->tHi))
      { /* Get the trace entry for time {t}. */
        nmsim_elem_neuron_trace_entry_t *tst = &(trne->ts[t - trne->tLo]); 
        tst->X = X; tst->I = I; tst->J = J;
      }
  }

void nmsim_elem_neuron_trace_write(FILE *wr, nmsim_elem_neuron_trace_t *trne)
  {
    char *ind1 = "";    /* Indentation for whole trace. */
    char *ind2 = "  ";  /* Indentation for each parameter and state. */
    
    /* Write the file header: */
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_elem_neuron_trace_FILE_TYPE, nmsim_elem_neuron_trace_VERSION);
    
    /* Write the neuron index: */
    fprintf(wr, "%sfirst_neuron = %d\n", ind2, trne->ineLo);
    fprintf(wr, "%slast_neuron = %d\n", ind2, trne->ineHi);
    
    /* Trim leading and trailing undefined entries: */
    nmsim_time_t tLo = trne->tLo;
    nmsim_time_t tHi = trne->tHi;
    while ((tLo <= tHi) && nmsim_elem_neuron_trace_entry_is_undefined(&(trne->ts[tLo - trne->tLo]))) { tLo++; }
    while ((tLo <= tHi) && nmsim_elem_neuron_trace_entry_is_undefined(&(trne->ts[tHi - trne->tLo]))) { tHi--; }
    if (tLo > tHi) { tLo = 0; tHi = -1; }
    
    /* Write the time range: */
    fprintf(wr, "%sinitial_time = %ld\n", ind2, tLo);
    fprintf(wr, "%sfinal_time = %ld\n", ind2, tHi);
    
    /* Write the defined entries: */
    bool_t single = (trne->ineLo == trne->ineHi);
    for (nmsim_time_t t = tLo; t <= tHi; t++)
      { nmsim_elem_neuron_trace_entry_t *tst = &(trne->ts[t - trne->tLo]);
        fputs(ind2, wr);
        fprintf(wr, "%10ld ", t);
        nmsim_elem_neuron_trace_entry_write(wr, tst, single); 
        fputs("\n", wr);
      }

    /* Write the file footer: */
    fputs(ind1, wr);
    filefmt_write_footer(wr, nmsim_elem_neuron_trace_FILE_TYPE);

    fflush(wr);
  }

nmsim_elem_neuron_trace_t *nmsim_elem_neuron_trace_read(FILE *rd)
  { 
    /* Read header line: */
    filefmt_read_header(rd, nmsim_elem_neuron_trace_FILE_TYPE, nmsim_elem_neuron_trace_VERSION);
    
    /* Read the neuron indices: */
    nmsim_elem_neuron_ix_t ineLo = (nmsim_elem_neuron_ix_t)
      nmsim_read_int64_param(rd, "first_neuron", 0, nmsim_elem_neuron_ix_MAX);
    nmsim_elem_neuron_ix_t ineHi = (nmsim_elem_neuron_ix_t)
      nmsim_read_int64_param(rd, "last_neuron", ineLo, nmsim_elem_neuron_ix_MAX);
   
    /* Read the time range: */
    nmsim_time_t tLo = (nmsim_time_t)
      nmsim_read_int64_param(rd, "initial_time", nmsim_time_MIN, nmsim_time_MAX);
    nmsim_time_t tHi = (nmsim_time_t)
      nmsim_read_int64_param(rd, "final_time", tLo, nmsim_time_MAX);
    demand(tLo <= tHi, "invalid time range");

    /* Allocate trace: */
    nmsim_elem_neuron_trace_t *trne = nmsim_elem_neuron_trace_new(ineLo, ineHi, tLo, tHi); 
    
    /* Read the states: */
     for (nmsim_time_t t = tLo; t <= tHi; t++)
      { nmsim_elem_neuron_trace_entry_t *tst = &(trne->ts[t - trne->tLo]);
        fget_skip_formatting_chars(rd);
        /* Check time: */
        (void)nmsim_read_int64_value(rd, "time", t, t);
        /* Read the state: */
        nmsim_elem_neuron_trace_entry_read(rd, tst); 
        fget_eol(rd);
      }

    /* Read footer line: */
    filefmt_read_footer(rd, nmsim_elem_neuron_trace_FILE_TYPE);

    return trne;
  }

void nmsim_elem_neuron_trace_compare
  ( nmsim_elem_neuron_trace_t *trne_read, 
    nmsim_elem_neuron_trace_t *trne_orig
  )
  {
    nmsim_compare_int64_param("first_neuron",  trne_read->ineLo, trne_orig->ineLo);
    nmsim_compare_int64_param("first_neuron",  trne_read->ineHi, trne_orig->ineHi);
    /* Compute union {tLo..tHi} of the time spans: */
    fprintf(stderr, "  read time range     %ld ... %ld\n", trne_read->tLo, trne_read->tHi);
    fprintf(stderr, "  original time range %ld ... %ld\n", trne_orig->tLo, trne_orig->tHi);
    nmsim_time_t tLo, tHi;
    nmsim_int64_range_unite(trne_read->tLo, trne_read->tHi, trne_orig->tLo, trne_orig->tHi, &tLo, &tHi);
    fprintf(stderr, "  combined time range %ld ... %ld\n", tLo, tHi);
    for (nmsim_time_t t = tLo; t <= tHi; t++)
      { /* fprintf(stderr, "    time %ld entries [%ld] and [%ld]\n", t, t-trne_read->tLo,  t-trne_orig->tLo); */
        /* Get the enties for time {t}: */
        nmsim_elem_neuron_trace_entry_t *tst_read = 
          ((t >= trne_read->tLo) && (t <= trne_read->tHi) ? &(trne_read->ts[t - trne_read->tLo]) : NULL);
        nmsim_elem_neuron_trace_entry_t *tst_orig = 
          ((t >= trne_orig->tLo) && (t <= trne_orig->tHi) ? &(trne_orig->ts[t - trne_orig->tLo]) : NULL);
        nmsim_elem_neuron_trace_entry_compare(tst_read, tst_orig);
      }
  }
