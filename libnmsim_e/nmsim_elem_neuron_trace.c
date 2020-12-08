/* See {nmsim_elem_neuron_trace.h} */
/* Last edited on 2020-12-07 14:44:31 by jstolfi */

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
  ( nmsim_elem_neuron_ix_t ine,       /* Index of the neuron in the network. */
    nmsim_time_t tlo,                /* Discrete time of first recorded state. */
    nmsim_time_t thi                 /* Discrete time of last recorded state. */
  )
  { 
    demand((ine >= 0) && (ine <= nmsim_elem_neuron_ix_MAX), "invalid neuron index");
    demand((tlo >= nmsim_time_MIN) && (tlo <= thi) && (thi <= nmsim_time_MAX), "invalid time range");
    
    /* Allocate and initialize the array of full states: */
    nmsim_step_count_t nns = thi-tlo+1; /* Number of discrete times monitored. */
    demand(nns <= nmsim_elem_neuron_trace_entries_MAX, "too many states");
    nmsim_elem_neuron_trace_entry_t *ts = notnull(malloc(nns*sizeof(nmsim_elem_neuron_trace_entry_t)), "no mem");
    for (nmsim_time_t t = tlo; t <= thi; t++)
       { nmsim_elem_neuron_trace_entry_clear(&(ts[t - tlo])); }
       
    /* Allocate and fill the head record: */
    nmsim_elem_neuron_trace_t *trne = notnull(malloc(sizeof(nmsim_elem_neuron_trace_t)), "no mem");
    (*trne) = (nmsim_elem_neuron_trace_t)
      { .ine = ine, .tlo = tlo, .thi = thi, .ts = ts };
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
    nmsim_time_t t,              /* Discrete time. */
    double V,                    /* Membrane potential of each neuron at {t}. */
    nmsim_step_count_t age,      /* Firing age of each neuron at time {t}. */
    double M,                    /* Recharge modulator of each neuron at time {t}. */
    double H                     /* Output modulator of each neuron at time {t}. */
  )
  {
    if ((t >= trne->tlo) && (t <= trne->thi))
      { /* Get the trace entry for time {t}. */
        nmsim_elem_neuron_trace_entry_t *tst = &(trne->ts[t - trne->tlo]); 
        tst->V = V; tst->age = age;
        tst->M = M; tst->H = H;
      }
  }
  
void nmsim_elem_neuron_trace_set_X_I_J
  ( nmsim_elem_neuron_trace_t *trne,
    nmsim_time_t t,              /* Time at start of step. */
    bool_t X,                    /* Firing indicators between {t} and {t+1}. */
    double I,                    /* External input of each neuron from {t} to {t+1} (mV). */
    double J                     /* Total input of each neuron from {t} to {t+1} (mV). */
  )  
  {
    if ((t >= trne->tlo) && (t <= trne->thi))
      { /* Get the trace entry for time {t}. */
        nmsim_elem_neuron_trace_entry_t *tst = &(trne->ts[t - trne->tlo]); 
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
    fprintf(wr, "%sneuron_elem = %d\n", ind2, trne->ine);
    
    /* Trim leading and trailing undefined entries: */
    nmsim_time_t tlo = trne->tlo;
    nmsim_time_t thi = trne->thi;
    while ((tlo <= thi) && nmsim_elem_neuron_trace_entry_is_undefined(&(trne->ts[tlo - trne->tlo]))) { tlo++; }
    while ((tlo <= thi) && nmsim_elem_neuron_trace_entry_is_undefined(&(trne->ts[thi - trne->tlo]))) { thi--; }
    if (tlo > thi) { tlo = 0; thi = -1; }
    
    /* Write the time range: */
    fprintf(wr, "%sinitial_time = %ld\n", ind2, tlo);
    fprintf(wr, "%sfinal_time = %ld\n", ind2, thi);
    
    /* Write the defined entries: */
    for (nmsim_time_t t = tlo; t <= thi; t++)
      { nmsim_elem_neuron_trace_entry_t *tst = &(trne->ts[t - trne->tlo]);
        fputs(ind2, wr);
        fprintf(wr, "%10ld ", t);
        nmsim_elem_neuron_trace_entry_write(wr, tst); 
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
    nmsim_elem_neuron_ix_t ine = (nmsim_elem_neuron_ix_t)
      nmsim_read_int64_param(rd, "neuron_elem", 0, nmsim_elem_neuron_ix_MAX);
   
    /* Read the time range: */
    nmsim_time_t tlo = (nmsim_time_t)
      nmsim_read_int64_param(rd, "initial_time", nmsim_time_MIN, nmsim_time_MAX);
    nmsim_time_t thi = (nmsim_time_t)
      nmsim_read_int64_param(rd, "final_time", tlo, nmsim_time_MAX);
    demand(tlo <= thi, "invalid time range");

    /* Allocate trace: */
    nmsim_elem_neuron_trace_t *trne = nmsim_elem_neuron_trace_new(ine, tlo, thi); 
    
    /* Read the states: */
     for (nmsim_time_t t = tlo; t <= thi; t++)
      { nmsim_elem_neuron_trace_entry_t *tst = &(trne->ts[t - trne->tlo]);
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
    nmsim_compare_int64_param("neuron_elem_index",  trne_read->ine, trne_orig->ine);
    /* Compute union {tlo..thi} of the time spans: */
    fprintf(stderr, "  read time range     %ld ... %ld\n", trne_read->tlo, trne_read->thi);
    fprintf(stderr, "  original time range %ld ... %ld\n", trne_orig->tlo, trne_orig->thi);
    nmsim_time_t tlo, thi;
    nmsim_int64_range_unite(trne_read->tlo, trne_read->thi, trne_orig->tlo, trne_orig->thi, &tlo, &thi);
    fprintf(stderr, "  combined time range %ld ... %ld\n", tlo, thi);
    for (nmsim_time_t t = tlo; t <= thi; t++)
      { /* fprintf(stderr, "    time %ld entries [%ld] and [%ld]\n", t, t-trne_read->tlo,  t-trne_orig->tlo); */
        /* Get the enties for time {t}: */
        nmsim_elem_neuron_trace_entry_t *tst_read = 
          ((t >= trne_read->tlo) && (t <= trne_read->thi) ? &(trne_read->ts[t - trne_read->tlo]) : NULL);
        nmsim_elem_neuron_trace_entry_t *tst_orig = 
          ((t >= trne_orig->tlo) && (t <= trne_orig->thi) ? &(trne_orig->ts[t - trne_orig->tlo]) : NULL);
        nmsim_elem_neuron_trace_entry_compare(tst_read, tst_orig);
      }
  }
