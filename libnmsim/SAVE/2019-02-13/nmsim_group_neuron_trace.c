/* See {nmsim_group_neuron_trace.h} */
/* Last edited on 2019-02-13 18:25:36 by jstolfi */

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
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_neuron_state.h>

#include <nmsim_group_neuron_trace.h>

nmsim_group_neuron_trace_t *nmsim_group_neuron_trace_new
  ( nmsim_group_neuron_ix_t ing,      /* Index of the neuron group's group in the network. */
    nmsim_class_neuron_ix_t inc,      /* Index of the neuron group's class in the network. */
    nmsim_time_t tini,                /* Discrete time of first recorded state. */
    nmsim_time_t tfin                 /* Discrete time of last recorded state. */
  )
  { 
    demand((ing >= 0) && (ing <= nmsim_group_neuron_ix_MAX), "invalid neuron group group");
    demand((inc >= 0) && (inc <= nmsim_class_neuron_ix_MAX), "invalid neuron class");
    demand((tini >= 0) && (tini <= tfin) && (tfin <= nmsim_time_MAX), "invalid time range");
    
    nmsim_step_count_t nns = tfin-tini+1; /* Number of discrete times monitored. */
    demand(nns <= nmsim_group_neuron_trace_states_MAX, "too many states");
    
    nmsim_group_neuron_trace_t *trng = notnull(malloc(sizeof(nmsim_group_neuron_trace_t)), "no mem");
    nmsim_group_neuron_state_t *st = notnull(malloc(nns*sizeof(nmsim_group_neuron_state_t)), "no mem");
    (*trng) = (nmsim_group_neuron_trace_t)
      { .ing = ing, .inc = inc,
        .tini = tini, .tfin = tfin,
        .st = st
      };
    return trng;
  }

void nmsim_group_neuron_trace_free(nmsim_group_neuron_trace_t *trng)    
  {
    if (trng != NULL)
      { if (trng->st != NULL) 
          { nmsim_time_t tini = trng->tini;
            nmsim_time_t tfin = trng->tfin;
            for (nmsim_time_t t = tini; t <= tfin; t++)
              { nmsim_group_neuron_state_t *st = &(trng->st[t-tini]);
                nmsim_group_neuron_state_free(st);
              }
            free(trng->st);
          }
        free(trng);
      }
  }
  
void nmsim_group_neuron_trace_write(FILE *wr, nmsim_group_neuron_trace_t *trng)
  {
    char *ind1 = "";    /* Indentation for whole trace. */
    char *ind2 = "  ";  /* Indentation for each parameter and state. */
    
    /* Write the file header: */
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_group_neuron_trace_FILE_TYPE, nmsim_group_neuron_trace_VERSION);
    
    /* Write the neuron group indices: */
    fprintf(wr, "%s_group = %d\n", ind2, trng->ing);
    fprintf(wr, "%sneuron_class = %d\n", ind2, trng->inc);
    
    /* Write the time range: */
    fprintf(wr, "%sinitial_time = %ld\n", ind2, trng->tini);
    fprintf(wr, "%sfinal_time = %ld\n", ind2, trng->tfin);
    
    /* Write the states: */
    for (nmsim_time_t t = trng->tini; t < trng->tfin; t++)
      { nmsim_group_neuron_state_t *st = &(trng->st[t-trng->tini]);
        fputs(ind2, wr);
        nmsim_group_neuron_state_write(wr, t, st); 
        fputs("\n", wr);
      }

    /* Write the file footer: */
    fputs(ind1, wr);
    filefmt_write_footer(wr, nmsim_group_neuron_trace_FILE_TYPE);

    fflush(wr);
  }

nmsim_group_neuron_trace_t *nmsim_group_neuron_trace_read(FILE *rd)
  { 
    /* Read header line: */
    filefmt_read_header(rd, nmsim_group_neuron_trace_FILE_TYPE, nmsim_group_neuron_trace_VERSION);
    
    /* Read the neuron group indices: */
    nmsim_group_neuron_ix_t ing = (nmsim_group_neuron_ix_t)
      read_int64_parm(rd, "neuron_group", 0, nmsim_group_neuron_ix_MAX);
    nmsim_class_neuron_ix_t inc = (nmsim_class_neuron_ix_t)
      read_int64_parm(rd, "neuron_class", 0, nmsim_class_neuron_ix_MAX);
    
    /* Read the time range: */
    nmsim_time_t tini = (nmsim_time_t)
      read_int64_parm(rd, "initial_time", 0, nmsim_time_MAX);
    nmsim_time_t tfin = (nmsim_time_t)
      read_int64_parm(rd, "final_time", 0, nmsim_time_MAX);
    demand(tini <= tfin, "invalid time range");

    /* Allocate trace: */
    nmsim_group_neuron_trace_t *trng = nmsim_group_neuron_trace_new(ing, inc, tini, tfin); 
    
    /* Read the states: */
     for (nmsim_time_t t = tini; t < tfin; t++)
      { nmsim_group_neuron_state_t *st = &(trng->st[t-trng->tini]);
        nmsim_group_neuron_state_read(rd, t, st); 
        fget_eol(rd);
      }

    /* Read footer line: */
    filefmt_read_footer(rd, nmsim_group_neuron_trace_FILE_TYPE);

    return trng;
  }

