/* See {nmsim_elem_net_trace.h} */
/* Last edited on 2020-12-16 00:23:03 by jstolfi */

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
#include <jsfile.h>
#include <jsrandom.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_neuron_trace.h>
#include <nmsim_elem_neuron_trace_stats.h>
#include <nmsim_elem_neuron_trace_entry.h>
 
#include <nmsim_elem_net_trace.h>

nmsim_elem_net_trace_t *nmsim_elem_net_trace_new 
  ( nmsim_time_t tlo,                 /* Discrete time of first recorded states. */
    nmsim_time_t thi,                 /* Discrete time of last recorded states. */
    nmsim_elem_neuron_count_t tne      /* Number of neurons to trace. */  
  )
  {
    demand((tlo >= nmsim_time_MIN) && (tlo <= thi) && (thi <= nmsim_time_MAX), "invalid time range");
    nmsim_elem_net_trace_t *etrace = notnull(malloc(sizeof(nmsim_elem_net_trace_t)), "no mem");
    nmsim_elem_neuron_trace_t **trne = notnull(malloc(tne*sizeof(nmsim_elem_neuron_trace_t*)), "no mem");
    for (nmsim_elem_neuron_ix_t k = 0; k < tne; k++) { trne[k] = NULL; }
    (*etrace) = (nmsim_elem_net_trace_t) { .tlo = tlo, .thi = thi, .tne = tne, .trne = trne };
    return etrace;
  }

void nmsim_elem_net_trace_free(nmsim_elem_net_trace_t *etrace)
  { if (etrace != NULL)
      { if (etrace->tne != 0) 
          { for (nmsim_elem_neuron_ix_t k = 0; k < etrace->tne; k++)
              { nmsim_elem_neuron_trace_t *trnek = etrace->trne[k];
                if (trnek != NULL) { nmsim_elem_neuron_trace_free(trnek); }
              }
            free(etrace->trne);
          }
        free(etrace);
      }
  }

void nmsim_elem_net_trace_set_V_age_M_H
  ( nmsim_elem_net_trace_t *etrace,
    nmsim_time_t t,                /* Discrete time. */
    nmsim_elem_neuron_count_t nne, /* Number of neurons in net. */
    double V[],                    /* Membrane potential of each neuron at {t}. */
    nmsim_step_count_t age[],      /* Firing age of each neuron at time {t}. */
    double M[],                    /* Recharge modulator of each neuron at time {t}. */
    double H[]                     /* Output modulator of each neuron at time {t}. */
  )
  { if (etrace != NULL)
      { for (nmsim_elem_neuron_ix_t k = 0; k < etrace->tne; k++)
          { nmsim_elem_neuron_trace_t *trnek = etrace->trne[k];
            if ((trnek != NULL) && (t >= trnek->tlo) && (t <= trnek->thi))
              { nmsim_elem_neuron_ix_t ine = trnek->ine; /* Index of monitored neuron. */
                assert((ine >= 0) && (ine < nne));
                nmsim_elem_neuron_trace_set_V_age_M_H(trnek, t, V[ine], age[ine], M[ine], H[ine]);
              }
          }
      }
  }

void nmsim_elem_net_trace_set_X_I_J
  ( nmsim_elem_net_trace_t *etrace,
    nmsim_time_t t,                /* Time at start of step. */
    nmsim_elem_neuron_count_t nne, /* Number of neurons in net. */
    bool_t X[],                    /* Firing indicators between {t} and {t+1}. */
    double I[],                    /* External input of each neuron from {t} to {t+1} (mV). */
    double J[]                     /* Total input of each neuron from {t} to {t+1} (mV). */
  )
  { if ((etrace != NULL) && (t >= etrace->tlo) && (t <= etrace->thi))
      { for (nmsim_elem_neuron_ix_t k = 0; k < etrace->tne; k++)
          { nmsim_elem_neuron_trace_t *trnek = etrace->trne[k];
            if ((trnek != NULL) && (t >= trnek->tlo) && (t <= trnek->thi))
              { nmsim_elem_neuron_ix_t ine = trnek->ine; /* Index of monitored neuron. */
                assert((ine >= 0) && (ine < nne));
                nmsim_elem_neuron_trace_set_X_I_J(trnek, t, X[ine], I[ine], J[ine]);
              }
          }
      }
  }

void nmsim_elem_net_trace_write(char *prefix, nmsim_elem_net_trace_t *etrace)
  {
    for (nmsim_elem_neuron_ix_t k = 0; k < etrace->tne; k++)
      { nmsim_elem_neuron_trace_t *trnek = etrace->trne[k];
        if (trnek != NULL)
          { nmsim_elem_neuron_ix_t ine = trnek->ine;
            nmsim_time_t tlo = trnek->tlo;
            nmsim_time_t thi = trnek->thi;
            { /* Trace data by time step: */
              char *fname = NULL;
              asprintf(&fname, "%s_n%010d_trace.txt", prefix, ine);
              FILE *wr = open_write(fname, TRUE);
              nmsim_elem_neuron_trace_write(wr, trnek);
              fclose(wr);
              free(fname);
            }
            { /* Statistical summary: */
              char *fname = NULL;
              asprintf(&fname, "%s_n%010d_stats.txt", prefix, ine);
              FILE *wr = open_write(fname, TRUE);
              nmsim_elem_net_sim_stats_t *trS = nmsim_elem_net_sim_stats_new(ine, ine, tlo, thi);
              nmsim_elem_neuron_trace_stats_compute(trnek, trS);
              nmsim_elem_net_sim_stats_write(wr, trS);
              fclose(wr);
              free(trS);
              free(fname); 
            }
          }
      }
  }

nmsim_elem_net_trace_t *nmsim_elem_net_trace_throw
  ( nmsim_elem_neuron_count_t nne,  /* Number of neurons in network. */
    nmsim_time_t tlo,              /* Discrete time of first neuron trace entries. */
    nmsim_time_t thi,              /* Discrete time of last neuron trace entries. */
    nmsim_elem_neuron_count_t tne   /* Number of neurons to monitor. */
  )
  { demand((tne >= 0) & (tne <= nne), "invalid monitored neuron count");
  
    /* Create trace structure: */
    nmsim_elem_net_trace_t *etrace = nmsim_elem_net_trace_new(tlo, thi, tne);
    
    /* Select {tne} random neurons out of {nne}, put in {perm[0..tne-1]}: */
    uint64_t *perm = uint64_choose(nne, tne, NULL);
    
    /* Create traces for those neurons with random time intervals: */
    for (nmsim_elem_neuron_ix_t k = 0; k < tne; k++)
      { nmsim_elem_neuron_ix_t inek = (nmsim_elem_neuron_ix_t)perm[k];
        nmsim_time_t tlok, thik;
        nmsim_throw_time_range(tlo, thi, &tlok, &thik);
        nmsim_elem_neuron_trace_t *trnek = nmsim_elem_neuron_trace_new(inek, tlok, thik);
        etrace->trne[k] = trnek;
      }
    return etrace;
  }
    
  
