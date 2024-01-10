/* See {nmsim_elem_net_trace.h} */
/* Last edited on 2019-02-13 18:22:57 by jstolfi */

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

#include <nmsim_basic.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_neuron_trace.h>
 
#include <nmsim_elem_net_trace.h>

nmsim_elem_net_trace_t *nmsim_elem_net_trace_new 
  ( nmsim_time_t tini,                   /* Discrete time of first recorded states. */
    nmsim_time_t tfin,                   /* Discrete time of last recorded states. */
    nmsim_elem_neuron_count_t nne_tr     /* Number of neurons to monitor. */
  )
  {
    demand((nne_tr >= 0) && (nne_tr <= nmsim_elem_neuron_count_MAX), "invalid neuron count");
    demand((tini >= 0) && (tini <= tfin) && (tfin <= nmsim_time_MAX), "invalid time range");
    
    nmsim_elem_neuron_trace_t **trne = notnull(malloc(nne_tr*sizeof(nmsim_elem_neuron_trace_t*)), "no mem");
    for (nmsim_elem_neuron_ix_t i = 0; i < nne_tr; i++) { trne[i] = NULL; }
    nmsim_elem_net_trace_t *etrace = notnull(malloc(sizeof(nmsim_elem_net_trace_t)), "no mem");
    (*etrace) = (nmsim_elem_net_trace_t)
      { .tini = tini, .tfin = tfin,
        .nne = nne_tr, .trne = trne
      };
    return etrace;
  }

void nmsim_elem_net_trace_free(nmsim_elem_net_trace_t *etrace)
  { if (etrace != NULL)
      { if (etrace->trne != NULL) 
          { for (nmsim_elem_neuron_ix_t i = 0; i < etrace->nne; i++)
              { nmsim_elem_neuron_trace_t *trne = etrace->trne[i];
                if (trne != NULL) { nmsim_elem_neuron_trace_free(trne); }
              }
            free(etrace->trne);
          }
        free(etrace);
      }
  }

void nmsim_elem_net_trace_set_V_X
  ( nmsim_elem_net_trace_t *etrace,
    nmsim_time_t t,                /* Time at start of step. */
    nmsim_elem_neuron_count_t nne, /* Number of neurons in net. */
    double V[],                    /* Membrane potential of each neuron at {t}. */
    bool_t X[]                     /* Firing indicators between {t} and {t+1}. */
  )
  { if (etrace != NULL)
      { for (nmsim_elem_neuron_ix_t i = 0; i < etrace->nne; i++)
          { nmsim_elem_neuron_trace_t *trne = etrace->trne[i];
            if ((trne != NULL) && (t >= trne->tini) && (t <= trne->tfin))
              { nmsim_elem_neuron_ix_t ine = trne->ine; /* Index of monitored neuron. */
                assert((ine >= 0) && (ine < nne));
                nmsim_elem_neuron_state_t *st = &(trne->st[t - trne->tini]); /* State rec for time {t}. */
                st->V = V[ine];
                st->X = X[ine];
              }
          }
      }
  }

void nmsim_elem_net_trace_set_age_H_I_J
  ( nmsim_elem_net_trace_t *etrace,
    nmsim_time_t t,                /* Time at start of step. */
    nmsim_elem_neuron_count_t nne, /* Number of neurons in net. */
    nmsim_step_count_t age[],      /* Firing age of each neuron at time {t}. */
    double H[],                    /* Output modulator of each neuron at time {t}. */
    double I[],                    /* External input of each neuron from {t} to {t+1}. */
    double J[]                     /* Total input of each neuron from {t} to {t+1}. */
  )
  { if ((etrace != NULL) && (t >= etrace->tini) && (t <= etrace->tfin))
      { for (nmsim_elem_neuron_ix_t i = 0; i < etrace->nne; i++)
          { nmsim_elem_neuron_trace_t *trne = etrace->trne[i];
            if ((trne != NULL) && (t >= trne->tini) && (t <= trne->tfin))
              { nmsim_elem_neuron_ix_t ine = trne->ine; /* Index of monitored neuron. */
                assert((ine >= 0) && (ine < nne));
                nmsim_elem_neuron_state_t *st = &(trne->st[t - trne->tini]); /* State rec for time {t}. */
                st->age = age[ine];
                st->H = H[ine];
                st->I = I[ine];
                st->J = J[ine];
              }
          }
      }
  }

void nmsim_elem_net_trace_set_M
  ( nmsim_elem_net_trace_t *etrace,
    nmsim_time_t t,                /* Time at start of step. */
    nmsim_elem_neuron_count_t nne, /* Number of neurons in net. */
    double M[]                     /* Potential decay modulator of each neuron at time {t}. */
  )
  { if ((etrace != NULL) && (t >= etrace->tini) && (t <= etrace->tfin))
      { for (nmsim_elem_neuron_ix_t i = 0; i < etrace->nne; i++)
          { nmsim_elem_neuron_trace_t *trne = etrace->trne[i];
            if ((trne != NULL) && (t >= trne->tini) && (t <= trne->tfin))
              { nmsim_elem_neuron_ix_t ine = trne->ine; /* Index of monitored neuron. */
                assert((ine >= 0) && (ine < nne));
                nmsim_elem_neuron_state_t *st = &(trne->st[t - trne->tini]); /* State rec for time {t}. */
                st->M = M[ine];
              }
          }
      }
  }

void nmsim_elem_net_trace_write(char *prefix, nmsim_elem_net_trace_t *etrace)
  {
    for (nmsim_elem_neuron_ix_t i = 0; i < etrace->nne; i++)
      { nmsim_elem_neuron_trace_t *trne = etrace->trne[i];
        if (trne != NULL)
          { char *fname = NULL;
            asprintf(&fname, "%s_n%010d.txt", prefix, trne->ine);
            FILE *wr = open_write(fname, TRUE);
            nmsim_elem_neuron_trace_write(wr, trne);
            fclose(wr);
            free(fname);
          }
      }
  }
