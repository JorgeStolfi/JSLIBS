/* See {nmsim_group_net_trace.h} */
/* Last edited on 2019-02-13 18:26:42 by jstolfi */

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
#include <nmsim_group_neuron_trace.h>
 
#include <nmsim_group_net_trace.h>

nmsim_group_net_trace_t *nmsim_group_net_trace_new 
  ( nmsim_time_t tini,                   /* Discrete time of first recorded states. */
    nmsim_time_t tfin,                   /* Discrete time of last recorded states. */
    nmsim_group_neuron_count_t nng_tr     /* Number of neuron groups to monitor. */
  )
  {
    demand((nng_tr >= 0) && (nng_tr <= nmsim_group_neuron_count_MAX), "invalid neuron group count");
    demand((tini >= 0) && (tini <= tfin) && (tfin <= nmsim_time_MAX), "invalid time range");
    
    nmsim_group_neuron_trace_t **trng = notnull(malloc(nng_tr*sizeof(nmsim_group_neuron_trace_t*)), "no mem");
    for (nmsim_group_neuron_ix_t i = 0; i < nng_tr; i++) { trng[i] = NULL; }
    nmsim_group_net_trace_t *etrace = notnull(malloc(sizeof(nmsim_group_net_trace_t)), "no mem");
    (*etrace) = (nmsim_group_net_trace_t)
      { .tini = tini, .tfin = tfin,
        .nng = nng_tr, .trng = trng
      };
    return etrace;
  }

void nmsim_group_net_trace_free(nmsim_group_net_trace_t *etrace)
  { if (etrace != NULL)
      { if (etrace->trng != NULL) 
          { for (nmsim_group_neuron_ix_t i = 0; i < etrace->nng; i++)
              { nmsim_group_neuron_trace_t *trng = etrace->trng[i];
                if (trng != NULL) 
                  { nmsim_group_neuron_trace_free(trng); }
              }
            free(etrace->trng);
          }
        free(etrace);
      }
  }

void nmsim_group_net_trace_write(char *prefix, nmsim_group_net_trace_t *etrace)
  {
    for (nmsim_group_neuron_ix_t i = 0; i < etrace->nng; i++)
      { nmsim_group_neuron_trace_t *trng = etrace->trng[i];
        if (trng != NULL)
          { char *fname = NULL;
            asprintf(&fname, "%s_g%010d.txt", prefix, trng->ing);
            FILE *wr = open_write(fname, TRUE);
            nmsim_group_neuron_trace_write(wr, trng);
            fclose(wr);
            free(fname);
          }
      }
  }
