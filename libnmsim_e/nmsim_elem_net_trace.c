/* See {nmsim_elem_net_trace.h} */
/* Last edited on 2020-12-17 13:44:05 by jstolfi */

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
#include <jsmath.h>
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

nmsim_elem_net_trace_t *nmsim_elem_net_trace_new(nmsim_elem_neuron_count_t ntr)
  {
    nmsim_elem_net_trace_t *etrace = notnull(malloc(sizeof(nmsim_elem_net_trace_t)), "no mem");
    nmsim_time_t tLo = nmsim_time_MAX;
    nmsim_time_t tHi = nmsim_time_MIN;
    nmsim_elem_neuron_trace_t **trne = notnull(malloc(ntr*sizeof(nmsim_elem_neuron_trace_t*)), "no mem");
    for (nmsim_elem_neuron_ix_t ktr = 0; ktr < ntr; ktr++) { trne[ktr] = NULL; }
    (*etrace) = (nmsim_elem_net_trace_t) { .tLo = tLo, .tHi = tHi, .ntr = ntr, .trne = trne };
    return etrace;
  }

void nmsim_elem_net_trace_set
  ( nmsim_elem_net_trace_t *etrace, 
    int32_t ktr,
    nmsim_elem_neuron_trace_t *trne
  )
  { demand((ktr >= 0) && (ktr < etrace->ntr), "invalid index {ktr}");
    demand(etrace->trne[ktr] == NULL, "entry {ktr} already set");
    etrace->trne[ktr] = trne;
    if (etrace->tLo > etrace->tHi)
      { etrace->tLo = trne->tLo; etrace->tHi = trne->tHi; }
    else
      { etrace->tLo = (nmsim_time_t)imin(etrace->tLo, trne->tLo);
        etrace->tHi = (nmsim_time_t)imin(etrace->tHi, trne->tHi);
      }
  }

void nmsim_elem_net_trace_free(nmsim_elem_net_trace_t *etrace)
  { if (etrace != NULL)
      { if (etrace->ntr != 0) 
          { for (nmsim_elem_neuron_ix_t ktr = 0; ktr < etrace->ntr; ktr++)
              { nmsim_elem_neuron_trace_t *trnek = etrace->trne[ktr];
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
      { for (nmsim_elem_neuron_ix_t ktr = 0; ktr < etrace->ntr; ktr++)
          { nmsim_elem_neuron_trace_t *trnek = etrace->trne[ktr];
            if ((trnek != NULL) && (t >= trnek->tLo) && (t <= trnek->tHi))
              { nmsim_elem_neuron_ix_t ineLo = trnek->ineLo; /* First neuron in monitored subset. */
                nmsim_elem_neuron_ix_t ineHi = trnek->ineHi; /* Last neuron in monitored subset. */
                assert((ineLo >= 0) && (ineLo <= ineHi) && (ineHi < nne));
                double V_tr, age_tr, M_tr, H_tr;
                if (ineLo == ineHi)
                  { nmsim_elem_neuron_ix_t ine = ineLo;
                    V_tr = V[ine]; age_tr = (double)age[ine]; M_tr = M[ine]; H_tr = H[ine];
                  }
                else
                  { V_tr = 0.0; age_tr = 0.0; M_tr = 0.0; H_tr = 0.0;
                    for (nmsim_elem_neuron_ix_t ine = ineLo; ine <= ineHi; ine++)
                      { V_tr += V[ine]; age_tr += (double)age[ine]; M_tr += M[ine]; H_tr += H[ine]; }
                    double dn = (double)(ineHi - ineLo + 1);
                    V_tr /= dn; age_tr /= dn; M_tr /= dn; H_tr /= dn;
                  }
                nmsim_elem_neuron_trace_set_V_age_M_H(trnek, t, V_tr, age_tr, M_tr, H_tr);
                
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
  { if ((etrace != NULL) && (t >= etrace->tLo) && (t <= etrace->tHi))
      { for (nmsim_elem_neuron_ix_t ktr = 0; ktr < etrace->ntr; ktr++)
          { nmsim_elem_neuron_trace_t *trnek = etrace->trne[ktr];
            if ((trnek != NULL) && (t >= trnek->tLo) && (t <= trnek->tHi))
              { nmsim_elem_neuron_ix_t ineLo = trnek->ineLo; /* First neuron in monitored subset. */
                nmsim_elem_neuron_ix_t ineHi = trnek->ineHi; /* Last neuron in monitored subset. */
                assert((ineLo >= 0) && (ineLo <= ineHi) && (ineHi < nne));
                double X_tr, I_tr, J_tr;
                if (ineLo == ineHi)
                  { nmsim_elem_neuron_ix_t ine = ineLo;
                    X_tr = (double)X[ine]; I_tr = I[ine]; J_tr = J[ine];
                  }
                else
                  { X_tr = 0.0; I_tr = 0.0; J_tr = 0.0;
                    for (nmsim_elem_neuron_ix_t ine = ineLo; ine <= ineHi; ine++)
                      { X_tr += (double)X[ine]; I_tr += I[ine]; J_tr += J[ine]; }
                    double dn = (double)(ineHi - ineLo + 1);
                    X_tr /= dn; I_tr /= dn; J_tr /= dn;
                  }
                nmsim_elem_neuron_trace_set_X_I_J(trnek, t, X_tr, I_tr, J_tr);
              }
          }
      }
  }

void nmsim_elem_net_trace_write(char *prefix, nmsim_elem_net_trace_t *etrace)
  {
    for (nmsim_elem_neuron_ix_t ktr = 0; ktr < etrace->ntr; ktr++)
      { nmsim_elem_neuron_trace_t *trnek = etrace->trne[ktr];
        if (trnek != NULL)
          { nmsim_elem_neuron_ix_t ineLo = trnek->ineLo; /* First neuron in monitored subset. */
            nmsim_elem_neuron_ix_t ineHi = trnek->ineHi; /* Last neuron in monitored subset. */
            nmsim_time_t tLo = trnek->tLo;
            nmsim_time_t tHi = trnek->tHi;
            { /* Trace data by time step: */
              char *fname = NULL;
              asprintf(&fname, "%s_ne%010d--%010d_trace.txt", prefix, ineLo, ineHi);
              FILE *wr = open_write(fname, TRUE);
              nmsim_elem_neuron_trace_write(wr, trnek);
              fclose(wr);
              free(fname);
            }
            { /* Statistical summary: */
              char *fname = NULL;
              asprintf(&fname, "%s_ne%010d--%010d_stats.txt", prefix, ineLo, ineHi);
              FILE *wr = open_write(fname, TRUE);
              nmsim_elem_net_sim_stats_t *trS = nmsim_elem_net_sim_stats_new(ineLo, ineHi, tLo, tHi);
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
    nmsim_time_t tLo,               /* Discrete time of first neuron trace entries. */
    nmsim_time_t tHi,               /* Discrete time of last neuron trace entries. */
    nmsim_elem_neuron_count_t ntr   /* Number of neurons to monitor. */
  )
  { demand((ntr >= 0) & (ntr <= nne), "invalid monitored neuron count");
  
    /* Max number of neurons in each monitored set: */
    nmsim_elem_neuron_count_t nne_tr_MAX = 4;
    
    /* Create trace structure: */
    nmsim_elem_net_trace_t *etrace = nmsim_elem_net_trace_new(ntr);
    
    /* Select {ntr} random neurons out of {nne}, put in {perm[0..ntr-1]}: */
    uint64_t *perm = uint64_choose(nne, ntr, NULL);
    
    /* Create traces for those neurons with random time intervals: */
    for (nmsim_elem_neuron_ix_t ktr = 0; ktr < ntr; ktr++)
      { nmsim_elem_neuron_ix_t ineLok = (nmsim_elem_neuron_ix_t)perm[ktr];
        nmsim_elem_neuron_ix_t ineHik = (nmsim_elem_neuron_ix_t)(ineLok + int64_abrandom(1, nne_tr_MAX) - 1);
        if (ineHik >= nne) { ineHik = nne-1; }
        assert((ineLok >= 0) && (ineLok <= ineHik) && (ineHik < nne));
        nmsim_time_t tLok, tHik;
        nmsim_throw_time_range(tLo, tHi, &tLok, &tHik);
        nmsim_elem_neuron_trace_t *trnek = nmsim_elem_neuron_trace_new(ineLok, ineHik, tLok, tHik);
        nmsim_elem_net_trace_set(etrace, ktr, trnek);
      }
    return etrace;
  }
    
  
