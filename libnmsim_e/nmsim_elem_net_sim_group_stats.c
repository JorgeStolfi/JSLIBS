/* See {nmsim_elem_net_sim_group_stats.h} */
/* Last edited on 2020-12-17 01:17:11 by jstolfi */

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
#include <jsfile.h>

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

#include <nmsim_elem_net_sim_group_stats.h>

nmsim_elem_net_sim_group_stats_t *nmsim_elem_net_sim_group_stats_new
  ( nmsim_group_net_t *gnet,   /* Data on neuron groups, including neuron index ranges. */
    nmsim_time_t tLo,         /* Discrete time of first state to consider. */
    nmsim_time_t tHi          /* Discrete time of last state to consider. */
  )
  {
    if (tLo > tHi) { return NULL; }
    nmsim_group_neuron_count_t nng = gnet->nng;
    nmsim_elem_net_sim_group_stats_t *gstats = notnull(malloc(sizeof(nmsim_elem_net_sim_group_stats_t)), "no mem");
    gstats->nng = nng;
    gstats->tLo = tLo;
    gstats->tHi = tHi;
    gstats->S = notnull(malloc(nng*sizeof(nmsim_elem_net_sim_stats_t *)), "no mem");
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_elem_neuron_ix_t ineLo = gnet->ngrp[ing].ine_start;
        nmsim_elem_neuron_ix_t ineHi = ineLo + gnet->ngrp[ing].nne - 1;
        nmsim_elem_net_sim_stats_t *S_g = nmsim_elem_net_sim_stats_new(ineLo, ineHi, tLo, tHi);
        gstats->S[ing] = S_g;
      }
    return gstats;
  }

void nmsim_elem_net_sim_group_stats_write(char *outPrefix, nmsim_elem_net_sim_group_stats_t *gstats)
  {
    if (gstats == NULL) { return; }
    nmsim_group_neuron_count_t nng = gstats->nng;
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_elem_net_sim_stats_t *S_g = gstats->S[ing];
        char *fname = jsprintf("%s_ng%010d_stats.txt", outPrefix, ing);
        FILE *wr = open_write(fname, TRUE);
        nmsim_elem_net_sim_stats_write(wr, S_g);
        fclose(wr);
        free(fname);
      }
  }

void nmsim_elem_net_sim_group_stats_initialize(nmsim_elem_net_sim_group_stats_t *gstats)
  { 
    if (gstats == NULL) { return; }
    nmsim_group_neuron_count_t nng = gstats->nng;
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_elem_net_sim_stats_t *S_g = gstats->S[ing];
        nmsim_elem_net_sim_stats_initialize(S_g);
      }
  }
   
void nmsim_elem_net_sim_group_stats_accumulate_V_age_M_H
  ( nmsim_elem_net_sim_group_stats_t *gstats,
    nmsim_time_t t,
    nmsim_elem_neuron_count_t nne,
    double V[],                /* Neuron potentials at time {t}. */
    nmsim_step_count_t age[],  /* Ages of neurons at time {t}. */
    double M[],                /* Recharge modulators at time {t}. */
    double H[]                 /* Output modulators at time {t}. */
  )
  {
    if (gstats == NULL) { return; }
    nmsim_group_neuron_count_t nng = gstats->nng;
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_elem_net_sim_stats_t *S_g = gstats->S[ing];
        nmsim_elem_net_sim_stats_accumulate_V_age_M_H(S_g, t, nne, V, age, M, H);
      }
  }
    
void nmsim_elem_net_sim_group_stats_accumulate_VF_AF_X_I_J
  ( nmsim_elem_net_sim_group_stats_t *gstats,
    nmsim_time_t t,
    nmsim_elem_neuron_count_t nne,
    double V[],                /* Neuron potentials at time {t}. */
    nmsim_step_count_t age[],  /* Ages of neurons at time {t}. */
    bool_t X[],                /* Firing indicators for step {t} to {t+1}. */
    double I[],                /* External neuron inputs for step {t} to {t+1}. */
    double J[]                 /* Total inputs for step {t} to {t+1}. */
  )
  {
    if (gstats == NULL) { return; }
    nmsim_group_neuron_count_t nng = gstats->nng;
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_elem_net_sim_stats_t *S_g = gstats->S[ing];
        nmsim_elem_net_sim_stats_accumulate_VF_AF_X_I_J(S_g, t, nne, V, age, X, I, J);
      }
  }

void nmsim_elem_net_sim_group_stats_finalize(nmsim_elem_net_sim_group_stats_t *gstats)
  { 
    if (gstats == NULL) { return; }
    nmsim_group_neuron_count_t nng = gstats->nng;
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_elem_net_sim_stats_t *S_g = gstats->S[ing];
        nmsim_elem_net_sim_stats_finalize(S_g);
      }
  }
