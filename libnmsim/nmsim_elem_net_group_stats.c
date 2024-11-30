/* See {nmsim_elem_net_group_stats.h} */
/* Last edited on 2021-01-06 13:42:33 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <sign.h>

#include <nmsim_group_neuron.h>
#include <nmsim_group_synapse.h>
#include <nmsim_group_net.h>

#include <nmsim_elem_neuron.h>
#include <nmsim_elem_synapse.h>
#include <nmsim_elem_net.h>

#include <nmsim_group_synapse_stats.h>

#include <nmsim_elem_net_group_stats.h>

nmsim_group_synapse_stats_t *nmsim_elem_net_group_stats_get
  ( nmsim_elem_net_t *enet, 
    nmsim_group_synapse_ix_t isg,
    sign_t sgn
  )
  {
    bool_t debug = FALSE;
    
    nmsim_elem_neuron_count_t nne = enet->nne; /* Total neurons in network. */
    nmsim_elem_synapse_count_t nse = enet->nse; /* Total synapses in network. */
    
    nmsim_group_net_t *gnet = enet->gnet;
    nmsim_group_neuron_count_t nng = gnet->nng; /* Total neuron groups in network. */
    nmsim_group_synapse_count_t nsg = gnet->nsg; /* Total synapse bundles in network. */
    
    assert((isg >= 0) && (isg  <nsg));
    
    nmsim_group_synapse_stats_t *st = notnull(malloc(sizeof(nmsim_group_synapse_stats_t)), "no mem");
    st->isg = isg;

    /* Get the synapse group {isg}: */
    nmsim_group_synapse_t *sgrp = &(gnet->sgrp[isg]);
    
    st->nse = sgrp->nse;

    /* Get the pre-synaptic neuron group {ing_pre} and its elem index ranges: */
    nmsim_group_neuron_ix_t ing_pre = sgrp->ing_pre;
    assert((ing_pre >= 0) && (ing_pre < nng));
    nmsim_group_neuron_t *ngrp_pre = &(gnet->ngrp[ing_pre]);
    nmsim_elem_neuron_count_t nne_pre = ngrp_pre->nne;
    assert((nne_pre >= 1) && (nne_pre <= nne)); /* Must have at least 1 neuron per group. */
    nmsim_elem_neuron_ix_t ine_pre_start = ngrp_pre->ine_start;
    nmsim_elem_neuron_ix_t ine_pre_lim = ine_pre_start + nne_pre;
    assert((ine_pre_start >= 0) && (ine_pre_start < ine_pre_lim) && (ine_pre_lim <= nne));
    st->ing_pre = ing_pre;
    st->nne_pre = nne_pre;

    /* Get the post-synaptic neuron group {ing_pos,neu_pos} and its elem index ranges: */
    nmsim_group_neuron_ix_t ing_pos = sgrp->ing_pos;
    assert((ing_pos >= 0) && (ing_pos < nng));
    nmsim_group_neuron_t *ngrp_pos = &(gnet->ngrp[ing_pos]);
    nmsim_elem_neuron_count_t nne_pos = ngrp_pos->nne;
    assert((nne_pos >= 1) && (nne_pos <= nne)); /* Must have at least 1 neuron per group. */
    nmsim_elem_neuron_ix_t ine_pos_start = ngrp_pos->ine_start;
    nmsim_elem_neuron_ix_t ine_pos_lim = ine_pos_start + nne_pos;
    assert((ine_pos_start >= 0) && (ine_pos_start < ine_pos_lim) && (ine_pos_lim <= nne));
    st->ing_pos = ing_pos;
    st->nne_pos = nne_pos;
    
    if (debug)
      { fprintf(stderr, "{nmsim_elem_net_group_stats_get}: ");
        fprintf(stderr, "isg = %d  ing_pre = %d ing_pos = %d", isg, ing_pre, ing_pos);
        fprintf(stderr, " ine_pre_start = %d ine_pre_lim = %d", ine_pre_start, ine_pre_lim);
        fprintf(stderr, " ine_pos_start = %d ine_pos_lim = %d", ine_pos_start, ine_pos_lim);
        fprintf(stderr, "\n");
      }

    /* Allocate tables {nse_in} with count of good synapses */
    /* and {Wtot_in} with total input synaptic weight, for each post-synaptic neuron, */ 
    nmsim_elem_synapse_count_t *nse_in = 
      notnull(malloc(nne_pos*sizeof(nmsim_elem_synapse_count_t)), "no mem");
    /* Allocate table : */ 
    double *Wtot_in = notnull(malloc(nne_pos*sizeof(double)), "no mem");
    /* Clear the tables for the post-synaptic neurons: */
    for (nmsim_elem_neuron_ix_t kne = 0;  kne < nne_pos; kne++)
      { nse_in[kne] = 0; 
        Wtot_in[kne] = 0.0;
      }

    /* Clear statistics record. . */
    st->nse = sgrp->nse;  /* Total number of synapses examined. */
    st->nse_cons = 0;     /* Number of synapses considered in the statistics. */
    // st->Wmax = -INF;      /* Weight with max absolute value. */
    // st->Wmin = +INF;      /* Weight with min absolute value. */

    double Wtot = 0.0;       /* Sum of weights of all considered synapses. */
    double Wtot2 = 0.0;      /* Sum of those weights squared. */

    nmsim_elem_synapse_count_t nse_proc = 0; /* Num synapses of this group examined so far. */
    nmsim_elem_synapse_count_t nse_zero = 0; /* Num synapses with zero weight. */
    /* Scan the neurons of the pre-synaptic group {ing_pre}: */
    for (nmsim_elem_neuron_ix_t ine_pre = ine_pre_start; ine_pre < ine_pre_lim; ine_pre++)
      { /* Get the neuron {ine_pre}: */
        nmsim_elem_neuron_t *neu_pre = &(enet->neu[ine_pre]);
        assert(neu_pre->ing == ing_pre);
        /* Scan the synapses out of that neuron, and process those with group {isg}: */
        nmsim_elem_synapse_ix_t ise_start = neu_pre->ise_out_start;
        nmsim_elem_synapse_ix_t ise_lim = ise_start + neu_pre->nse_out;
        assert((ise_start >= 0) && (ise_start <= ise_lim) && (ise_lim <= nse));
        for (nmsim_elem_synapse_ix_t ise = ise_start; ise < ise_lim; ise++)
          { /* Get the synapse {ise}: */
            nmsim_elem_synapse_t *syn = &(enet->syn[ise]);
            assert(syn->ine_pre == ine_pre);
            if (syn->isg == isg) 
              { /* Synapse belong to selected group. */
                if (syn->W == 0.0) 
                  { nse_zero++; }
                else if ((sgn == 0) || (((float)sgn)*syn->W > 0.0))
                  { /* Weight is nonzero and has the proper sign: */
                    st->nse_cons++;
                    // if (fabs(syn->W) > fabs(st->Wmax)) { st->Wmax = syn->W; }
                    // if (fabs(syn->W) < fabs(st->Wmin)) { st->Wmin = syn->W; }
                    Wtot += syn->W;
                    Wtot2 += syn->W*syn->W;
                    /* Accumulate weight and count of inputs of postsynaptic neuron: */
                    nmsim_elem_neuron_ix_t ine_pos = syn->ine_pos;
                    nmsim_elem_neuron_t *neu_pos = &(enet->neu[ine_pos]);
                    assert(neu_pos->ing == ing_pos);
                    if ((ine_pos < ine_pos_start) || (ine_pos >= ine_pos_lim))
                      { fprintf(stderr, "** ");
                        fprintf(stderr, "isg = %d  ing_pre = %d ing_pos = %d", isg, ing_pre, ing_pos);
                        fprintf(stderr, " ise = %d  ine_pre = %d ine_pos = %d", ise, ine_pre, ine_pos);
                        fprintf(stderr, " ine_pos_start = %d ine_pos_lim = %d", ine_pos_start, ine_pos_lim);
                        fprintf(stderr, "\n");
                      }
                    assert((ine_pos >= ine_pos_start) && (ine_pos < ine_pos_lim));
                    nmsim_elem_neuron_ix_t kne_pos = ine_pos - ine_pos_start;
                    Wtot_in[kne_pos] += syn->W;
                    nse_in[kne_pos]++;
                  }
                nse_proc++;
              }
          }
      }
    assert(nse_proc == sgrp->nse);
    if (sgn == 0) { assert(st->nse_cons + nse_zero == sgrp->nse); }
    
    /* Convert the sums {Wtot,Wtot2} into average and deviation: */
    nmsim_elem_synapse_count_t n = st->nse_cons; /* Synapses consdered in stats. */
    double dn = (double)n;
    st->Wavg = (n < 1 ? 0.0 : Wtot/dn);
    st->Wdev = (n < 2 ? 0.0: sqrt(fmax(0.0, Wtot2 - dn*st->Wavg*st->Wavg)/(dn-1)));
    
    /* Compute statistics of total input synaptic weight in post-synaptic neuron group: */
    double Wtot_tot = 0.0;  /* Sum of {Wtot_in[i]} */
    double Wtot_tot2 = 0.0; /* Sum of {Wtot_in[i]^2} */
    nmsim_elem_synapse_count_t nin_tot = 0; /* Total inputs in neuron group seen. */
    for (nmsim_elem_neuron_ix_t kne = 0; kne < nne_pos; kne++)
      { nmsim_elem_synapse_count_t nin_ine = nse_in[kne];
        double Wtot_ine = Wtot_in[kne];
        Wtot_tot += Wtot_ine;
        Wtot_tot2 += Wtot_ine*Wtot_ine;
        nin_tot += nin_ine;
      }
    assert(nin_tot == st->nse_cons);
    nmsim_elem_synapse_count_t m = nne_pos; /* Neurons in post-synaptic group. */
    double dm = (double)m;
    st->navg = dn/dm;
    st->Wina = Wtot_tot/dm;
    st->Wind = (m < 2 ? 0.0: sqrt(fmax(0.0, Wtot_tot2 - dm*st->Wina*st->Wina)/(dm-1)));

    free(Wtot_in);
    free(nse_in);
    
    return st;
  }
    
void nmsim_elem_net_group_stats_print_one
  ( FILE *wr, 
    nmsim_elem_net_t *enet,
    nmsim_group_synapse_ix_t isg,
    char *pref,
    int32_t hdr,
    char *suff
  )
  {
    if (pref != NULL) { fputs(pref, wr); }
    if (hdr == 0)
      { assert((isg >= 0) && (isg < enet->gnet->nsg));
        
        /* Gather statistics of all synapses, excitatory only, and inhibitory only: */
        nmsim_group_synapse_stats_t *st_all = nmsim_elem_net_group_stats_get(enet, isg, 00);
        nmsim_group_synapse_stats_t *st_exc = nmsim_elem_net_group_stats_get(enet, isg, +1);
        nmsim_group_synapse_stats_t *st_inh = nmsim_elem_net_group_stats_get(enet, isg, -1);
        
        /* Print relevant data for each subset: */
        nmsim_group_synapse_stats_print_gen(wr, st_all, 0); fputs(" ", wr);
        nmsim_group_synapse_stats_print_uns(wr, st_all, 0); fputs(" ", wr);
        nmsim_group_synapse_stats_print_sgn(wr, st_exc, +1, 0); fputs(" ", wr);
        nmsim_group_synapse_stats_print_sgn(wr, st_inh, -1, 0);
        
        free(st_all);
        free(st_exc);
        free(st_inh);
        
      }
    else
      { 
        nmsim_group_synapse_stats_print_gen(wr, NULL, hdr); fputs(" ", wr);
        nmsim_group_synapse_stats_print_uns(wr, NULL, hdr); fputs(" ", wr);
        nmsim_group_synapse_stats_print_sgn(wr, NULL, +1, hdr); fputs(" ", wr);
        nmsim_group_synapse_stats_print_sgn(wr, NULL, -1, hdr);
      }
    if (suff != NULL) { fputs(suff, wr); }
  }
  
void nmsim_elem_net_group_stats_print_all
  ( FILE *wr, 
    nmsim_elem_net_t *enet,
    char *pref,
    char *suff
  )
  {
    fprintf(wr, "The network has %d neurons in %d groups\n", enet->nne, enet->gnet->nng);
    fprintf(wr, "The network has %d synapses in %d groups\n", enet->nse, enet->gnet->nsg);
    fprintf(wr, "\n\n");

    /* Print table header: */
    for (uint32_t hdr = 1;  hdr <= 4; hdr++)
      { nmsim_elem_net_group_stats_print_one(wr, NULL, -1, pref, hdr, suff); }

    /* Print table body: */
    nmsim_group_synapse_count_t nsg = enet->gnet->nsg;
    for (nmsim_group_synapse_ix_t isg = 0; isg < nsg; isg++)
      { nmsim_group_synapse_t *sgrp = &(enet->gnet->sgrp[isg]);
        if (sgrp->nse > 0)
          { nmsim_elem_net_group_stats_print_one(wr, enet, isg, pref, 0, suff); }
      }
      
    /* Print line of dashes: */
    nmsim_elem_net_group_stats_print_one(wr, NULL, -1, pref, 4, suff);
    
    /* Print table legend: */
    fputs("\n\n", wr);
    fputs(nmsim_group_synapse_stats_print_INFO, wr);
    fputs("\n", wr);
    fflush(wr);
  }
