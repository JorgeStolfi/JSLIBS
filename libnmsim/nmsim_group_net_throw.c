/* See {nmsim_group_net_throw.h} */
/* Last edited on 2020-12-11 19:55:41 by jstolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>

#include <nmsim_basic.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_synapse.h>
#include <nmsim_group_net.h>

#include <nmsim_group_net_throw.h>

nmsim_group_net_t *nmsim_group_net_throw
  ( nmsim_class_net_t *cnet,
    nmsim_group_neuron_count_t nng, 
    nmsim_group_synapse_count_t nsg,
    nmsim_elem_neuron_count_t nne, 
    nmsim_elem_synapse_count_t nse
  )
  {
    bool_t debug = FALSE;
    
    nmsim_class_neuron_count_t nnc = cnet->nnc;
    nmsim_class_synapse_count_t nsc = cnet->nsc;
    if (debug) { fprintf(stderr, "{nmsim_group_net_throw}: nnc = %d nsc = %d ...\n", nnc, nsc); }
    
    demand((nng >= 1) && (nng <= nmsim_group_neuron_count_MAX), "invalid neuron group count");
    demand((nsg >= 0) && (nsg <= nmsim_group_synapse_count_MAX), "invalid synapse group count");

    demand((nne >= 1) && (nne <= nmsim_elem_neuron_count_MAX), "invalid neuron count");
    demand((nse >= 0) && (nse <= nmsim_elem_synapse_count_MAX), "invalid synapse count");
    
    nmsim_group_net_t *gnet = nmsim_group_net_new(cnet, nng, nsg);

    /* Generate {nng} random neuron groups: */
    /* Leave {ngrp.{nne,ine_start,nsg_out}} at 0 for now. */
    if (debug) { fprintf(stderr, "generating %d neuron groups ...\n", nng); }
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_class_neuron_ix_t inc = (nmsim_class_neuron_ix_t)int64_abrandom(0, nnc-1);
        nmsim_group_neuron_t ngrp = (nmsim_group_neuron_t)
          { .inc = inc, .nne = 0, .ine_start = 0, .nsg_out = 0 };
        gnet->ngrp[ing] = ngrp;
      }

    /* Choose the number {ngrp.nsg_out} of synaptic bundles out */
    /* of each neuron group by throwing bundles at random: */
    if (debug) { fprintf(stderr, "choosing {.nsg_out} for each neuron group ...\n"); }
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++) 
      { gnet->ngrp[ing].nsg_out = 0; }
    for (nmsim_group_synapse_ix_t ksg = 0; ksg < nsg; ksg++)
      { nmsim_group_neuron_ix_t ing = (nmsim_group_neuron_ix_t)int64_abrandom(0, nng-1);
        gnet->ngrp[ing].nsg_out++;
      }
    
    /* Create those synaptic bundles, leaving {sgrp.{nse,ise_start}} at 0 for now. */
    if (debug) { fprintf(stderr, "creating the %d synaptic bundles ...\n", nsg); }
    nmsim_group_synapse_ix_t isg = 0; /* Index of next bundle. */
    for (nmsim_group_neuron_ix_t ing_pre = 0; ing_pre < nng; ing_pre++)
      { /* Get a neuron group {ngrp_pre}: */
        nmsim_group_neuron_t *ngrp_pre = &(gnet->ngrp[ing_pre]);
        /* Add its outgoing synapse groups: */
        nmsim_group_synapse_count_t nsg_out = ngrp_pre->nsg_out;
        for (nmsim_group_synapse_ix_t ksg = 0; ksg < nsg_out; ksg++)
          { nmsim_class_synapse_ix_t isc = (nmsim_class_synapse_ix_t)int64_abrandom(0, nsc-1);
            nmsim_group_neuron_ix_t ing_pos = (nmsim_class_synapse_ix_t)int64_abrandom(0, nng-1);
            nmsim_group_synapse_t sgrp = (nmsim_group_synapse_t)
              { .isc = isc, .ing_pre = ing_pre, .ing_pos = ing_pos, .nse = 0 };
            assert(isg < nsg);
            gnet->sgrp[isg] = sgrp;
            isg++;
          }
      }
    assert(isg == nsg);
 
    /* Choose the number of neurons {ngrp.nne} in each neuron group */
    /* (at least one on each group) by throwing neurons at random, */
    /* still leaving {ngrp.ine_start} at 0. */
    if (debug) { fprintf(stderr, "choosing the number of neurons {ngrp.nne} per group ...\n"); }
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++) 
      { gnet->ngrp[ing].nne = 1; }
    for (nmsim_elem_neuron_ix_t kne = 0; kne < nne - nng; kne++)
      { nmsim_group_neuron_ix_t ing = (nmsim_group_neuron_ix_t)int64_abrandom(0, nng-1);
        gnet->ngrp[ing].nne++;
      }

    /* Choose the number of synapses {sgrp.nse} in each synaptic bundle */
    /* by throwing synapses at random, still leaving {sgrp.ise_start} at 0. */
    if (debug) { fprintf(stderr, "choosing the number of synapses {sgrp.nse} per bundle ...\n"); }
    for (nmsim_group_synapse_ix_t isg = 0; isg < nsg; isg++) 
      { gnet->sgrp[isg].nse = 0; }
    for (nmsim_elem_synapse_ix_t kse = 0; kse < nse; kse++)
      { nmsim_group_synapse_ix_t isg = (nmsim_group_synapse_ix_t)int64_abrandom(0, nsg-1);
        gnet->sgrp[isg].nse++;
      }

    return gnet;
  }
