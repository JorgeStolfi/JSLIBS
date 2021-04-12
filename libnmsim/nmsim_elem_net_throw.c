/* See {nmsim_elem_net_throw.h} */
/* Last edited on 2020-12-16 01:32:16 by jstolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsrandom.h>

#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_synapse.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_synapse.h>
#include <nmsim_elem_net.h>

#include <nmsim_elem_net_throw.h>

void nmsim_elem_net_throw_elems(nmsim_elem_net_t *enet);
  /* Adds to the network {enet} the neuron and synapse elements
    according to group parameters from {enet->gnet}.
    
    Namely, for each neuron group {ngrp} in the network {enet}, adds
    {ngrp.nne} neurons, assigned to that group. 
    
    Then, for each synapse group {sgrp} in {enet},
    adds synapses from neurons of group {A = sgrp.ign_pre}
    to neurons of group {B = sgrp.ing_pos}. Specifically, the procedure will add on
    average {sgrp.K} incoming synapses to each neuron in {B}. The
    strengths of those synapses will be drawn from the distrution
    specified in the synapse class of that bundle, {enet->pnet->cnet->sclass[sgrp.isc]}.
    
    The presynaptic neuron or each synapse is chosen among the {A} neurons
    at random, with replacement.  Therefore, there may be more than one
    synapse from that bundle between a given pair of neurons. 
    
    The strength of each synapse will be chosen at random from a Gaussian
    distribution, whose mean {W_avg} and deviation {W_dev} are defined
    in the synapse class of the bundle.  */

void nmsim_elem_net_expand_neuron_group
  ( nmsim_elem_net_t *enet,
    nmsim_group_neuron_ix_t ing,
    nmsim_elem_neuron_ix_t *ineP
  );
  /* Stores into the vector {enet->neu} a set of {nne_g}
    neurons, all from the neuron group {ing}.   
    
    The neurons will be in conseuctive positions starting at index {*ineP}.
    The field {enet->gnet->bgrp[ing].ine_start} is set to {*ineP}.
    Increments {*ineP} with the number of neurons created.*/

void nmsim_elem_net_expand_synapse_group
  ( nmsim_elem_net_t *enet,
    nmsim_group_synapse_ix_t isg,
    nmsim_elem_neuron_ix_t ine0,
    nmsim_elem_neuron_ix_t ine1,
    nmsim_elem_neuron_ix_t jne0,
    nmsim_elem_neuron_ix_t jne1
  );
  /* Appends to the vector {enet->syn[0..nse-1]} a number of random synapses
    from synapse group {isg}, from neurons {ine0..ine1} to neurons {jne0..jne1}.
    
    Each neuron {jne} in {jne0..jne1} will have {K} input synapses, where
    {K} is {enet->gnet->sgrp[isg].K}.  The source neuron {ine} of each synapse will be
    selected randomly from the neurons {ine0..ine1}, with replacement
    (that is, the procedure my add two or more synapses coming from the
    same neuron {ine}. 
    
    The strength of each synapse will be chosen at random from a Gaussian
    distribution, whose mean {W_avg} and deviation {W_dev} are defined
    in the synapse class record {enet->gnet->cnet->sclass[isc]}, where
    {isc} is {enet->gnet->sgrp[isg].isc}. 
    
    The vector {enet->syn} will be extended by the function if and when needed. */

/* IMPLEMENTATIONS */

nmsim_elem_net_t *nmsim_elem_net_throw(nmsim_group_net_t *gnet)
  { 
    bool_t debug = FALSE;
    
    nmsim_group_neuron_count_t nng = gnet->nng;
    nmsim_group_synapse_count_t nsg = gnet->nsg;
    
    if (debug) { fprintf(stderr, "begin {nmsim_elem_net_throw}\n"); }
    if (debug) { fprintf(stderr, "groups: neurons = %d synapses = %d\n", nng, nsg); }
    
    /* Compute the total number of neurons {nne}: */
    nmsim_elem_neuron_count_t nne = 0;
    for (nmsim_group_neuron_ix_t ing = 0; ing < gnet->nng; ing++) { nne +=  gnet->ngrp[ing].nne; }
      
    /* Compute the total number of synapses {nse}: */
    nmsim_elem_synapse_count_t nse = 0;
    for (nmsim_group_synapse_ix_t isg = 0; isg < gnet->nsg; isg++) { nse +=  gnet->sgrp[isg].nse; }
    
    if (debug) { fprintf(stderr, "elements: neurons = %d synapses = %d\n", nne, nse); }

    /* Allocate the element-level network: */
    nmsim_elem_net_t *enet = nmsim_elem_net_new(gnet, nne, nse);
    
    /* Create the neuron elements for each neuron group, */
    /* and define the fields {neu.ine_start}: */
    nmsim_elem_neuron_ix_t ine = 0; /* Index of next neuron in table. */
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_group_neuron_t *ngrp = &(gnet->ngrp[ing]);
        nmsim_elem_neuron_count_t nne_g = ngrp->nne;
        ngrp->ine_start = ine;
        if (debug) { fprintf(stderr, "generating neurons of group %d (%d..%d)\n", ing, ine, ine+nne_g-1); }
        /* Generate the neurons of this group, leaving */
        /* {neu.nse_out} and {neu.ise_ut_start} at zero: */
        for (nmsim_elem_neuron_ix_t kne = 0; kne < nne_g; kne++)
          { nmsim_elem_neuron_t neu = (nmsim_elem_neuron_t)
              { .ing = ing, .nse_out = 0, .ise_out_start = 0 };
            assert(ine < nne);
            enet->neu[ine] = neu;
            ine++;
          }
      }
    assert(ine == nne);
    
    /* Generate the synapse elements in any order, temporarily, */
    /* by throwing {sgrp.nne} synapses at random in each synapse group. */
    /* Define the number {neu.nse_out} of synapses out of each neuron. */
    
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++)
      { enet->neu[ine].nse_out = 0; }
    nmsim_elem_synapse_ix_t ise = 0; /* Index of next free slot in {enet->syn}. */
    for (nmsim_group_synapse_ix_t isg = 0; isg < nsg; isg++) 
      { nmsim_group_synapse_t *sgrp = &(gnet->sgrp[isg]);
        nmsim_elem_synapse_count_t nse_g = sgrp->nse; /* Synapses in group. */
            
        /* Get the pre-synaptic neuron group of this synapse group: */
        nmsim_group_neuron_ix_t ing_pre = sgrp->ing_pre;
        assert((ing_pre >= 0) && (ing_pre < nng));
        nmsim_group_neuron_t *ngrp_pre = &(gnet->ngrp[ing_pre]);
        nmsim_elem_neuron_ix_t ine_pre_start = ngrp_pre->ine_start;
        nmsim_elem_neuron_ix_t ine_pre_lim = ine_pre_start + ngrp_pre->nne;

        /* Get the post-synaptic neuron group of this synapse group: */
        nmsim_group_neuron_ix_t ing_pos = sgrp->ing_pos;
        assert((ing_pos >= 0) && (ing_pos < nng));
        nmsim_group_neuron_t *ngrp_pos = &(gnet->ngrp[ing_pos]);
        nmsim_elem_neuron_ix_t ine_pos_start = ngrp_pos->ine_start;
        nmsim_elem_neuron_ix_t ine_pos_lim = ine_pos_start + ngrp_pos->nne;

        if (debug) 
          { fprintf(stderr, "throwing %d sinapses in bundle %d ", nse_g, isg);
            fprintf(stderr, " from neuron group %d (%d..%d)\n", ing_pre, ine_pre_start, ine_pre_lim-1);
            fprintf(stderr, " to neuron group %d (%d..%d)\n", ing_pos, ine_pos_start, ine_pos_lim-1);
          }
        
        /* Throw {nse_g} synapses between those neurons: */
        for (nmsim_elem_synapse_ix_t kse = 0; kse < nse_g; kse++)
          { 
            /* Choose the pre synaptic neuron index: */
            nmsim_elem_neuron_ix_t ine_pre = (nmsim_elem_neuron_ix_t)
              int64_abrandom(ine_pre_start, ine_pre_lim-1);
            assert((ine_pre >= 0) && (ine_pre < nne));
            nmsim_elem_neuron_t *neu_pre = &(enet->neu[ine_pre]);
            
            /* Choose the post-synaptic neuron index: */
            nmsim_elem_neuron_ix_t ine_pos = (nmsim_elem_neuron_ix_t)
              int64_abrandom(ine_pos_start, ine_pos_lim-1);
            assert((ine_pos >= 0) && (ine_pos < nne));
              
            if (debug) { fprintf(stderr, "threw 1 sinapse from %d to %d\n", ine_pre, ine_pos); }

            /* Choose the synapse weight: */
            nmsim_class_synapse_ix_t isc = sgrp->isc;
            nmsim_class_synapse_t *sclass = gnet->cnet->sclass[isc];
            float W = (float)dloggaussrand(sclass->W_avg, sclass->W_dev);
            
            /* Pack and store: */
            assert(ise < nse);
            enet->syn[ise] = (nmsim_elem_synapse_t)
              { .isg = isg, .ine_pre = ine_pre, .ine_pos = ine_pos, .W = W };
            ise++;
              
            /* Count it as one more synapse out of neuron {ine_pre}: */
            neu_pre->nse_out++;
          }
       }
    assert(ise == nse);
    
    /* Define the fields {neu.ise_out_start} and reset temporarily */
    /* the fields {neu.nse_out} to zero: */
    if (debug) { fprintf(stderr, "computing fields {neu.ise_out_start} ...\n"); }
    nmsim_elem_synapse_count_t nse_out_tot = 0; /* Sum of {neu.nse_out}. */
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++)
      { nmsim_elem_neuron_t *neu = &(enet->neu[ine]);
        if (debug) { fprintf(stderr, "threw total %d sinapses out of neuron %d\n", neu->nse_out, ine); }
        neu->ise_out_start = nse_out_tot;
        nse_out_tot += neu->nse_out;
        neu->nse_out = 0;
      }
    assert(nse_out_tot == nse);   

    /* Sort the synapses using the fields {neu.ise_out_start,neu.nse_out}: */
    if (debug) { fprintf(stderr, "sorting synapses ...\n"); }
    nmsim_elem_synapse_t *syn_sorted = notnull(malloc(nse*sizeof(nmsim_elem_synapse_t)),"no mem");
    for (nmsim_elem_synapse_ix_t ise = 0; ise < nse; ise++)
      { nmsim_elem_synapse_t *syn = &(enet->syn[ise]);
        nmsim_elem_neuron_ix_t ine_pre = syn->ine_pre;
        assert((ine_pre >= 0) && (ine_pre < nne));
        nmsim_elem_neuron_t *neu_pre = &(enet->neu[ine_pre]);
        nmsim_elem_synapse_ix_t jse = neu_pre->ise_out_start + neu_pre->nse_out;
        assert((jse >= 0) && (jse < nse));
        syn_sorted[jse] = (*syn);
        neu_pre->nse_out++;
      }
    free(enet->syn);
    enet->syn = syn_sorted;
        
    /* Paranoia: check consistency of {neu.nse_out,neu.ise_out_pre}: */
    ise = 0; /* Index of next undefined synapse. */
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++)
      { nmsim_elem_neuron_t *neu = &(enet->neu[ine]);
        assert(neu->ise_out_start == ise);
        ise += neu->nse_out;
      }
    assert(ise == nse); 
    
    return enet;
  }
