/* See {nmsim_elem_net_throw.h} */
/* Last edited on 2020-12-11 19:55:29 by jstolfi */

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
    
    /* Compute the total number of neurons {nne}: */
    nmsim_elem_neuron_count_t nne = 0;
    for (nmsim_group_neuron_ix_t ing = 0; ing < gnet->nng; ing++) { nne +=  gnet->ngrp[ing].nne; }
      
    /* Compute the total number of synapses {nse}: */
    nmsim_elem_synapse_count_t nse = 0;
    for (nmsim_group_synapse_ix_t isg = 0; isg < gnet->nsg; isg++) { nse +=  gnet->sgrp[isg].nse; }
    
    /* Allocate the element-level network: */
    nmsim_elem_net_t *enet = nmsim_elem_net_new(gnet, nne, nse);
    
    /* Create the neuron elements for each neuron group, */
    /* and define the fields {neu.ine_start}: */
    nmsim_elem_neuron_ix_t ine = 0; /* Index of next neuron in table. */
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_group_neuron_t *ngrp = &(gnet->ngrp[ing]);
        nmsim_elem_neuron_count_t nne_g = ngrp->nne;
        ngrp->ine_start = ine;
        /* Generate the neurons of this group, leaving */
        /* {neu.nse_out} and {neu.ise_ut_start} at zero: */
        for (nmsim_elem_neuron_ix_t kne = 0; kne< nne_g; kne++)
          { nmsim_elem_neuron_t neu = (nmsim_elem_neuron_t)
              { .ing = ing, .nse_out = 0, .ise_out_start = 0 };
            enet->neu[ine] = neu;
            ine++;
          }
      }
    assert(ine == nne);
    
    /* Choose the number {neu.nse_out} of synapses out of each neuron, */
    /* by throwing {sgrp.nne} synapses at random in each synapse group: */
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++)
      { enet->neu[ine].nse_out = 0; }
    for (nmsim_group_synapse_ix_t isg = 0; isg < nsg; isg++) 
      { nmsim_group_synapse_t *sgrp = &(gnet->sgrp[isg]);
        nmsim_elem_synapse_count_t nse_g = sgrp->nse; /* Synapses in group. */
        if (debug) { fprintf(stderr, "throwing %d sinapses in bundle %d\n", nse_g, isg); }
            
        /* Get the pre-synaptic neuron group of this synapse group: */
        nmsim_group_neuron_ix_t ing_pre = sgrp->ing_pre;
        nmsim_group_neuron_t *ngrp_pre = &(gnet->ngrp[ing_pre]);
        
        for (nmsim_elem_synapse_ix_t kse = 0; kse < nse_g; kse++)
          { /* Pick a random neuron in neuron group {ing_pre}: */
            nmsim_elem_neuron_ix_t ine_pre_start = ngrp_pre->ine_start;
            nmsim_elem_neuron_ix_t ine_pre_lim = ine_pre_start + ngrp_pre->nne;
            nmsim_elem_neuron_ix_t ine_pre = 
              (nmsim_elem_neuron_ix_t)int64_abrandom(ine_pre_start, ine_pre_lim-1);
            if (debug) { fprintf(stderr, "threw 1 sinapse out of %d\n", ine_pre); }
            /* Count it as one more synapse out of neuron {ine_pre}: */
            enet->neu[ine_pre].nse_out++;
          }
       }
        
    /* Define the fields {neu.ise_out_start} and reset temporarily */
    /* the fields {neu.nse_out} to zero: */
    nmsim_elem_synapse_ix_t ise = 0; /* Index of next undefined synapse. */
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++)
      { nmsim_elem_neuron_t *neu = &(enet->neu[ine]);
        if (debug) { fprintf(stderr, "threw %d sinapses out of neuron %d\n", neu->nse_out, ine); }
        neu->ise_out_start = ise;
        ise += neu->nse_out;
        neu->nse_out = 0;
      }
    assert(ise == nse);   

    /* Create the synapse elements for each neuron, */
    /* and restore the fields {neu.nse_out}: */
    for (nmsim_group_synapse_ix_t isg = 0; isg < nsg; isg++) 
      { nmsim_group_synapse_t *sgrp = &(gnet->sgrp[isg]);
        nmsim_elem_synapse_count_t nse_g = sgrp->nse; /* Synapses in group. */
        if (debug) { fprintf(stderr, "generating %d sinapses in bundle %d\n", nse_g, isg); }
        
        /* Get the class of this synaptic group: */
        nmsim_class_synapse_ix_t isc = sgrp->isc;
        nmsim_class_synapse_t *sclass = gnet->cnet->sclass[isc];
        
        /* Get the pre- and post-synaptic neuron groups of this synapse group: */
        nmsim_group_neuron_ix_t ing_pre = sgrp->ing_pre;
        nmsim_group_neuron_t *ngrp_pre = &(gnet->ngrp[ing_pre]);
        nmsim_elem_neuron_ix_t ine_pre_start = ngrp_pre->ine_start;
        nmsim_elem_neuron_ix_t ine_pre_lim = ine_pre_start + ngrp_pre->nne;

        nmsim_group_neuron_ix_t ing_pos = sgrp->ing_pos;
        nmsim_group_neuron_t *ngrp_pos = &(gnet->ngrp[ing_pos]);
        nmsim_elem_neuron_ix_t ine_pos_start = ngrp_pos->ine_start;
        nmsim_elem_neuron_ix_t ine_pos_lim = ine_pos_start + ngrp_pos->nne;
        
        /* Generate the synapse elements in this group and insert them */
        /* in the proper places, incrementinng the fields {.nse_out}: */
        for (nmsim_elem_synapse_ix_t kse = 0; kse < nse_g; kse++)
          { /* Choose the pre and post synaptic indices: */
            nmsim_elem_neuron_ix_t ine_pre = (nmsim_elem_neuron_ix_t)
              int64_abrandom(ine_pre_start, ine_pre_lim-1);
            nmsim_elem_neuron_t *neu_pre = &(enet->neu[ine_pre]);
            
            nmsim_elem_neuron_ix_t ine_pos = (nmsim_elem_neuron_ix_t)
              int64_abrandom(ine_pos_start, ine_pos_lim-1);
              
            /* Decide where the synapse will go: */
            nmsim_elem_synapse_ix_t ise = neu_pre->ise_out_start + neu_pre->nse_out;
            
            /* Choose the synapse weight: */
            float W = (float)dloggaussrand(sclass->W_avg, sclass->W_dev);
            
            if (debug) { fprintf(stderr, "generating sinapse %d from %d %d\n", ise, ine_pre, ine_pos); }

            /* Pack and store: */
            enet->syn[ise] = (nmsim_elem_synapse_t)
              { .isg = isg, .ine_pre = ine_pre, .ine_pos = ine_pos, .W = W };
              
            /* Count one more for neuron {ine_pre}: */
            neu_pre->nse_out++;
         }
      }
    
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
