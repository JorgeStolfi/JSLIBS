/* See {nmsim_elem_net.h} */
/* Last edited on 2020-12-07 15:50:45 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsrandom.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>

#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_synapse.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_synapse.h>

#include <nmsim_elem_net.h>

/* LOW-LEVEL EXPANSION PROCEDURES */

nmsim_elem_net_t *nmsim_elem_net_new(nmsim_group_net_t *gnet);
  /* Allocates an {nmsim_elem_net_t} structure with the group-level
    network {gnet}, initially with zero neurons and zero synapses.
    The description is invalid, until the synapse and
    neuron groups have been expanded -- e.g. with
    {nmsim_elem_net_expand_groups}. */
    
void nmsim_elem_net_expand_groups(nmsim_elem_net_t *enet);
  /* Expands the neuron and synapse groups of {enet->gnet}
    into neuron and synapse elements
    
    Namely, for each neuron group {ngrp} in the network {enet}, adds
    {ngrp.nne} neurons, assigned to that group. 
    
    Then, for each synapse group {sgrp} in {enet},
    adds synapses from neurons of group {A = sgrp.ign_pre}
    to neurons of group {B = sgrp.ing_pos}. Specifically, the procedure will add on
    average {sgrp.K} incoming synapses to each neuron in {B}. The
    strengths of those synapses will be drawn from the distrution
    specified in the synapse class of that bundle, {enet->pnet->cnet->sclass.e[sgrp.isc]}.
    
    The presynaptic neuron or each synapse is chosen among the {A} neurons
    at random, with replacement.  Therefore, there may be more than one
    synapse from that bundle between a given pair of neurons. 
    
    The strength of each synapse will be chosen at random from a Gaussian
    distribution, whose mean {W_avg} and deviation {W_dev} are defined
    in the synapse class of the bundle.  */

void nmsim_elem_net_expand_neuron_group
  ( nmsim_elem_net_t *enet,
    nmsim_group_neuron_ix_t ing,
    nmsim_elem_neuron_count_t nne_g
  );
  /* Appends to the vector {enet->neu} a list of {nne_g}
    neurons, all from the neuron group {ing}. */

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
    {K} is {enet->gnet->sgrp.e[isg].K}.  The source neuron {ine} of each synapse will be
    selected randomly from the neurons {ine0..ine1}, with replacement
    (that is, the procedure my add two or more synapses coming from the
    same neuron {ine}. 
    
    The strength of each synapse will be chosen at random from a Gaussian
    distribution, whose mean {W_avg} and deviation {W_dev} are defined
    in the synapse class record {enet->gnet->cnet->sclass.e[isc]}, where
    {isc} is {enet->gnet->sgrp.e[isg].isc}. 
    
    The vector {enet->syn} will be extended by the function if and when needed. */

nmsim_elem_neuron_ix_t nmsim_elem_net_add_neuron
  ( nmsim_elem_net_t *enet, 
    nmsim_elem_neuron_t *neu
  );
  /* Adds to the network {enet} a new neuron with attributes
    {*neu}. The neuron is stored in the {enet->neu} vector, which
    is expanded if necessary, and incrementng {enet->nne}. Returns the
    neuron index {ine}. */
    
nmsim_elem_synapse_ix_t nmsim_elem_net_add_synapse
  ( nmsim_elem_net_t *enet,
    nmsim_elem_synapse_t *syn
  );
  /* Adds to the network {enet} a new synapse with attributes
    {*syn}. The synapse is stored in the {enet->syn} vector, which is
    expanded if necessary, and incrementng {enet->nse}. Returns the
    synapse index {ise}. */

nmsim_elem_net_t *nmsim_elem_net_new(nmsim_group_net_t *gnet)
  { 
    nmsim_elem_net_t *enet = notnull(malloc(sizeof(nmsim_elem_net_t)), "no mem");
    (*enet) = (nmsim_elem_net_t)
      { .gnet = gnet,
        .nne = 0, .nse = 0,
        .neu = nmsim_elem_neuron_vec_new(20),
        .syn = nmsim_elem_synapse_vec_new(20)
      };
    return enet;
  }

nmsim_elem_net_t *nmsim_elem_net_throw(nmsim_group_net_t *gnet)
  { 
    nmsim_elem_net_t *enet = nmsim_elem_net_new(gnet);
    nmsim_elem_net_expand_groups(enet);
    return enet;
  }
  
void nmsim_elem_net_free(nmsim_elem_net_t *enet)    
  {
    if (enet != NULL)
      { if (enet->gnet != NULL) { nmsim_group_net_free(enet->gnet); }
        if (enet->neu.e != NULL) { free(enet->neu.e); }
        if (enet->syn.e != NULL) { free(enet->syn.e); }
        free(enet);
      }
  }
     
void nmsim_elem_net_expand_groups(nmsim_elem_net_t *enet)
  {
    nmsim_group_net_t *gnet = enet->gnet; /* Neuron and synapse groups. */
    nmsim_group_neuron_count_t nng = gnet->nng;
    nmsim_group_synapse_count_t nsg = gnet->nsg;
    
    /* The neurons of a group {ing} in {0..nng-1} will be */
    /* {ine_start[ing]} to  {ine_start[ing+1]-1}: */
    nmsim_elem_neuron_ix_t *ine_start = notnull(malloc((nng+1)*sizeof(nmsim_elem_neuron_ix_t)), "not null");
    ine_start[0] = enet->nne; /* Start of group 0 neurons. */
    
    /* Expand the neuron groups: */
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_elem_neuron_count_t nne_g = gnet->ngrp.e[ing].nne; /* Num of neurons in group. */ 
        nmsim_elem_net_expand_neuron_group(enet, ing, nne_g);
        ine_start[ing+1] = enet->nne; /* First neuron not in group. */
      }

    /* Expand synapse groups: */
    for (nmsim_group_synapse_ix_t isg = 0; isg < nsg; isg++)
      { /* Get the indices of pre- and post-synaptic neuron groups {ing_pre,ing_pos}: */
        nmsim_group_neuron_ix_t ing_pre = gnet->sgrp.e[isg].ing_pre;
        nmsim_group_neuron_ix_t ing_pos = gnet->sgrp.e[isg].ing_pos;
        assert((ing_pre >= 0) && (ing_pre < nng));
        assert((ing_pos >= 0) && (ing_pos < nng));

        /* Get the index ranges of pre- and post-synaptic neurons: */
        nmsim_elem_neuron_ix_t ine_pre0 = ine_start[ing_pre];
        nmsim_elem_neuron_ix_t ine_pre1 = ine_start[ing_pre+1] - 1;
        nmsim_elem_neuron_ix_t ine_pos0 = ine_start[ing_pos];
        nmsim_elem_neuron_ix_t ine_pos1 = ine_start[ing_pos+1] - 1;
        
        /* Now expand the synapse group: */
        nmsim_elem_net_expand_synapse_group(enet, isg, ine_pre0, ine_pre1, ine_pos0, ine_pos1); 
      }
  }

nmsim_elem_neuron_ix_t nmsim_elem_net_add_neuron
  ( nmsim_elem_net_t *enet, 
    nmsim_elem_neuron_t *neu
  )
  { demand((neu->ing >= 0) && (neu->ing < enet->gnet->nng), "invalid neuron group");
    nmsim_elem_neuron_ix_t ine = enet->nne;
    demand(ine <= nmsim_elem_neuron_ix_MAX, "too many neurons");
    nmsim_elem_neuron_vec_expand(&(enet->neu), ine);
    enet->neu.e[ine] = (*neu);
    (enet->nne)++;
    return ine;
  }
  
nmsim_elem_synapse_ix_t nmsim_elem_net_add_synapse
  ( nmsim_elem_net_t *enet,
    nmsim_elem_synapse_t *syn
  )
  { demand((syn->isg >= 0) && (syn->isg < enet->gnet->nsg), "invalid synapse group index");
    demand((syn->ine_pre >= 0) && (syn->ine_pre < enet->nne), "invalid pre-synaptic neuron index");
    demand((syn->ine_pos >= 0) && (syn->ine_pos < enet->nne), "invalid post_synaptic neuron index");
    nmsim_elem_synapse_ix_t ise = enet->nse;
    demand(ise <= nmsim_elem_synapse_ix_MAX, "too many synapses");
    nmsim_elem_synapse_vec_expand(&(enet->syn), ise);
    enet->syn.e[ise] = (*syn);
    (enet->nse)++;
    return ise;
  }

void nmsim_elem_net_expand_neuron_group
  ( nmsim_elem_net_t *enet,
    nmsim_group_neuron_ix_t ing,
    nmsim_elem_neuron_count_t nne_g
  )
  {
    nmsim_elem_neuron_ix_t ine = enet->nne; /* Expected index of next neuron. */
    for (nmsim_elem_neuron_ix_t k = 0; k < nne_g; k++)
      { nmsim_elem_neuron_t neu = (nmsim_elem_neuron_t) { .ing = ing };
        nmsim_elem_neuron_ix_t jne = nmsim_elem_net_add_neuron(enet, &neu);
        assert(jne == ine);
        ine++;
      }
  }
  
void nmsim_elem_net_expand_synapse_group
  ( nmsim_elem_net_t *enet,
    nmsim_group_synapse_ix_t isg,
    nmsim_elem_neuron_ix_t ine0,
    nmsim_elem_neuron_ix_t ine1,
    nmsim_elem_neuron_ix_t jne0,
    nmsim_elem_neuron_ix_t jne1
  )
  {
    nmsim_group_synapse_t *sgrp = &(enet->gnet->sgrp.e[isg]);
    nmsim_class_synapse_t *sclass = enet->gnet->cnet->sclass.e[sgrp->isc];
    double K = sgrp->K; /* Average indegre for this synapse group. */
    double W_avg = sclass->W_avg;
    double W_dev = sclass->W_dev;
    if ((K > 0.0) && (jne0 <= jne1))
      { demand(ine0 <= ine1, "presynaptic neuron group has zero neurons");
        double q = K/(K + 1);  /* Probability of adding one more synapse. */
        nmsim_elem_synapse_ix_t ise = enet->nse; /* Expected index of next synapse. */
        for (nmsim_elem_neuron_ix_t jne = jne0; jne <= jne1; jne++)
          { /* Generate synapses into neuron {jne}: */
            while (drandom() <= q)
              { /* Add one more synapse into {jne}: */
                nmsim_elem_synapse_t syn = nmsim_elem_synapse_throw
                  ( isg, isg, ine0, ine1, jne, jne, W_avg, W_dev);
                nmsim_elem_synapse_ix_t jse = nmsim_elem_net_add_synapse(enet, &syn);
                assert(jse == ise);
                ise++;
              }
          }
      }
  }
  
void nmsim_elem_net_write(FILE *wr, nmsim_elem_net_t *enet, double timeStep)
  {
    char *ind1 = "";    /* Indentation for whole description. */
    char *ind2 = "  ";  /* Indentation for items. */
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_elem_net_FILE_TYPE, nmsim_elem_net_VERSION);
    
    /* Write the neuron and synapse groups: */
    nmsim_group_net_write(wr, enet->gnet, timeStep);
    
    /* Write the neurons: */
    fprintf(wr, "%sneuron_elems = %d\n", ind2, enet->nne);
    for (nmsim_elem_neuron_ix_t ine = 0; ine < enet->nne; ine++)
      { nmsim_elem_neuron_write(wr, ine, &(enet->neu.e[ine])); fputs("\n", wr); }
    
    /* Write the synapses: */
    fprintf(wr, "%ssynapse_elems = %d\n", ind2, enet->nse);
    for (nmsim_elem_synapse_ix_t ise = 0; ise < enet->nse; ise++)
      { nmsim_elem_synapse_write(wr, ise, &(enet->syn.e[ise])); fputs("\n", wr); }
    fputs(ind1, wr);
    filefmt_write_footer(wr, nmsim_elem_net_FILE_TYPE);

    fflush(wr);
  }

nmsim_elem_net_t *nmsim_elem_net_read(FILE *rd, double timeStep)
  { 
    /* Read header line: */
    filefmt_read_header(rd, nmsim_elem_net_FILE_TYPE, nmsim_elem_net_VERSION);
    
    /* Read the neuron and synapse groups:*/
    nmsim_group_net_t *gnet = nmsim_group_net_read(rd, timeStep);
    
    /* Allocate empty net: */
    nmsim_elem_net_t *enet = nmsim_elem_net_new(gnet); 
      
    /* Vector of neuron counts per group, to check with declared ones: */
    nmsim_group_neuron_count_t nng = gnet->nng;
    nmsim_elem_neuron_count_t *nne_g = notnull(malloc(nng*sizeof(nmsim_elem_neuron_count_t)), "no mem");
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++){ nne_g[ing] = 0; }
    
    /* Read the neurons: */
    nmsim_elem_neuron_count_t nne = (nmsim_elem_neuron_count_t)
      nmsim_read_int64_param(rd, "neuron_elems", 0, nmsim_elem_neuron_count_MAX);
    fprintf(stderr, "neuron_elems = %d\n", nne);
      
    /* Preallocate the neuron vector for faster reading: */
    nmsim_elem_neuron_vec_expand(&(enet->neu), nne-1);
    
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++)
      { nmsim_elem_neuron_t neu = nmsim_elem_neuron_read(rd, ine, nng - 1);
        fget_eol(rd);
        nmsim_elem_neuron_ix_t jne = nmsim_elem_net_add_neuron(enet, &neu);
        assert(jne == ine);
        /* Count neuron as part of its group: */
        nmsim_group_neuron_ix_t ing = neu.ing;  /* Neuron's group index. */
        demand(nne_g[ing] < nmsim_elem_neuron_count_MAX, "too many neurons in group");
        nne_g[ing]++;
      }
      
    /* Check neuron counts per group: */
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { demand (gnet->ngrp.e[ing].nne == nne_g[ing], "neuron count in group does not match"); }
    
    /* Read the synapses: */
    nmsim_elem_synapse_count_t nse = (nmsim_elem_synapse_count_t)
      nmsim_read_int64_param(rd, "synapse_elems", 0, nmsim_elem_synapse_count_MAX);
    fprintf(stderr, "synapse_elems = %d\n", nse);

    if (nse > 0)
      { 
        /* Preallocate the synapse vector for faster reading: */
        nmsim_elem_synapse_vec_expand(&(enet->syn), nse-1);

        for (nmsim_elem_synapse_ix_t ise = 0; ise < nse; ise++)
          { nmsim_elem_synapse_t syn = nmsim_elem_synapse_read(rd, ise, gnet->nsg - 1, enet->nne - 1);
            fget_eol(rd);
            nmsim_elem_synapse_ix_t jse = nmsim_elem_net_add_synapse(enet, &syn);
            assert(jse == ise);
            /* Check if groups of neuron match the synapse group attributes: */
            nmsim_group_synapse_t *sgrp = &(gnet->sgrp.e[syn.isg]);
            demand(enet->neu.e[syn.ine_pre].ing == sgrp->ing_pre, "pre-synaptic neuron group mismatch");
            demand(enet->neu.e[syn.ine_pos].ing == sgrp->ing_pos, "post-synaptic neuron group mismatch");
          }
       }

    /* Read footer line: */
    filefmt_read_footer(rd, nmsim_elem_net_FILE_TYPE);

    return enet;
  }

void nmsim_elem_net_compare(nmsim_elem_net_t *enet_read, nmsim_elem_net_t *enet_orig)
  {
    assert(enet_read->gnet != NULL);
    assert(enet_orig->gnet != NULL);
    if (enet_read->gnet != enet_orig->gnet) { nmsim_group_net_compare(enet_read->gnet, enet_orig->gnet); }
    
    /* Compare neuron elems: */
    nmsim_compare_int64_param("neuron_elems", enet_read->nne, enet_orig->nne);
    assert(enet_read->nne <= enet_read->neu.ne);
    for (nmsim_elem_neuron_ix_t ine = 0; ine < enet_orig->nne; ine++)
      { nmsim_elem_neuron_compare(&(enet_read->neu.e[ine]), &(enet_orig->neu.e[ine])); }

    /* Compare synapse elems: */
    nmsim_compare_int64_param("synapse_elems", enet_read->nse, enet_orig->nse);
    assert(enet_read->nse <= enet_read->syn.ne);
    for (nmsim_elem_synapse_ix_t ise = 0; ise < enet_orig->nse; ise++)
      {  nmsim_elem_synapse_compare(&(enet_read->syn.e[ise]), &(enet_orig->syn.e[ise])); }
  }
