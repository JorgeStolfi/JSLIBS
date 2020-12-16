/* See {nmsim_elem_net.h} */
/* Last edited on 2020-12-13 16:07:06 by jstolfi */

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

nmsim_elem_net_t *nmsim_elem_net_new
  ( nmsim_group_net_t *gnet, 
    nmsim_elem_neuron_count_t nne,
    nmsim_elem_synapse_count_t nse
  )
  { 
    nmsim_elem_net_t *enet = notnull(malloc(sizeof(nmsim_elem_net_t)), "no mem");
    (*enet) = (nmsim_elem_net_t)
      { .gnet = gnet,
        .nne = nne, 
        .nse = nse,
        .neu = notnull(malloc(nne*sizeof(nmsim_elem_neuron_t)), "no mem"),
        .syn = (nse == 0 ? NULL : notnull(malloc(nse*sizeof(nmsim_elem_synapse_t)), "no mem"))
      };
    return enet;
  }

void nmsim_elem_net_free(nmsim_elem_net_t *enet)    
  {
    if (enet != NULL)
      { if (enet->gnet != NULL) { nmsim_group_net_free(enet->gnet); }
        if (enet->neu != NULL) { free(enet->neu); }
        if (enet->syn != NULL) { free(enet->syn); }
        free(enet);
      }
  }
     
void nmsim_elem_net_write(FILE *wr, nmsim_elem_net_t *enet, double timeStep)
  {
    char *ind1 = "";    /* Indentation for whole description. */
    char *ind2 = "  ";  /* Indentation for items. */
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_elem_net_FILE_TYPE, nmsim_elem_net_VERSION);
    
    /* Write the neuron and synapse counts: */
    fprintf(wr, "%sneuron_elems = %d\n", ind2, enet->nne);
    fprintf(wr, "%ssynapse_elems = %d\n", ind2, enet->nse);

    /* Write the group-level description: */
    nmsim_group_net_write(wr, enet->gnet, timeStep);
 
     /* Write the neurons: */
    for (nmsim_elem_neuron_ix_t ine = 0; ine < enet->nne; ine++)
      { nmsim_elem_neuron_write(wr, ine, &(enet->neu[ine])); fputs("\n", wr); }
    
    /* Write the synapses: */
    for (nmsim_elem_synapse_ix_t ise = 0; ise < enet->nse; ise++)
      { nmsim_elem_synapse_write(wr, ise, &(enet->syn[ise])); fputs("\n", wr); }
   
    fputs(ind1, wr);
    filefmt_write_footer(wr, nmsim_elem_net_FILE_TYPE);

    fflush(wr);
  }

nmsim_elem_net_t *nmsim_elem_net_read(FILE *rd, double timeStep)
  { 
    bool_t debug = FALSE;
    
    fprintf(stderr, "{nmsim_elem_net_read} begin\n");
    
    /* Read header line: */
    filefmt_read_header(rd, nmsim_elem_net_FILE_TYPE, nmsim_elem_net_VERSION);
    
    fprintf(stderr, "{nmsim_elem_net_read} reading element counts ...\n");

    /* Read the number of neurons:  */
    nmsim_elem_neuron_count_t nne = (nmsim_elem_neuron_count_t)
      nmsim_read_int64_param(rd, "neuron_elems", 0, nmsim_elem_neuron_count_MAX);
    
    /* Read the number of synapses:  */
    nmsim_elem_synapse_count_t nse = (nmsim_elem_synapse_count_t)
      nmsim_read_int64_param(rd, "synapse_elems", 0, nmsim_elem_synapse_count_MAX);
    
    fprintf(stderr, "{nmsim_elem_net_read} nne = %d nse = %d\n", nne, nse);

    /* Read the neuron and synapse groups:*/
    nmsim_group_net_t *gnet = nmsim_group_net_read(rd, nne, nse, timeStep);
    
    /* Allocate network: */
    nmsim_elem_net_t *enet = nmsim_elem_net_new(gnet, nne, nse); 
      
    /* Vector of neuron counts per group, to check with declared ones: */
    nmsim_group_neuron_count_t nng = gnet->nng;
    nmsim_elem_neuron_count_t *nne_g = notnull(malloc(nng*sizeof(nmsim_elem_neuron_count_t)), "no mem");
    
    /* Read the neurons: */
    fprintf(stderr, "{nmsim_elem_net_read} reading neuron elements ...\n");
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++) { nne_g[ing] = 0; }
    nmsim_elem_synapse_count_t nse_out_tot = 0; /* Sum of {.nse_out} so far. */
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++)
      { if (debug) { fprintf(stderr, "reading neuron %d - %d synapses so far\n", ine, nse_out_tot); }
        nmsim_elem_neuron_t neu = nmsim_elem_neuron_read(rd, ine, nng - 1, nse - nse_out_tot);
        fget_eol(rd);
        if (debug) { fprintf(stderr, "read neuron %d - %d outputs\n", ine, neu.nse_out); }
        /* Store in neuron table: */
        neu.ise_out_start = nse_out_tot;
        nse_out_tot += neu.nse_out;
        demand(nse_out_tot <= enet->nse, "too many output synapses");
        enet->neu[ine] = neu;
        /* Count neuron as part of its group: */
        nmsim_group_neuron_ix_t ing = neu.ing;  /* Neuron's group index. */
        demand(nne_g[ing] < nmsim_elem_neuron_count_MAX, "too many neurons in group");
        nne_g[ing]++;
      }
      
    /* Check if the sum of output synapses per neuron matches the stated total: */
    demand(nse_out_tot == nse, "number of output synapses does not match");
      
    /* Check neuron counts per group: */
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { demand (gnet->ngrp[ing].nne == nne_g[ing], "neuron count in group does not match"); }
    
    if (nse > 0)
      { 
        fprintf(stderr, "{nmsim_elem_net_read} reading synapse elements ...\n");
        
        /* Vector of output synapses per neuron, to sort the synapses: */
        nmsim_elem_synapse_count_t *nse_out_n = 
          notnull(malloc(nne*sizeof(nmsim_elem_synapse_count_t)), "no mem");
        for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) { nse_out_n[ine] = 0; }

        /* Read the synapses: */
        for (nmsim_elem_synapse_ix_t ise = 0; ise < nse; ise++)
          { nmsim_elem_synapse_t syn = nmsim_elem_synapse_read
              ( rd, ise, gnet->nsg - 1, enet->nne - 1 );
            fget_eol(rd);
            
            /* Choose the new synapse index {jse} to properly sort the synapses: */
            nmsim_elem_neuron_ix_t ine_pre = syn.ine_pre;
            nmsim_elem_neuron_t *neu_pre = &(enet->neu[ine_pre]);
            demand(nse_out_n[ine_pre] < neu_pre->nse_out, "too many synapses for neuron");
            nmsim_elem_synapse_ix_t jse = neu_pre->ise_out_start + nse_out_n[ine_pre];
            assert((jse >= 0) && (jse < nse));
            enet->syn[jse] = syn;
            nse_out_n[ine_pre]++;
            
            /* Check if groups of neuron match the synapse group attributes: */
            nmsim_group_synapse_t *sgrp = &(gnet->sgrp[syn.isg]);
            demand
              ( enet->neu[syn.ine_pre].ing == sgrp->ing_pre, 
                "pre-synaptic neuron group mismatch"
              );
            demand
              ( enet->neu[syn.ine_pos].ing == sgrp->ing_pos, 
                "post-synaptic neuron group mismatch"
              );
          }
          
        /* Paranoia: check number of output synapses per neuron: */
        for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++)
          { nmsim_elem_neuron_t *neu = &(enet->neu[ine]);
            assert(nse_out_n[ine] == neu->nse_out);
          }
      }

    /* Read footer line: */
    fprintf(stderr, "{nmsim_elem_net_read} reading footer ...\n");
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
    for (nmsim_elem_neuron_ix_t ine = 0; ine < enet_orig->nne; ine++)
      { nmsim_elem_neuron_compare(&(enet_read->neu[ine]), &(enet_orig->neu[ine])); }

    /* Compare synapse elems: */
    nmsim_compare_int64_param("synapse_elems", enet_read->nse, enet_orig->nse);
    for (nmsim_elem_synapse_ix_t ise = 0; ise < enet_orig->nse; ise++)
      {  nmsim_elem_synapse_compare(&(enet_read->syn[ise]), &(enet_orig->syn[ise])); }
  }
