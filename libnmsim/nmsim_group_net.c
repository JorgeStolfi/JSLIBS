/* See {nmsim_group_net.h} */
/* Last edited on 2020-12-11 19:55:35 by jstolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <jsmath.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_synapse.h>

#include <nmsim_group_net.h>

nmsim_group_net_t *nmsim_group_net_new
  ( nmsim_class_net_t *cnet, 
    nmsim_group_neuron_count_t nng,
    nmsim_group_synapse_count_t nsg
  )
  {
    nmsim_group_net_t *gnet = notnull(malloc(sizeof(nmsim_group_net_t)), "no mem");
    demand(cnet != NULL, "class-level description is {NULL}");
    demand((nng >= 1) && (nng <= nmsim_group_neuron_count_MAX), "invalid neuron group count");
    demand((nsg >= 0) && (nng <= nmsim_group_synapse_count_MAX), "invalid synapse group count");
    (*gnet) = (nmsim_group_net_t)
      { .cnet = cnet,
        .nng = nng, 
        .nsg = nsg,
        .ngrp = notnull(malloc(nng*sizeof(nmsim_group_neuron_t)), "no mem"),
        .sgrp = (nsg == 0 ? NULL : notnull(malloc(nsg*sizeof(nmsim_group_synapse_t)), "no mem"))
      };
    return gnet;
  }

void nmsim_group_net_free(nmsim_group_net_t *gnet)
  {
    if (gnet != NULL)
      { if (gnet->cnet != NULL) { nmsim_class_net_free(gnet->cnet); }
        if (gnet->ngrp != NULL) { free(gnet->ngrp); }
        if (gnet->sgrp != NULL) { free(gnet->sgrp); }
        free(gnet);
      }
  }
  
void nmsim_group_net_write(FILE *wr, nmsim_group_net_t *gnet, double timeStep)
  { 
    char *ind1 = "  ";    /* Indentation for whole description. */
    char *ind2 = "    ";  /* Indentation for items. */
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_group_net_FILE_TYPE, nmsim_group_net_VERSION);
    
    /* Write the neuron and synapse group counts: */
    fprintf(wr, "%sneuron_groups = %d\n", ind2, gnet->nng);
    fprintf(wr, "%ssynapse_groups = %d\n", ind2, gnet->nsg);

    /* Write the neuron and synapse classes: */
    nmsim_class_net_write(wr, gnet->cnet, timeStep);
    
    /* Write the neuron groups: */
    for (nmsim_group_neuron_ix_t ing = 0; ing < gnet->nng; ing++)
      { nmsim_group_neuron_write(wr, ing, &(gnet->ngrp[ing])); fputs("\n", wr); }
    
    /* Write the synapse groups: */
    for (nmsim_group_synapse_ix_t isg = 0; isg < gnet->nsg; isg++)
      { nmsim_group_synapse_write(wr, isg, &(gnet->sgrp[isg])); fputs("\n", wr); }

    fputs(ind1, wr);
    filefmt_write_footer(wr, nmsim_group_net_FILE_TYPE);

    fflush(wr);
  }

nmsim_group_net_t *nmsim_group_net_read
  ( FILE *rd, 
    nmsim_elem_neuron_count_t nne_g_max,
    nmsim_elem_neuron_count_t nse_g_max,
    double timeStep
  )
  { 
    bool_t debug = FALSE;
    /* Read header line: */
    filefmt_read_header(rd, nmsim_group_net_FILE_TYPE, nmsim_group_net_VERSION);
    
    /* Read the number of neuron groups: */
    nmsim_group_neuron_count_t nng = (nmsim_group_neuron_count_t)
      nmsim_read_int64_param(rd, "neuron_groups", 1, nmsim_group_neuron_count_MAX);
      
    /* Read the number of synapse groups: */
    nmsim_group_synapse_count_t nsg = (nmsim_group_synapse_count_t)
      nmsim_read_int64_param(rd, "synapse_groups", 0, nmsim_group_synapse_count_MAX);

    if (debug) { fprintf(stderr, "begin {nmsim_group_net_read} nng = %d nsg = %d\n", nng, nsg); }

    /* Read the neuron and synapse classes:*/
    nmsim_class_net_t *cnet = nmsim_class_net_read(rd, timeStep);
    
    /* Allocate group-level net: */
    nmsim_group_net_t *gnet = nmsim_group_net_new(cnet, nng, nsg); 
    
    /* Read the neuron groups: */
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_group_neuron_t ngrp = nmsim_group_neuron_read
          (rd, ing, cnet->nnc - 1, nne_g_max, nsg);
        fget_eol(rd);
        gnet->ngrp[ing] = ngrp;
      }
    
    /* Read the synapse groups: */
    if (nsg > 0)
      { assert(cnet->nsc > 0);
        for (nmsim_group_synapse_ix_t isg = 0; isg < nsg; isg++)
          { nmsim_group_synapse_t sgrp = nmsim_group_synapse_read
              ( rd, isg, cnet->nsc - 1, gnet->nng - 1, nse_g_max);
            fget_eol(rd);
            gnet->sgrp[isg] = sgrp;
          }
      }

    /* Read footer line: */
    filefmt_read_footer(rd, nmsim_group_net_FILE_TYPE);

    return gnet;
  }
    
void nmsim_group_net_compare(nmsim_group_net_t *gnet_read, nmsim_group_net_t *gnet_orig)
  { 
    assert(gnet_read->cnet != NULL);
    assert(gnet_orig->cnet != NULL);
    if (gnet_read->cnet != gnet_orig->cnet) 
      { nmsim_class_net_compare(gnet_read->cnet, gnet_orig->cnet); }
    
    /* Compare neuron groups: */
    nmsim_compare_int64_param("neuron_groups", gnet_read->nng, gnet_orig->nng);
    for (nmsim_group_neuron_ix_t ing = 0; ing < gnet_orig->nng; ing++)
      { nmsim_group_neuron_compare(&(gnet_read->ngrp[ing]), &(gnet_orig->ngrp[ing])); }

    /* Compare synapse groups: */
    nmsim_compare_int64_param("synapse_groups", gnet_read->nsg, gnet_orig->nsg);
    for (nmsim_group_synapse_ix_t isg = 0; isg < gnet_orig->nsg; isg++)
      {  nmsim_group_synapse_compare(&(gnet_read->sgrp[isg]), &(gnet_orig->sgrp[isg])); }
  }    
