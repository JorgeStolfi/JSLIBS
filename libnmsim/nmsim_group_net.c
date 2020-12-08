/* See {nmsim_group_net.h} */
/* Last edited on 2019-06-27 17:07:43 by jstolfi */

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
#include <jsrandom.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_synapse.h>

#include <nmsim_group_net.h>

nmsim_group_net_t *nmsim_group_net_new(nmsim_class_net_t *cnet)
  {
    nmsim_group_net_t *gnet = notnull(malloc(sizeof(nmsim_group_net_t)), "no mem");
    demand(cnet != NULL, "class-level description is {NULL}");
    (*gnet) = (nmsim_group_net_t)
      { .cnet = cnet,
        .nng = 0, .nsg = 0,
        .ngrp = nmsim_group_neuron_vec_new(20),
        .sgrp = nmsim_group_synapse_vec_new(20)
      };
    return gnet;
  }

void nmsim_group_net_free(nmsim_group_net_t *gnet)
  {
    if (gnet != NULL)
      { if (gnet->cnet != NULL) { nmsim_class_net_free(gnet->cnet); }
        if (gnet->ngrp.e != NULL) { free(gnet->ngrp.e); }
        if (gnet->sgrp.e != NULL) { free(gnet->sgrp.e); }
        free(gnet);
      }
  }
  
nmsim_group_neuron_ix_t nmsim_group_net_add_neuron_group
  ( nmsim_group_net_t *gnet, 
    nmsim_group_neuron_t *ngrp
  )
  { demand((ngrp->inc >= 0) && (ngrp->inc < gnet->cnet->nnc), "invalid neuron class");
    demand((ngrp->nne >= 0) && (ngrp->nne <= nmsim_elem_neuron_count_MAX), "invalid neuron count");
    nmsim_group_neuron_ix_t ing = gnet->nng;
    demand(ing <= nmsim_group_neuron_ix_MAX, "too many neuron  groups");
    nmsim_group_neuron_vec_expand(&(gnet->ngrp), ing);
    gnet->ngrp.e[ing] = (*ngrp);
    (gnet->nng)++;
    return ing;
  }
  
nmsim_group_synapse_ix_t nmsim_group_net_add_synapse_group
  ( nmsim_group_net_t *gnet,
    nmsim_group_synapse_t *sgrp
  )
  { demand((sgrp->isc >= 0) && (sgrp->isc < gnet->cnet->nsc), "invalid synapse class index");
    demand((sgrp->ing_pre >= 0) && (sgrp->ing_pre < gnet->nng), "invalid pre-synaptic neuron group index");
    demand((sgrp->ing_pos >= 0) && (sgrp->ing_pos < gnet->nng), "invalid post_synaptic neuron group index");
    demand(sgrp->K >= 0.0, "invalid in-degree");
    nmsim_group_synapse_ix_t isg = gnet->nsg;
    demand(isg <= nmsim_group_synapse_ix_MAX, "too many synapse groups");
    nmsim_group_synapse_vec_expand(&(gnet->sgrp), isg);
    gnet->sgrp.e[isg] = (*sgrp);
    (gnet->nsg)++;
    return isg;
  }

void nmsim_group_net_write(FILE *wr, nmsim_group_net_t *gnet, double timeStep)
  { 
    char *ind1 = "  ";    /* Indentation for whole description. */
    char *ind2 = "    ";  /* Indentation for items. */
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_group_net_FILE_TYPE, nmsim_group_net_VERSION);
    
    /* Write the neuron and synapse classes: */
    nmsim_class_net_write(wr, gnet->cnet, timeStep);
    
    /* Write the neuron groups: */
    fprintf(wr, "%sneuron_groups = %d\n", ind2, gnet->nng);
    for (nmsim_group_neuron_ix_t ing = 0; ing < gnet->nng; ing++)
      { nmsim_group_neuron_write(wr, ing, &(gnet->ngrp.e[ing])); fputs("\n", wr); }
    
    /* Write the synapse groups: */
    fprintf(wr, "%ssynapse_groups = %d\n", ind2, gnet->nsg);
    for (nmsim_group_synapse_ix_t isg = 0; isg < gnet->nsg; isg++)
      { nmsim_group_synapse_write(wr, isg, &(gnet->sgrp.e[isg])); fputs("\n", wr); }
    fputs(ind1, wr);
    filefmt_write_footer(wr, nmsim_group_net_FILE_TYPE);

    fflush(wr);
  }

nmsim_group_net_t *nmsim_group_net_read(FILE *rd, double timeStep)
  { 
    /* Read header line: */
    filefmt_read_header(rd, nmsim_group_net_FILE_TYPE, nmsim_group_net_VERSION);
    
    /* Read the neuron and synapse classes:*/
    nmsim_class_net_t *cnet = nmsim_class_net_read(rd, timeStep);
    
    /* Allocate empty net: */
    nmsim_group_net_t *gnet = nmsim_group_net_new(cnet); 
    
    /* Read the neuron groups: */
    nmsim_group_neuron_count_t nng = (nmsim_group_neuron_count_t)
      nmsim_read_int64_param(rd, "neuron_groups", 0, nmsim_group_neuron_count_MAX);
      
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_group_neuron_t ngrp = nmsim_group_neuron_read(rd, ing, cnet->nnc - 1);
        fget_eol(rd);
        nmsim_group_neuron_ix_t jng = nmsim_group_net_add_neuron_group(gnet, &ngrp);
        assert(jng == ing);
      }
    
    /* Read the synapse groups: */
    nmsim_group_synapse_count_t nsg = (nmsim_group_synapse_count_t)
      nmsim_read_int64_param(rd, "synapse_groups", 0, nmsim_group_synapse_count_MAX);
    for (nmsim_group_synapse_ix_t isg = 0; isg < nsg; isg++)
      { nmsim_group_synapse_t sgrp = nmsim_group_synapse_read(rd, isg, cnet->nsc - 1, gnet->nng - 1);
        fget_eol(rd);
        nmsim_group_synapse_ix_t jsg = nmsim_group_net_add_synapse_group(gnet, &sgrp);
        assert(jsg == isg);
      }

    /* Read footer line: */
    filefmt_read_footer(rd, nmsim_group_net_FILE_TYPE);

    return gnet;
  }
    
nmsim_group_net_t *nmsim_group_net_throw
  ( nmsim_class_net_t *cnet,
    nmsim_group_neuron_count_t nng, 
    nmsim_group_synapse_count_t nsg,
    nmsim_elem_neuron_count_t nne_g_min,
    nmsim_elem_neuron_count_t nne_g_max,
    double K_min,
    double K_max
  )
  {
    demand((nng >= 0) && (nng <= nmsim_group_neuron_count_MAX), "invalid neuron group count");
    demand((nsg >= 0) && (nsg <= nmsim_group_synapse_count_MAX), "invalid synapse group count");

    demand((nne_g_min >= 0) && (nne_g_max <= nmsim_elem_neuron_count_MAX), "invalid min/max neuron counts per group");
    demand(nne_g_max >= nne_g_min, "invalid neuron count range per group");

    demand((K_min >= 0) && (K_max <= (double)nne_g_max), "invalid min/max avg indegree");
    demand(K_max >= K_min, "invalid avg indegree range");

    nmsim_group_net_t *gnet = nmsim_group_net_new(cnet);

    /* Add neuron groups: */
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++)
      { nmsim_group_neuron_t ngrp = nmsim_group_neuron_throw(cnet->nnc - 1, nne_g_min, nne_g_max);
        nmsim_group_neuron_ix_t jng = nmsim_group_net_add_neuron_group(gnet, &ngrp);
        assert(jng == ing);
      }
    
    /* Add synapse groups: */
    for (nmsim_group_synapse_ix_t isg = 0; isg < nsg; isg++)
      { nmsim_group_neuron_ix_t ing_pre = (nmsim_group_neuron_ix_t)int64_abrandom(0, nng-1);
        nmsim_group_neuron_ix_t ing_pos = (nmsim_group_neuron_ix_t)int64_abrandom(0, nng-1);
        /* Limit the indegree to the size of the pre-synaptic group: */
        nmsim_elem_neuron_count_t nne_pre = gnet->ngrp.e[ing_pre].nne; 
        double K_g_min = fmin(K_min, (double)nne_pre);
        double K_g_max = fmin(K_max, (double)nne_pre);
        nmsim_group_synapse_t sgrp = nmsim_group_synapse_throw(cnet->nsc - 1, ing_pre, ing_pos, K_g_min, K_g_max);
        nmsim_group_synapse_ix_t jsg = nmsim_group_net_add_synapse_group(gnet, &sgrp);
        assert(jsg == isg);
      }
    return gnet;
  }

void nmsim_group_net_compare(nmsim_group_net_t *gnet_read, nmsim_group_net_t *gnet_orig)
  { 
    assert(gnet_read->cnet != NULL);
    assert(gnet_orig->cnet != NULL);
    if (gnet_read->cnet != gnet_orig->cnet) { nmsim_class_net_compare(gnet_read->cnet, gnet_orig->cnet); }
    
    /* Compare neuron groups: */
    nmsim_compare_int64_param("neuron_groups", gnet_read->nng, gnet_orig->nng);
    assert(gnet_read->nng <= gnet_read->ngrp.ne);
    for (nmsim_group_neuron_ix_t ing = 0; ing < gnet_orig->nng; ing++)
      { nmsim_group_neuron_compare(&(gnet_read->ngrp.e[ing]), &(gnet_orig->ngrp.e[ing])); }

    /* Compare synapse groups: */
    nmsim_compare_int64_param("synapse_groups", gnet_read->nsg, gnet_orig->nsg);
    assert(gnet_read->nsg <= gnet_read->sgrp.ne);
    for (nmsim_group_synapse_ix_t isg = 0; isg < gnet_orig->nsg; isg++)
      {  nmsim_group_synapse_compare(&(gnet_read->sgrp.e[isg]), &(gnet_orig->sgrp.e[isg])); }
  }    
