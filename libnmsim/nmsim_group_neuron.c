/* See {nmsim_group_neuron.h} */
/* Last edited on 2020-12-12 16:32:45 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <fget.h>
#include <vec.h>
#include <jsmath.h>
#include <jsrandom.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_group_neuron.h>
 
void nmsim_group_neuron_write(FILE *wr, nmsim_group_neuron_ix_t ing, nmsim_group_neuron_t *ngrp)
  {
    char *ind1 = "    ";   /* Indent of group line. */
    fprintf(wr, "%s%d %d %d %d", ind1, ing, ngrp->inc, ngrp->nne, ngrp->nsg_out);
    /* Do not write {ngrp->ine_start}. */
  }

nmsim_group_neuron_t nmsim_group_neuron_read
  ( FILE *rd, 
    nmsim_group_neuron_ix_t ing, 
    nmsim_class_neuron_ix_t inc_max,
    nmsim_elem_neuron_count_t nne_g_max,
    nmsim_group_synapse_count_t nsg_out_max
  )
  {
    bool_t debug = FALSE;
    
    if (debug) { fprintf(stderr, "    reading neuron group %d ... ", ing); }
        
    (void)nmsim_read_int64_value(rd, "neuron group index", ing, ing);
    
    nmsim_class_neuron_ix_t inc = (nmsim_class_neuron_ix_t)
      nmsim_read_int64_value(rd, "neuron class", 0, inc_max);
    
    nmsim_elem_neuron_count_t nne_g = (nmsim_elem_neuron_count_t)
      nmsim_read_int64_value(rd, "neurons in group", 1, nne_g_max);

    nmsim_group_synapse_count_t nsg_out = (nmsim_group_synapse_count_t)
      nmsim_read_int64_value(rd, "synaptic bundles out of group", 0, nsg_out_max);

    nmsim_group_neuron_t ngrp = (nmsim_group_neuron_t)
      { .inc = inc, .nne = nne_g, .ine_start = 0, .nsg_out = nsg_out };

    if (debug) { nmsim_group_neuron_show(stderr, NULL, &ngrp, "\n"); }
        
    return ngrp;
  }

void nmsim_group_neuron_show(FILE *wr, char *pref, nmsim_group_neuron_t *ngrp, char *suff)
  { 
    if (pref != NULL) { fputs(pref, wr); }
    fprintf(wr, "{ class = %d  neurons = %d", ngrp->inc, ngrp->nne);
    if (ngrp->ine_start != 0) 
      { fprintf(wr, " (%d..%d)", ngrp->ine_start, ngrp->ine_start + ngrp->nne - 1); }
    fprintf(wr, " output bundles = %d", ngrp->nsg_out);
    fprintf(wr, " }");
    if (suff != NULL) { fputs(suff, wr); }
  }

nmsim_group_neuron_t nmsim_group_neuron_throw
  ( nmsim_class_neuron_ix_t inc_max,
    nmsim_elem_neuron_count_t nne_min, 
    nmsim_elem_neuron_count_t nne_max,
    nmsim_group_synapse_count_t nsg_out_min, 
    nmsim_group_synapse_count_t nsg_out_max
  )
  { /* Pick a class: */
    demand((inc_max >= 0) && (inc_max <= nmsim_class_neuron_ix_MAX), "invalid neuron class range");
    nmsim_class_neuron_ix_t inc = (nmsim_class_neuron_ix_t)int64_abrandom(0, inc_max);
    
    /* Can we have at least one synapse per synaptic group? */
    demand(nsg_out_min <= nne_max, "too few neurons for this many synaptic bundles");
    
    /* Pick a neuron count: */
    demand((nne_min >= 1) && (nne_max <= nmsim_elem_neuron_count_MAX), "invalid neuron count in group");
    demand((nne_min <= nne_max), "invalid neuron count range in group");
    nmsim_elem_neuron_count_t nne_g = (nmsim_elem_neuron_count_t)int64_abrandom(nne_min, nne_max);
    
    /* Pick a number of output synaptic bundles: */
    demand
      ( (nsg_out_min >= 0) && (nsg_out_max <= nmsim_group_synapse_count_MAX), 
        "invalid output bundle count in group"
      );
    demand((nsg_out_min <= nsg_out_max), "invalid output synaptic bundle count range in group");
    if (nsg_out_min > nne_g) { nsg_out_min = nne_g; }
    assert(nsg_out_min <= nne_g);
    assert(nsg_out_min <= nsg_out_max);
    nmsim_group_synapse_count_t nsg_out = (nmsim_group_synapse_count_t)
      int64_abrandom(nsg_out_min, nsg_out_max);
    
    /* Pack and deliver: */
    nmsim_group_neuron_t ngrp = (nmsim_group_neuron_t)
      { .inc = inc, .nne = nne_g, .ine_start = 0, .nsg_out = nsg_out };
    return ngrp;
  }

void nmsim_group_neuron_compare(nmsim_group_neuron_t *ngrp_read, nmsim_group_neuron_t *ngrp_orig)
  {
    nmsim_compare_int64_param("class",  ngrp_read->inc, ngrp_orig->inc);
    nmsim_compare_int64_param("neuron count",  ngrp_read->nne, ngrp_orig->nne);
    /* Ignore {.ine_start}. */
  }

vec_typeimpl(nmsim_group_neuron_vec_t,nmsim_group_neuron_vec,nmsim_group_neuron_t);
