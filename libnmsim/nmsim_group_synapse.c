/* See {nmsim_group_synapse.h} */
/* Last edited on 2020-12-11 17:47:33 by jstolfi */

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
#include <nmsim_group_synapse.h> 

void nmsim_group_synapse_write(FILE *wr, nmsim_group_synapse_ix_t isg, nmsim_group_synapse_t *sgrp)
  {
    char *ind1 = "    ";   /* Indent of group line. */
    fprintf(wr, "%s%d %d", ind1, isg, sgrp->isc);
    fprintf(wr, " %d %d", sgrp->ing_pre, sgrp->ing_pos);
    fprintf(wr, " %d", sgrp->nse);
  }

nmsim_group_synapse_t nmsim_group_synapse_read
  ( FILE *rd, 
    nmsim_group_synapse_ix_t isg, 
    nmsim_class_synapse_ix_t isc_max, 
    nmsim_group_neuron_ix_t ing_max, 
    nmsim_elem_synapse_count_t nse_g_max
  )
  { 
    (void)nmsim_read_int64_value(rd, "synapse group index", isg, isg);
    
    nmsim_class_synapse_ix_t isc = (nmsim_class_synapse_ix_t)
      nmsim_read_int64_value(rd, "synapse class index", 0, isc_max);
    
    nmsim_group_neuron_ix_t ing_pre = (nmsim_group_neuron_ix_t)
      nmsim_read_int64_value(rd, "pre-synaptic neuron group index", 0, ing_max);
    nmsim_group_neuron_ix_t ing_pos = (nmsim_group_neuron_ix_t)
      nmsim_read_int64_value(rd, "post-synaptic neuron group index", 0, ing_max);
    
    nmsim_elem_synapse_count_t nse_g = (nmsim_elem_synapse_count_t)
      nmsim_read_int64_value(rd, "number of synapses in group", 0, nse_g_max);
    
    nmsim_group_synapse_t sgrp = (nmsim_group_synapse_t)
      { .isc = isc, .ing_pre = ing_pre, .ing_pos = ing_pos, .nse = nse_g };

    return sgrp;
  }

void nmsim_group_synapse_show(FILE *wr, char *pref, nmsim_group_synapse_t *sgrp, char *suff)
  { 
    if (pref != NULL) { fputs(pref, wr); }
    fprintf(wr, "{ class = %d", sgrp->isc);
    fprintf(wr, " ing_pre = %d ing_pos = %d", sgrp->ing_pre, sgrp->ing_pos);
    fprintf(wr, " nse = %d", sgrp->nse);
    fprintf(wr, " }");
    if (suff != NULL) { fputs(suff, wr); }
  }
  
nmsim_group_synapse_t nmsim_group_synapse_throw
  ( nmsim_class_synapse_ix_t isc_max, 
    nmsim_group_neuron_ix_t ing_pre, 
    nmsim_group_neuron_ix_t ing_pos,
    nmsim_elem_synapse_count_t nse_min,
    nmsim_elem_synapse_count_t nse_max
  )
  {
    /* Choose the synapse class {isc}: */
    demand((isc_max >= 0) && (isc_max <= nmsim_class_synapse_ix_MAX), "invalid max synapse class index");
    nmsim_class_synapse_ix_t isc = (nmsim_class_synapse_ix_t)int64_abrandom(0, isc_max);
    
    /* Check the pre- and post-synaptic group indices: */
    demand((ing_pre >= 0) && (ing_pre <= nmsim_group_neuron_ix_MAX), "invalid pre-synaptic neuron group index");
    demand((ing_pos >= 0) && (ing_pos <= nmsim_group_neuron_ix_MAX), "invalid post-synaptic neuron group index");

    /* Choose the number of synapses {nse}: */
    demand((nse_min >= 0) && (nse_max <= nmsim_elem_synapse_count_MAX), "invalid min/max synapse count in group");
    demand((nse_min <= nse_max), "invalid synapse count range in group");
    nmsim_elem_synapse_count_t nse_g = (nmsim_class_synapse_count_t)int64_abrandom(nse_min, nse_max);

    /* Package it: */
    nmsim_group_synapse_t sgrp = (nmsim_group_synapse_t)
      { .isc = isc, .ing_pre = ing_pre, .ing_pos = ing_pos, .nse = nse_g };
    return sgrp;
  }

void nmsim_group_synapse_compare(nmsim_group_synapse_t *sgrp_read, nmsim_group_synapse_t *sgrp)
  {
    nmsim_compare_int64_param("isc",  sgrp_read->isc, sgrp->isc);
    nmsim_compare_int64_param("ing_pre",  sgrp_read->ing_pre, sgrp->ing_pre);
    nmsim_compare_int64_param("ing_pos",  sgrp_read->ing_pos, sgrp->ing_pos);
    nmsim_compare_int64_param("nse",  sgrp_read->nse, sgrp->nse);
  }

