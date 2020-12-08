/* See {nmsim_group_synapse.h} */
/* Last edited on 2019-05-28 14:41:21 by jstolfi */

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
    fprintf(wr, "%s%d %d  %d %d  ", ind1, isg, sgrp->isc, sgrp->ing_pre, sgrp->ing_pos);
    nmsim_write_double_value(wr, sgrp->K, nmsim_write_KW_PREC, FALSE, TRUE, FALSE);
  }

nmsim_group_synapse_t nmsim_group_synapse_read
  ( FILE *rd, 
    nmsim_group_synapse_ix_t isg, 
    nmsim_class_synapse_ix_t isc_max, 
    nmsim_group_neuron_ix_t ing_max
  )
  { 
    (void)nmsim_read_int64_value(rd, "synapse group index", isg, isg);
    nmsim_class_synapse_ix_t isc = (nmsim_class_synapse_ix_t)
      nmsim_read_int64_value(rd, "synapse class index", 0, isc_max);
    nmsim_group_neuron_ix_t ing_pre = (nmsim_group_neuron_ix_t)
      nmsim_read_int64_value(rd, "pre-synaptic neuron grou index", 0, ing_max);
    nmsim_group_neuron_ix_t ing_pos = (nmsim_group_neuron_ix_t)
      nmsim_read_int64_value(rd, "post-synaptic neuron group index", 0, ing_max);
    double Kmax = (double)nmsim_elem_neuron_count_MAX;
    double K = nmsim_read_double_value(rd, "mean in-degree of post-synaptic neurons", 0.0, Kmax);
    
    nmsim_group_synapse_t sgrp = (nmsim_group_synapse_t)
      { .isc = isc, .ing_pre = ing_pre, .ing_pos = ing_pos, .K = K };

    return sgrp;
  }

void nmsim_group_synapse_show(FILE *wr, char *pref, nmsim_group_synapse_t *sgrp, char *suff)
  { 
    if (pref != NULL) { fputs(pref, wr); }
    fprintf(wr, "{ class = %d", sgrp->isc);
    fprintf(wr, " ing_pre = %d ing_pos = %d", sgrp->ing_pre, sgrp->ing_pos);
    fprintf(wr, " K = %.3f", sgrp->K);
    fprintf(wr, " }");
    if (suff != NULL) { fputs(suff, wr); }
  }
  
nmsim_group_synapse_t nmsim_group_synapse_throw
  ( nmsim_class_synapse_ix_t isc_max, 
    nmsim_group_neuron_ix_t ing_pre, 
    nmsim_group_neuron_ix_t ing_pos,
    double K_min,
    double K_max
  )
  {
    demand((isc_max >= 0) && (isc_max <= nmsim_class_synapse_ix_MAX), "invalid max synapse class index");
    demand((ing_pre >= 0) && (ing_pre <= nmsim_group_neuron_ix_MAX), "invalid pre-synaptic neuron group index");
    demand((ing_pos >= 0) && (ing_pos <= nmsim_group_neuron_ix_MAX), "invalid post-synaptic neuron group index");

    demand((K_min >= 0) && (K_max <= nmsim_elem_neuron_count_MAX), "invalid min/max avg indegree in group");
    demand((K_min <= K_max), "invalid avg indegree range in group");
    
    nmsim_class_synapse_ix_t isc = (nmsim_class_synapse_ix_t)int64_abrandom(0, isc_max);
    double K = nmsim_throw_double(K_min, K_max);
    nmsim_group_synapse_t sgrp = (nmsim_group_synapse_t)
      { .isc = isc, .ing_pre = ing_pre, .ing_pos = ing_pos, .K = K };
    return sgrp;
  }

void nmsim_group_synapse_compare(nmsim_group_synapse_t *sgrp_read, nmsim_group_synapse_t *sgrp)
  {
    nmsim_compare_int64_param("isc",  sgrp_read->isc, sgrp->isc);
    nmsim_compare_int64_param("ing_pre",  sgrp_read->ing_pre, sgrp->ing_pre);
    nmsim_compare_int64_param("ing_pos",  sgrp_read->ing_pos, sgrp->ing_pos);
    nmsim_compare_double_param("K",  sgrp_read->K, sgrp->K, nmsim_write_KW_PREC, TRUE, FALSE);
  }

vec_typeimpl(nmsim_group_synapse_vec_t,nmsim_group_synapse_vec,nmsim_group_synapse_t);
