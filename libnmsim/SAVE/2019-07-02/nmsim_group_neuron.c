/* See {nmsim_group_neuron.h} */
/* Last edited on 2019-03-28 18:22:35 by jstolfi */

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
    fprintf(wr, "%s%d %d %d", ind1, ing, ngrp->inc, ngrp->nne);
  }

nmsim_group_neuron_t nmsim_group_neuron_read
  ( FILE *rd, 
    nmsim_group_neuron_ix_t ing, 
    nmsim_class_neuron_ix_t inc_max
  )
  {
    (void)nmsim_read_int64_value(rd, "neuron group index", ing, ing);
    nmsim_class_neuron_ix_t inc = (nmsim_class_neuron_ix_t)
      nmsim_read_int64_value(rd, "neuron class", 0, inc_max);
    nmsim_elem_neuron_count_t nne_g = (nmsim_elem_neuron_count_t)
      nmsim_read_int64_value(rd, "neurons in group", 0, nmsim_elem_neuron_count_MAX);

    nmsim_group_neuron_t ngrp = (nmsim_group_neuron_t){ .inc = inc, .nne = nne_g };
    return ngrp;
  }

void nmsim_group_neuron_show(FILE *wr, char *pref, nmsim_group_neuron_t *ngrp, char *suff)
  { 
    if (pref != NULL) { fputs(pref, wr); }
    fprintf(wr, "{ class = %d  neurons = %d", ngrp->inc, ngrp->nne);
    fprintf(wr, " }");
    if (suff != NULL) { fputs(suff, wr); }
  }

nmsim_group_neuron_t nmsim_group_neuron_throw
  ( nmsim_class_neuron_ix_t inc_max,
    nmsim_elem_neuron_count_t nne_g_min, 
    nmsim_elem_neuron_count_t nne_g_max
  )
  { demand((inc_max >= 0) && (inc_max <= nmsim_class_neuron_ix_MAX), "invalid neuron class range");
    demand((nne_g_min >= 0) && (nne_g_max <= nmsim_elem_neuron_count_MAX), "invalid neuron count in group");
    demand((nne_g_min <= nne_g_max), "invalid neuron count range in group");
    nmsim_class_neuron_ix_t inc = (nmsim_class_neuron_ix_t)int64_abrandom(0, inc_max);
    nmsim_elem_neuron_count_t nne_g = (nmsim_elem_neuron_count_t)int64_abrandom(nne_g_min, nne_g_max);
    nmsim_group_neuron_t ngrp = (nmsim_group_neuron_t){ .inc = inc, .nne = nne_g };
    return ngrp;
  }

void nmsim_group_neuron_compare(nmsim_group_neuron_t *ngrp_read, nmsim_group_neuron_t *ngrp_orig)
  {
    nmsim_compare_int64_param("class",  ngrp_read->inc, ngrp_orig->inc);
    nmsim_compare_int64_param("neuron count",  ngrp_read->nne, ngrp_orig->nne);
  }

vec_typeimpl(nmsim_group_neuron_vec_t,nmsim_group_neuron_vec,nmsim_group_neuron_t);
