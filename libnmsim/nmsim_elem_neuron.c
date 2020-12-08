/* See {nmsim_elem_neuron.h} */
/* Last edited on 2019-03-28 18:17:18 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <vec.h>
#include <affirm.h>
#include <jsrandom.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_elem_neuron.h>

void nmsim_elem_neuron_write(FILE *wr, nmsim_elem_neuron_ix_t ine, nmsim_elem_neuron_t *neu)
  {
    char *ind1 = "  ";   /* Indent of neuron line. */
    fprintf(wr, "%s%d %d", ind1, ine, neu->ing);
  }

nmsim_elem_neuron_t nmsim_elem_neuron_read
  ( FILE *rd, 
    nmsim_elem_neuron_ix_t ine, 
    nmsim_group_neuron_ix_t ing_max
  )
  { 
    (void)nmsim_read_int64_value(rd, "neuron index", ine, ine);
    nmsim_group_neuron_ix_t ing = (nmsim_group_neuron_ix_t)
      nmsim_read_int64_value(rd, "neuron group index", 0, ing_max);

    nmsim_elem_neuron_t neu = (nmsim_elem_neuron_t){ .ing = ing };

    return neu;
  }

void nmsim_elem_neuron_show(FILE *wr, char *pref, nmsim_elem_neuron_t *neu, char *suff)
  { 
    if (pref != NULL) { fputs(pref, wr); }
    fprintf(wr, "{ group = %d }", neu->ing);
    if (suff != NULL) { fputs(suff, wr); }
  }
  
nmsim_elem_neuron_t nmsim_elem_neuron_throw(nmsim_group_neuron_ix_t ing_max)
  {
    nmsim_group_neuron_ix_t ing = (nmsim_group_neuron_ix_t)int64_abrandom(0, ing_max);
    nmsim_elem_neuron_t neu = (nmsim_elem_neuron_t){ .ing = ing };
    return neu;
  }

void nmsim_elem_neuron_compare(nmsim_elem_neuron_t *neu_read, nmsim_elem_neuron_t *neu_orig)
  {
    nmsim_compare_int64_param("ing",  neu_read->ing, neu_orig->ing);
  }

vec_typeimpl(nmsim_elem_neuron_vec_t, nmsim_elem_neuron_vec, nmsim_elem_neuron_t);
