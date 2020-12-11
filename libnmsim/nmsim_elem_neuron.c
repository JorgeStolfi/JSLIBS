/* See {nmsim_elem_neuron.h} */
/* Last edited on 2020-12-11 02:52:28 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
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
    fprintf(wr, "%s%d %d %d", ind1, ine, neu->ing, neu->nse_out);
    /* Do not write {neu.ise_start}. */
  }

nmsim_elem_neuron_t nmsim_elem_neuron_read
  ( FILE *rd, 
    nmsim_elem_neuron_ix_t ine, 
    nmsim_group_neuron_ix_t ing_max, 
    nmsim_elem_synapse_count_t nse_out_max 
  )
  { 
    (void)nmsim_read_int64_value(rd, "neuron index", ine, ine);
    nmsim_group_neuron_ix_t ing = (nmsim_group_neuron_ix_t)
      nmsim_read_int64_value(rd, "neuron group index", 0, ing_max);
    nmsim_elem_synapse_count_t nse_out = (nmsim_elem_synapse_ix_t)
      nmsim_read_int64_value(rd, "number of output synapses", 0, nse_out_max);

    nmsim_elem_neuron_t neu = 
      (nmsim_elem_neuron_t){ .ing = ing, .nse_out = nse_out, .ise_out_start = 0 };

    return neu;
  }

void nmsim_elem_neuron_show(FILE *wr, char *pref, nmsim_elem_neuron_t *neu, char *suff)
  { 
    if (pref != NULL) { fputs(pref, wr); }
    fprintf
      ( wr, "{ group = %d nse_out = %d ise_start = %d }", 
        neu->ing, neu->nse_out, neu->ise_out_start 
      );
    if (suff != NULL) { fputs(suff, wr); }
  }
  
nmsim_elem_neuron_t nmsim_elem_neuron_throw(nmsim_group_neuron_ix_t ing_max)
  {
    nmsim_group_neuron_ix_t ing = (nmsim_group_neuron_ix_t)int64_abrandom(0, ing_max);
    nmsim_elem_neuron_t neu = 
      (nmsim_elem_neuron_t){ .ing = ing, .nse_out = 0, .ise_out_start = 0 };
    return neu;
  }

void nmsim_elem_neuron_compare
  ( nmsim_elem_neuron_t *neu_read, 
    nmsim_elem_neuron_t *neu_orig
  )
  { 
    nmsim_compare_int64_param("ing",     neu_read->ing,     neu_orig->ing);
    nmsim_compare_int64_param("nse_out", neu_read->nse_out, neu_orig->nse_out);
    /* Ignores {ise_start}. */
  }

