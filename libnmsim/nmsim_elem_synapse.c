/* See {nmsim_elem_synapse.h} */
/* Last edited on 2020-12-10 14:32:43 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsrandom.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_elem_synapse.h>

void nmsim_elem_synapse_write(FILE *wr, nmsim_elem_synapse_ix_t ise, nmsim_elem_synapse_t *syn)
  {
    char *ind1 = "  ";   /* Indent of synapse line. */
    fprintf(wr, "%s%d %d  %d %d  ", ind1, ise, syn->isg, syn->ine_pre, syn->ine_pos);
    nmsim_write_double_value(wr, (double)syn->W, nmsim_write_KW_PREC, TRUE, TRUE, FALSE);
  }

nmsim_elem_synapse_t nmsim_elem_synapse_read
  ( FILE *rd, 
    nmsim_elem_synapse_ix_t ise, 
    nmsim_group_synapse_ix_t isg_max, 
    nmsim_elem_neuron_ix_t ine_max
  )
  { 
    (void)nmsim_read_int64_value(rd, "synapse index", ise, ise);
    nmsim_group_synapse_ix_t isg = 
      (nmsim_group_synapse_ix_t)nmsim_read_int64_value(rd, "synapse group index", 0, isg_max);
    nmsim_elem_neuron_ix_t ine_pre = 
      (nmsim_elem_neuron_ix_t)nmsim_read_int64_value(rd, "pre-synaptic neuron index", 0, ine_max);
    nmsim_elem_neuron_ix_t ine_pos = 
      (nmsim_elem_neuron_ix_t)nmsim_read_int64_value(rd, "post-synaptic neuron index", 0, ine_max);
    float W = (float)nmsim_read_double_value(rd, "resting synaptic weight", -1000.0, +1000.0);
    
    nmsim_elem_synapse_t syn = (nmsim_elem_synapse_t)
      { .isg = isg, .ine_pre = ine_pre, .ine_pos = ine_pos, .W = W };

    return syn;
  }

void nmsim_elem_synapse_show(FILE *wr, char *pref, nmsim_elem_synapse_t *syn, char *suff)
  { 
    if (pref != NULL) { fputs(pref, wr); }
    fprintf(wr, "{ group = %d", syn->isg);
    fprintf(wr, " ine_pre = %d ine_pos = %d", syn->ine_pre, syn->ine_pos);
    fprintf(wr, " W = %.3f", syn->W);
    fprintf(wr, " }");
    if (suff != NULL) { fputs(suff, wr); }
  }
  
nmsim_elem_synapse_t nmsim_elem_synapse_throw
  ( nmsim_group_synapse_ix_t isg_min, 
    nmsim_group_synapse_ix_t isg_max, 
    nmsim_elem_neuron_ix_t ine_pre_min,
    nmsim_elem_neuron_ix_t ine_pre_max, 
    nmsim_elem_neuron_ix_t ine_pos_min,
    nmsim_elem_neuron_ix_t ine_pos_max,
    double W_avg,
    double W_dev
  )
  {
    demand((isg_min >= 0) && (isg_min <= isg_max), "invalid synapse group index range");
    demand(isg_max <= nmsim_group_synapse_ix_MAX, "invalid max synapse group index");
    demand((ine_pre_min >= 0) && (ine_pre_min <= ine_pre_max), "invalid presynaptic neuron index range");
    demand((ine_pre_max >= 0) && (ine_pre_max <= nmsim_elem_neuron_ix_MAX), "invalid max presynaptic neuron index");
    demand((ine_pos_min >= 0) && (ine_pos_min <= ine_pos_max), "invalid postsynaptic neuron index range");
    demand((ine_pos_max >= 0) && (ine_pos_max <= nmsim_elem_neuron_ix_MAX), "invalid max postsynaptic neuron index");

    nmsim_group_synapse_ix_t isg = (nmsim_group_synapse_ix_t)int64_abrandom(isg_min, isg_max);
    nmsim_elem_neuron_ix_t ine_pre = (nmsim_elem_neuron_ix_t)int64_abrandom(ine_pre_min, ine_pre_max);
    nmsim_elem_neuron_ix_t ine_pos = (nmsim_elem_neuron_ix_t)int64_abrandom(ine_pos_min, ine_pos_max);
    
    /* Pick a synaptic strength {W}: */
    float W = (float)dloggaussrand(W_avg, W_dev);
    if (W_dev != 0.0) 
      { fprintf(stderr, "  synapse weight %18.12f", W);
        double mod = nmsim_select_rounding_mod(-W_dev/100, +W_dev/100);
        W = (float)(mod*floor(W/mod + 0.5));
        fprintf(stderr, " rounding modulus = %18.12f rounded to %18.12f\n", mod, W);
      }

    nmsim_elem_synapse_t syn = (nmsim_elem_synapse_t)
      { .isg = isg, .ine_pre = ine_pre, .ine_pos = ine_pos, .W = W };
    return syn;
  }

void nmsim_elem_synapse_compare(nmsim_elem_synapse_t *syn_read, nmsim_elem_synapse_t *sin_orig)
  {
    nmsim_compare_int64_param("isg",  syn_read->isg, sin_orig->isg);
    nmsim_compare_int64_param("ine_pre",  syn_read->ine_pre, sin_orig->ine_pre);
    nmsim_compare_int64_param("ine_pos",  syn_read->ine_pos, sin_orig->ine_pos);
    nmsim_compare_double_param("W",  syn_read->W, sin_orig->W, nmsim_write_KW_PREC, TRUE, FALSE);
  }

