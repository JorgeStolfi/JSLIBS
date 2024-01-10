/* See {nmsim_group_neuron_state.h} */
/* Last edited on 2019-07-02 15:14:53 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


#include <bool.h>
#include <affirm.h>
#include <fget.h>

#include <nmsim_basic.h>
#include <nmsim_firing_func.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_neuron_state.h>
#include <nmsim_cohort_neuron_state.h>

nmsim_group_neuron_state_t nmsim_group_neuron_state_make(void)
  { nmsim_group_neuron_state_t st = (nmsim_group_neuron_state_t)
      { .nc = 0, .cs = nmsim_cohort_neuron_state_vec_new(0) };
    return st;
  }

void nmsim_group_neuron_state_free(nmsim_group_neuron_state_t *st)
  { 
    if (st->cs.e != NULL) { free(st->cs.e); }
  }

void nmsim_group_neuron_state_throw
  ( nmsim_group_neuron_state_t *st, 
    nmsim_class_neuron_t *nclass,
    nmsim_cohort_neuron_count_t nc
  )
  {
    /* Allocate and initialize the cohort vector: */
    st->nc = nc;
    if (st->cs.e == NULL) { st->cs.ne = 0; } /* Just in case. */
    assert((st->cs.ne == 0) == (st->cs.e == NULL));
    
    /* Ensure that it has storage for {nc+1} cohorts: */
    nmsim_cohort_neuron_state_vec_expand(&(st->cs), nc);
    double V_B = nclass->V_B;
    double c_B = nclass->c_B;
    double M = nclass->M_R;
    double H = nclass->H_R; 
    double V_avg = nclass->V_R;
    double V_dev = 0.0;
    for (nmsim_cohort_neuron_ix_t k = 0; k <= nc; k++)
      { /* These cohort parameters are irrelevant for {k < nc}. */
        double eta = (k == nc ? 1.0 : 0.0);
        nmsim_cohort_neuron_state_set(&(st->cs.e[k]), V_avg, V_dev, M, H, eta);
        /* Compute potentials and modulators for next cohort: */
        double c = c_B*M;
        V_avg = V_B + c*(V_avg - V_B);
        M = 1.0 + (M - 1.0)*nclass->M_mu;
        H = 1.0 + (H - 1.0)*nclass->H_mu; 
      }
  }

void nmsim_group_neuron_state_copy
  ( nmsim_group_neuron_state_t *src, 
    nmsim_group_neuron_state_t *dst
  )
  {
    nmsim_cohort_neuron_count_t nc = src->nc;
    dst->nc = nc;
    /* Make sure that {dst->cs} has space for cohorts {0..nc}: */
    nmsim_cohort_neuron_state_vec_expand(&(dst->cs), nc);
    /* Copy the cohort states: */
    for (nmsim_cohort_neuron_ix_t k = 0; k <= nc; k++)
      { dst->cs.e[k] = src->cs.e[k]; }
    /* Copy auxiliary data: */
  }

void nmsim_group_neuron_state_write(FILE *wr, nmsim_group_neuron_state_t *st)
  {
    fprintf(wr, "  %5d", st->nc);
    fprintf(wr, "  ");
    
    for (nmsim_cohort_neuron_ix_t k = 0; k <= st->nc; k++)
      { nmsim_cohort_neuron_state_t *cs = &(st->cs.e[k]);
        nmsim_cohort_neuron_state_write(wr, k, cs);
      }
  }

void nmsim_group_neuron_state_read
  ( FILE *rd,
    nmsim_group_neuron_state_t *st
  )
  {
    fget_skip_spaces(rd);
    /* Read cohort states: */
    st->nc = (nmsim_cohort_neuron_count_t)nmsim_read_int64_value(rd, "nc", 0, nmsim_cohort_neuron_count_MAX);
    st->cs = nmsim_cohort_neuron_state_vec_new(st->nc + 1);
    for (nmsim_cohort_neuron_ix_t k = 0; k <= st->nc; k++)
      { nmsim_cohort_neuron_state_t *cs = &(st->cs.e[k]);
        nmsim_cohort_neuron_state_read(rd, k, cs);
      }
  }

void nmsim_group_neuron_state_compare(nmsim_group_neuron_state_t *st_read, nmsim_group_neuron_state_t *st_orig)
  {
    nmsim_compare_int64_param("cohort_count",  st_read->nc, st_orig->nc);
    
    for (nmsim_cohort_neuron_ix_t k = 0; k <= st_orig->nc; k++)
      { nmsim_cohort_neuron_state_t *cs_read = &(st_read->cs.e[k]);
        nmsim_cohort_neuron_state_t *cs_orig = &(st_orig->cs.e[k]);
        nmsim_cohort_neuron_state_compare(cs_read, cs_orig);
      }
  }    


