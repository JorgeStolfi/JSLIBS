/* See {msm_test_tools.h} */
/* Last edited on 2022-10-20 10:27:56 by stolfi */

#define msm_test_tools_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include <affirm.h>
#include <bool.h>

#include <msm_basic.h>
#include <msm_seq_desc.h>
#include <msm_double_vec.h>
#include <msm_cand_refine.h>
#include <msm_ps_tools.h>

#include <msm_test_tools.h>

void msm_debug_double_vec(double *v, int32_t nv, char *fmt)
  { fprintf(stderr, "[");
    int32_t j;
    for (j = 0; j < nv; j++)
      { fprintf(stderr, " "); fprintf(stderr, fmt, v[j]); }
    fprintf(stderr, " ]");
  }
    
void msm_debug_int32_vec(int32_t *v, int32_t nv, char *fmt)
  { fprintf(stderr, "[");
    int32_t j;
    for (j = 0; j < nv; j++)
      { fprintf(stderr, " "); fprintf(stderr, fmt, v[j]); }
    fprintf(stderr, " ]");
  }

msm_cand_t msm_test_tools_get_optimum_pairing
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1, 
    msm_rung_t gini,
    msm_rung_t gfin, 
    int32_t delta, 
    msm_rung_step_score_proc_t *step_score, 
    bool_t verbose,
    msm_dyn_tableau_t *tb,
    int32_t maxIter
  ) 
  { /* Get some parameters: */
    int32_t n0 = seq0->size;
    int32_t n1 = seq1->size;
    int32_t nmax = (n0 > n1 ? n0 : n1);
    /* Build a Bresenham pairing {gv[0..ng-1]} betwen {seq0} and {seq1}: */
    int32_t ng = 0; 
    msm_rung_vec_t gv = msm_rung_vec_new(nmax);
    /* Store the initial rung {gini}: */
    msm_rung_vec_expand(&gv,ng); gv.e[ng] = gini; ng++;
    /* Interpolate between {gini} and {gfin}: */
    msm_rung_interpolate(gini, gfin, &ng, &gv);
    msm_rung_vec_trim(&gv, ng);
    /* Make it into a candidate: */
    msm_pairing_t *pr = msm_pairing_from_rung_vec(&gv);
    msm_cand_t cd = msm_cand_from_pairing(seq0, seq1, pr, /*score*/ 0.0);
    /* Refine it: */
    int32_t iter = 0;
    while (iter < maxIter)
      { int32_t n_steps = 0;
        int32_t n_entries = 0;
        msm_cand_t cdref = msm_cand_refine(&cd, delta, 0, FALSE, FALSE, 6, step_score, verbose, tb, &n_steps, &n_entries);
        if (msm_cand_equivalent(&cd, &cdref, FALSE)) 
          { msm_pairing_free(cdref.pr); break; }
        else
          { msm_pairing_free(cd.pr); cd = cdref; }
      }
    return cd; 
  }

void msm_test_seq_write_and_plot_named
  ( msm_seq_desc_t *seq,
    double_vec_t *smp, 
    char *title,
    char *name,
    char *tag,
    double fontSize,
    bool_t orig
  )
  {
    fprintf(stderr, "sequence = "); 
    msm_seq_desc_write(stderr, " ( ", seq, 1, 1, 4, " )"); 
    fprintf(stderr, "  samples = %4d\n", smp->ne);

    /* Estimate size of scale labels: */
    int32_t maxXLabChars = (int32_t)ceil(log(fmax(smp->ne,2)*1.05)/M_LN10);
    int32_t maxYLabChars = 5; /* E.g. "+1.00" */
    
    /* Ignore the {seq} for now: */
    msm_double_vec_write_named(smp, name, tag);
    
    /* Plotting: */
    msm_ps_tools_t *mps = msm_ps_tools_new_graph
      ( NULL, name, tag, 
        /*hGraphSize:*/ 165.0, 
        /*vGraphSize:*/ 25.0, 
        /*scaleL:*/ FALSE, /*titleL:*/ FALSE, 
        /*scaleR:*/ FALSE, /*titleR:*/ FALSE,
        /*scaleB:*/ FALSE, /*titleB:*/ FALSE,
        /*scaleT:*/ FALSE, /*titleT:*/ FALSE, 
        /*fontSize:*/ fontSize,
        /*maxXLabChars:*/ maxXLabChars,
        /*maxYLabChars:*/ maxYLabChars,
        /*mrg:*/ 1.0
      );
    double skip = (orig ? seq->skip : 0.0);
    double step = (orig ? pow(2.0, seq->estep) : 1.0);
    msm_ps_tools_draw_graphs
      ( mps, 1, smp->ne, NULL, skip, step, smp->e, +INF, -INF );
  }

