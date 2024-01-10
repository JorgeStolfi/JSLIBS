/* See {msm_test_tools.h} */
/* Last edited on 2008-01-12 12:12:51 by stolfi */

#define msm_test_tools_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <string.h>
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

void msm_debug_double_vec(double *v, int nv, char *fmt)
  { fprintf(stderr, "[");
    int j;
    for (j = 0; j < nv; j++)
      { fprintf(stderr, " "); fprintf(stderr, fmt, v[j]); }
    fprintf(stderr, " ]");
  }
    
void msm_debug_int_vec(int *v, int nv, char *fmt)
  { fprintf(stderr, "[");
    int j;
    for (j = 0; j < nv; j++)
      { fprintf(stderr, " "); fprintf(stderr, fmt, v[j]); }
    fprintf(stderr, " ]");
  }

msm_cand_t msm_test_tools_get_optimum_pairing
  ( msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp, 
    int delta, 
    msm_rung_step_score_proc_t *step_score, 
    msm_dyn_tableau_t *tb
  ) 
  { /* Get some parameters: */
    int nx = xp->nbas;
    int ny = yp->nbas;
    int nmax = (nx > ny ? nx : ny);
    /* Build a Bresenham pairing {gv[0..ng-1]} betwen {xp} and {yp}: */
    int ng = 0; 
    msm_rung_vec_t gv = msm_rung_vec_new(nmax);
    /* Store the initial rung {gini}: */
    msm_rung_t gini = (msm_rung_t){{ 0, 0 }};
    msm_rung_vec_expand(&gv,ng); gv.el[ng] = gini; ng++;
    /* Interpolate between {gini} and {gfin}: */
    msm_rung_t gfin = (msm_rung_t){{ nx-1, ny-1 }};
    msm_rung_interpolate(gini, gfin, &ng, &gv);
    msm_rung_vec_trim(&gv, ng);
    /* Make it into a candidate: */
    msm_pairing_t *pr = msm_pairing_from_rung_vec(&gv, FALSE);
    msm_cand_t cdraw = msm_cand_from_pairing(xp, yp, pr, /*score*/ 0.0);
    /* Refine it: */
    msm_cand_t cdopt = msm_cand_refine(&cdraw, delta, 0, FALSE, 6, step_score, tb);
    return cdopt; 
  }

void msm_test_seq_write_and_plot_named
  ( msm_seq_desc_t *seq,
    int den,
    double_vec_t *smp, 
    char *title,
    char *name,
    char *tag
  )
  {
    fprintf(stderr, "sequence = "); 
    msm_seq_desc_write(stderr, " ( ", seq, 1, 1, 4, " )"); 
    fprintf(stderr, "  samples = %4d  den = %2d\n", smp->nel, den);

    /* Check whether the sample vector is compatible with the seq descriptor: */
    int ns = smp->nel;
    int nbas = (seq->circ ? den*ns : den*(ns-1)+1);
    demand(seq->nbas == nbas, "inconsistent sequence descriptor");
    
    /* Ignore the {seq} for now: */
    msm_double_vec_write_named(smp, name, tag);
    
    /* Plotting: */
    msm_ps_tools_t *mps = msm_ps_tools_new_graph
      ( NULL, name, tag, 
        /*hGraphSize:*/ 150.0, 
        /*vGraphSize:*/ 25.0, 
        /*scaleL:*/ FALSE, /*titleL:*/ FALSE, 
        /*scaleR:*/ FALSE, /*titleR:*/ FALSE,
        /*scaleB:*/ FALSE, /*titleB:*/ FALSE,
        /*scaleT:*/ FALSE, /*titleT:*/ FALSE, 
        /*fontSize:*/ 8.0,
        /*maxLabChars:*/ 4,
        /*mrg:*/ 1.0
      );
    msm_ps_tools_draw_graphs(mps, 1, smp->nel, seq->circ, NULL, smp->el, +INF, -INF);
  }

