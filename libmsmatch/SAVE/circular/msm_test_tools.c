/* See {msm_test_tools.h} */
/* Last edited on 2008-04-19 16:23:38 by stolfi */

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
    msm_dyn_tableau_t *tb,
    int maxIter
  ) 
  { /* Get some parameters: */
    int nx = xp->npos;
    int ny = yp->npos;
    int nmax = (nx > ny ? nx : ny);
    /* Build a Bresenham pairing {gv[0..ng-1]} betwen {xp} and {yp}: */
    int ng = 0; 
    msm_rung_vec_t gv = msm_rung_vec_new(nmax);
    /* Store the initial rung {gini}: */
    msm_rung_t gini = (msm_rung_t){{ 0, 0 }};
    msm_rung_vec_expand(&gv,ng); gv.e[ng] = gini; ng++;
    /* Interpolate between {gini} and {gfin}: */
    msm_rung_t gfin = (msm_rung_t){{ nx-1, ny-1 }};
    msm_rung_interpolate(gini, gfin, &ng, &gv);
    msm_rung_vec_trim(&gv, ng);
    /* Make it into a candidate: */
    msm_pairing_t *pr = msm_pairing_from_rung_vec(&gv, FALSE);
    msm_cand_t cd = msm_cand_from_pairing(xp, yp, pr, /*score*/ 0.0);
    /* Refine it: */
    int iter = 0;
    while (iter < maxIter)
      { msm_cand_t cdref = msm_cand_refine(&cd, delta, 0, FALSE, FALSE, 6, step_score, tb);
        if (msm_cand_equivalent(&cd, &cdref, FALSE)) 
          { msm_pairing_free(cdref.pr); break; }
        else
          { msm_pairing_free(cd.pr); cd = cdref; }
      }
    return cd; 
  }

void msm_test_seq_write_and_plot_named
  ( msm_seq_desc_t *seq,
    int den,
    double_vec_t *smp, 
    char *title,
    char *name,
    char *tag,
    double fontSize
  )
  {
    fprintf(stderr, "sequence = "); 
    msm_seq_desc_write(stderr, " ( ", seq, 1, 1, 4, " )"); 
    fprintf(stderr, "  samples = %4d  den = %2d\n", smp->ne, den);

    /* Check whether the sample vector is compatible with the seq descriptor: */
    int ns = smp->ne;
    int npos = (seq->circ ? den*ns : den*(ns-1)+1);
    demand(seq->npos == npos, "inconsistent sequence descriptor");
    
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
        /*maxLabChars:*/ 4,
        /*mrg:*/ 1.0
      );
    msm_ps_tools_draw_graphs
      ( mps, 1, smp->ne, seq->circ, NULL, smp->e, +INF, -INF );
  }

