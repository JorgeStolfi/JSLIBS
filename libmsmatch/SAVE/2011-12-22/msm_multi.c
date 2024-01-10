/* See msm_multi.h */
/* Last edited on 2008-04-19 16:23:15 by stolfi */ 

#define msm_multi_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>

#include <msm_basic.h>
#include <msm_seq_desc.h>
#include <msm_cand.h>
#include <msm_cand_vec.h>
#include <msm_multi.h>
#include <msm_cand_refine.h>

#include <affirm.h>

msm_cand_vec_t msm_multi_find_matches
  ( msm_cand_vec_t *cdvini,
    msm_seq_desc_t a[], 
    msm_seq_desc_t b[],
    int maxMatches,
    int maxLevel,
    int minCover,
    int nwtb0,
    int nwtb1,
    double frac,
    bool_t refine,
    int delta,
    int kappa,
    bool_t expand,
    bool_t shrink,
    int maxUnp,
    msm_rung_step_score_proc_t *step_score,
    msm_multi_report_proc_t *report
  )
  { /* Working matrices: */
    msm_dyn_tableau_t tb = msm_dyn_tableau_new(); 

    /* Current candidate set: */
    msm_cand_vec_t cdv = msm_cand_vec_new(0);
    
    /* Multiscale search and refine at all scales: */
    int level;
    for (level = maxLevel; level >= 0; level--)
      { 
        /* Scale reduction factor from original strings: */
        int scale = 1 << level;

        fprintf(stderr, "\n");
        fprintf(stderr, "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n");
        fprintf(stderr, "LEVEL %d - SCALE %d\n", level, scale);

        /* Strings that are being matched at this level: */
        msm_seq_desc_t *ap = &(a[level]);
        msm_seq_desc_t *bp = &(b[level]);

        demand(ap->level == level, "inconsistent {level} of {ap}");
        demand(bp->level == level, "inconsistent {level} of {bp}");

        /* Get the unrefined candidate list {cdvraw} for this level: */
        msm_cand_vec_t cdvraw = msm_cand_vec_new(0); 
        if (level == maxLevel)
          { /* Let {cdvraw} be an alias of {cdvini}: */
            cdvraw = (*cdvini);
          }
        else
          { /* Let {cdvraw} be {*cdv} mapped to the current level: */
            fprintf(stderr, "mapping %d candidates", cdv.ne);
            fprintf(stderr, " from level %d to level %d ...\n", level+1, level);
            /* Map the candidates: */
            int nwtb = (level == 0 ? nwtb0 : nwtb1); /* Width of filter table. */
            cdvraw = msm_cand_vec_map_to_finer
              ( &cdv, ap->npos, bp->npos, ap->circ, bp->circ, nwtb );
          }

        /* Reclaim the previous {cdv} candidate list, if not shared: */
        if ((cdv.e != cdvraw.e) && (cdv.e != cdvini->e)) { free(cdv.e); }

        /* Report the unrefined candidates: */
        if (report != NULL) { report(ap, bp, &cdvraw, FALSE); }

        /* Convert {cdvraw} to the definitive list, save in {cdv}: */
        if (refine)
          { /* Refine, sort, and prune {cdvraw}: */
            fprintf(stderr, "refining %d candidates\n", cdvraw.ne);
            int nprune = maxMatches * scale;
            fprintf(stderr, "  pruning to at most %d candidates\n", nprune);
            int mincov = (minCover + scale - 1)/scale;
            fprintf(stderr, "  min sequence coverage %d\n", mincov);
            cdv = msm_cand_vec_refine
              ( &cdvraw, ap, bp, 
                delta, kappa, expand, shrink, maxUnp,
                step_score, &tb, 
                mincov, nprune, frac
              );
          }
        else
          { /* Let {cdv} be an alias for {cdvraw}: */
            cdv = cdvraw;
          }

        /* Reclaim the {cdvraw} candidate list, if not shared: */
        if ((cdvraw.e != cdv.e) && (cdvraw.e != cdvini->e)) { free(cdvraw.e); }

        /* Report the refined candidates: */
        if (report != NULL) { report(ap, bp, &cdv, TRUE); }

        fprintf(stderr, "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n");
        fprintf(stderr, "\n");

      }

    fprintf(stderr, "mapped and refined %d candidates at finest level\n", cdv.ne);
    return cdv;
  }
