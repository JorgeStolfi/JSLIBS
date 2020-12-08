/* See msm_multi.h */
/* Last edited on 2017-04-28 10:45:57 by stolfilocal */ 

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
    msm_seq_desc_t sv0[], 
    msm_seq_desc_t sv1[],
    int maxMatches,
    int maxLevel,
    int minSamples[],
    double frac,
    bool_t refine,
    int delta,
    int kappa,
    int expand,
    int shrink,
    int maxUnp,
    msm_rung_step_score_proc_t *step_score,
    msm_multi_report_proc_t *report,
    bool_t verbose
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
        msm_seq_desc_t *seq0 = &(sv0[level]);
        msm_seq_desc_t *seq1 = &(sv1[level]);

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
            cdvraw = msm_cand_vec_map(&cdv, seq0, seq1);
          }

        /* Reclaim the previous {cdv} candidate list, if not shared: */
        if ((cdv.e != cdvraw.e) && (cdv.e != cdvini->e)) { free(cdv.e); }

        /* Report the unrefined candidates: */
        if (report != NULL) { report(level, seq0, seq1, &cdvraw, FALSE); }

        /* Convert {cdvraw} to the definitive list, save in {cdv}: */
        if (refine)
          { /* Refine, sort, and prune {cdvraw}: */
            fprintf(stderr, "refining %d candidates\n", cdvraw.ne);
            int nprune = maxMatches * scale;
            fprintf(stderr, "  pruning to at most %d candidates\n", nprune);
            int mincov = minSamples[level];
            fprintf(stderr, "  min covered samples in each sequence %d\n", mincov);
            assert(mincov > 0);
            cdv = msm_cand_vec_refine 
              ( &cdvraw, seq0, seq1, 
                delta, kappa, expand, shrink, maxUnp,
                step_score, verbose, &tb, 
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
        if (report != NULL) { report(level, seq0, seq1, &cdv, TRUE); }
        
        fprintf(stderr, "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n");
        fprintf(stderr, "\n");

      }
    
    msm_dyn_tableau_free(&tb);

    fprintf(stderr, "mapped and refined %d candidates at finest level\n", cdv.ne);
    return cdv;
  }
