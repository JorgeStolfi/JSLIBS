/* See {dnae_test_tools.h} */
/* Last edited on 2022-10-31 11:21:54 by stolfi */

#define dnae_test_tools_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <jsrandom.h>

#include <affirm.h>
#include <bool.h>

#include <msm_basic.h>
#include <msm_seq_desc.h>
#include <msm_double_vec.h>
#include <msm_ps_tools.h>
#include <msm_test_tools.h>

#include <dnae_nucleic.h>
#include <dnae_seq.h>

#include <dnae_test_tools.h>

void dnae_test_tools_seq_write_and_plot_named
  ( dnae_seq_t *seq,
    char *title,
    char *name,
    char *tag,
    bool_t plot,
    double hSize,
    double vSize,
    double fontSize,
    int32_t maxXLabChars
  )
  {
    bool_t debug = TRUE;
    
    fprintf(stderr, "sequence = "); 
    msm_seq_desc_write(stderr, " ( ", &(seq->sd), 1, 1, 4, " )"); 
    fprintf(stderr, "  samples = %d\n", seq->dv.ne);

    /* Check whether the sample vector is compatible with the seq descriptor: */
    
    int32_t size = dnae_seq_num_datums(seq);
    assert(size == seq->dv.ne);
    
    /* Ignore the {seq} for now: */
    dnae_seq_write_named(name, tag, ".eqs", seq);
    
    if (plot)
      { /* Plot the seqence: */
        if (maxXLabChars < 0) { maxXLabChars = (int32_t)ceil(log(fmax(size,9)*1.05)/M_LN10); }
        int32_t maxYLabChars = 5; 

        msm_ps_tools_t *mps = msm_ps_tools_new_graph
          ( NULL, name, tag, 
            /*hGraphSize:*/ hSize, 
            /*vGraphSize:*/ vSize, 
            /*scaleL:*/ FALSE, /*titleL:*/ FALSE, 
            /*scaleR:*/ FALSE, /*titleR:*/ FALSE,
            /*scaleB:*/ FALSE, /*titleB:*/ FALSE,
            /*scaleT:*/ FALSE, /*titleT:*/ FALSE, 
            /*fontSize:*/ fontSize,
            /*maxXLabChars:*/ maxXLabChars,
            /*maxYLabChars:*/ maxYLabChars,
            /*mrg:*/ 1.0
          );

        /* Repack the samples as a single vector {smp[]}, for {msm_ps_tools_draw_graphs}: */
        int32_t nplt = size; /* Number of samples for plotting. */
        double smp[dnae_CHANNELS*nplt]; 
        int32_t i, c;
        for (i = 0; i < nplt; i++)
          { if (debug) { fprintf(stderr, "datum[%6d]", i); }
            for (c = 0; c < dnae_CHANNELS; c++)
              { dnae_sample_enc_t s_enc = dnae_seq_get_sample_enc(seq, i, c);
                double s_fac = seq->sfac.f[c];
                double s_dec = dnae_sample_decode(s_enc, s_fac);
                if (debug) { fprintf(stderr, " %+6d = %+11.7f", s_enc, s_dec); }
                smp[nplt*c + i] = s_dec;
              }
            if (debug) { fprintf(stderr, "\n"); }
          }

        /* Indices in original sequence: */
        double start = seq->sd.skip;  /* Actually always an integer. */
        double step = pow(2.0, seq->sd.estep); /* may be fracional if sequence was interpolated. */
        /* Assumes that the decoded samples do not exceed 1 in magnitude: */
        double yMin = -1.2;
        double yMax = +1.2;
        msm_ps_tools_draw_graphs
          ( mps, dnae_CHANNELS, size, 
            NULL, start, step,
            smp, yMin, yMax
          );
        msm_ps_tools_close(mps);
      }
  }

void dnae_test_tools_seq_multi_write_and_plot_named
  ( dnae_seq_t seq[],
    int32_t maxLevel,
    char *title,
    char *name,
    char *tag,
    bool_t plot,
    double hSize,
    double vSize,
    double fontSize,
    int32_t maxXLabChars
  )
  { 
    /* Compute max label chars for all plots, if necessary: */
    if (plot & (maxXLabChars < 0))
      { int32_t nposMax = dnae_seq_num_datums(&(seq[0]));
        maxXLabChars = (int32_t)ceil(log(fmax(nposMax,2)*1.05)/M_LN10);
      }
    int32_t level;
    for (level = 0; level <= maxLevel; level++)
      { char *titlei = NULL; 
        asprintf(&titlei, "%s level %02d", title, level);
        char *tagi = NULL; 
        asprintf(&tagi, "%s-%02d", tag, level);
        dnae_test_tools_seq_write_and_plot_named
          ( &(seq[level]), titlei, name, tagi, 
            plot, hSize, vSize, fontSize, maxXLabChars
          );
        free(titlei); free(tagi);
      }
  }

#define dnae_rung_HUGE (msm_rung_t){{ INT32_MAX, INT32_MAX }}

void dnae_test_tools_make_seq_pair
  ( char *borg, 
    double mutProb, 
    double delProb,
    dnae_seq_id_t xid, 
    char *xtag,
    dnae_seq_t *xP,
    dnae_seq_id_t yid, 
    char *ytag,
    dnae_seq_t *yP,
    msm_pairing_t **prP,
    char *outName
  )
  { /* Generate two mutated copies of {b}, save the indices of copied chars: */
    char *xdna, *ydna; /* The mutated copies. */
    msm_rung_vec_t xgv, ygv;  /* Index pairs for letters that were copied. */
    dnae_nucleic_string_mutate(borg, TRUE, TRUE, mutProb, delProb, &xdna, &xgv);
    dnae_nucleic_string_mutate(borg, TRUE, TRUE, mutProb, delProb, &ydna, &ygv);
    /* Merge the rung vectors {xgv,ygv} by component 0 to create the pairing {*prP}: */
    int32_t minEach = 1; /* Each side of every step must increase. */
    int32_t minSum = 2;  /* The sum must increase by 2. */
    int32_t ngmax = (xgv.ne < ygv.ne ? xgv.ne : ygv.ne);
    msm_rung_vec_t gv = msm_rung_vec_new(ngmax);
    int32_t ng = 0; /* Output rungs will be {gv[0..ng-1]}. */
    int32_t kx = 0, ky = 0;
    while ((kx < xgv.ne) && (ky < ygv.ne))
      { msm_rung_t xg = (kx >= xgv.ne ? dnae_rung_HUGE : xgv.e[kx]);
        msm_rung_t yg = (ky >= ygv.ne ? dnae_rung_HUGE : ygv.e[ky]);
        if (xg.c[0] < yg.c[0]) 
          { kx++; }
        else if (xg.c[0] > yg.c[0])
          { ky++; }
        else
          { msm_rung_t g = (msm_rung_t){{ xg.c[1], yg.c[1] }};
            if (ng > 0) { (void)msm_rung_step_is_increasing(gv.e[ng-1], g, minEach, minSum, /*die*/ TRUE); }
            msm_rung_vec_expand(&gv, ng);
            gv.e[ng] = g; ng++; 
            kx++; ky++; 
          }
      }
    demand(ng > 0, "empty output pairing");
    msm_rung_vec_trim(&gv, ng);
    /* Generate pairing from rungs: */
    msm_pairing_t *pr = msm_pairing_from_rung_vec(&gv);
    /* Generate comment: */
    char *cmt = NULL;
    asprintf(&cmt, "mutated mp = %8.6f dp = %8.6f\n", mutProb, delProb);
    (*xP) = dnae_seq_from_nucleic_string(xid, xtag, FALSE, cmt, xdna);
    (*yP) = dnae_seq_from_nucleic_string(yid, ytag, FALSE, cmt, ydna);
    (*prP) = pr;
    /* Write sequences to disk: */
    dnae_test_tools_write_generated_sequence(outName, xdna, xP);
    dnae_test_tools_write_generated_sequence(outName, ydna, yP);
    free(xdna); free(xgv.e);
    free(ydna); free(ygv.e);
  }

void dnae_test_tools_write_generated_sequence(char *outName, char *b, dnae_seq_t *s)
  { /* Write sequence to disk: */
    char *fileTag = NULL;
    asprintf(&fileTag, "-%s", s->sd.name);
    dnae_nucleic_string_write_named(outName, fileTag, ".bas", b, s->cmt);
    dnae_seq_write_named(outName, fileTag, ".egs", s);
    free(fileTag);
  }


