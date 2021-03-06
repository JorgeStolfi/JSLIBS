/* See dnae_seq_ps.h */
/* Last edited on 2014-06-10 10:44:10 by stolfilocal */

#define dnae_seq_C_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include <vec.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <jsfile.h>
#include <pswr.h>

#include <msm_ps_tools.h>

#include <dnae_seq.h>
#include <dnae_nucleic.h>
#include <dnae_sample.h>
#include <dnae_datum.h>

#include <dnae_seq_ps.h>

/* IMPLEMENTATIONS */

void dnae_seq_ps_plot_named
  ( dnae_seq_t *seq, 
    double hSize, 
    double vSize,
    double fontSize,
    char *name,
    char *tag
  )
  { /* Create the {msm_ps_tools_t} object stream from it: */
    msm_ps_tools_t *dp = msm_ps_tools_new_graph
      ( NULL, name, tag, 
        /*hGraphSize:*/ hSize, 
        /*vGraphSize:*/ vSize, 
        /*scaleL:*/  TRUE,  /*titleL:*/ FALSE,
        /*scaleR:*/  FALSE, /*titleR:*/ FALSE,
        /*scaleB:*/  TRUE,  /*titleB:*/ FALSE,
        /*scaleT:*/  FALSE, /*titleT:*/ FALSE,
        /*fontSize:*/ fontSize,
        /*maxXLabChars:*/ 7,
        /*maxYLabChars:*/ 5,
        /*mrg:*/ msm_EPS_MARGIN_MM
      );
    
    /* Plot the sequence: */
    dnae_seq_ps_plot(dp, seq);
    /* Close and cleanup: */
    msm_ps_tools_close(dp);
  }
  
void dnae_seq_ps_plot(msm_ps_tools_t *dp, dnae_seq_t *seq)
  { 
    /* Get sequence size: */
    int nd = dnae_seq_num_datums(seq); /* Number of datums. */
    int nc = dnae_CHANNELS; /* Number of channels per datum. */
    if (nd == 0) { return; }

    /* Get scaling factors: */
    double *sfac = seq->sfac.f; /* Sample scaling factor per channel. */
    
    /* Choose the nominal Y range {[yMin_yMax]}: */
    int c;
    double sfMax = 1.0e-100; /* Max sample scale factor among channels. */
    for (c = 0; c < nc; c++) { if (sfac[c] > sfMax) { sfMax = sfac[c]; } }
    double vMax = dnae_sample_decode(dnae_sample_enc_VALID_MAX, sfMax);
    double ySkosh = 0.0;
    double yMin = -vMax - ySkosh;
    double yMax = +vMax + ySkosh;
    
    /* Extract the data to plot: */
    int ny = nd*nc;
    double *y = (double*) malloc(sizeof(double)*ny);
    int i;
    for (i = 0; i < nd; i++)
      { for (c = 0; c < nc; c++)
          { dnae_sample_enc_t s = dnae_seq_get_sample_enc(seq, i, c);
            y[c*nd + i] = dnae_sample_decode(s, sfac[c]);
          }
      }
       
    /* Plot the graphs: */
    int start = (int)(seq->sd.skip);
    int step = (1 << seq->sd.estep);
    msm_ps_tools_draw_graphs(dp, nc, nd,  NULL, start, step, y, yMin, yMax); 
    free(y);
  }

