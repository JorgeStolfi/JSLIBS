/* See dnae_spectrum.h */
/* Last edited on 2014-08-26 16:14:58 by stolfilocal */

#define dnae_spectrum_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#include <fftw3.h>

#include <pswr.h>
#include <jsrandom.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <vec.h>

#include <msm_image.h>
#include <msm_ps_tools.h>

#include <dnae_seq.h>
#include <dnae_datum.h>
#include <dnae_sample.h>
#include <dnae_nucleic.h>

#include <dnae_spectrum.h>

/* IMPLEMENTATIONS */

double_vec_t dnae_spectrum_from_seq(dnae_seq_t *seqp)
  {
    /* Get sequence size: */
    int nc = dnae_CHANNELS;
    int nd = dnae_seq_num_datums(seqp); /* Number of sample datums. */

    /* Get scaling factors: */
    double *sfac = seqp->sfac.f; /* Sample scaling factor per channel. */
    
    /* Allocate vectors {in,ot} for {fftw3.h} routines: */
    fftw_complex *in = fftw_malloc(sizeof(fftw_complex)*nd); /* Input vector. */
    fftw_complex *ot = fftw_malloc(sizeof(fftw_complex)*nd); /* Output vector. */

    /* Precompute parameters for {fftw3.h} routines: */
    fftw_plan plan = fftw_plan_dft_1d(nd, in, ot, FFTW_FORWARD, FFTW_ESTIMATE);

    /* Compute the total power spectrum {pwr[0..fMax]}: */
    int fMax = nd/2; /* Maximum frequency in power spectrum. */
    int nf = fMax+1; /* Number of enries in power spectrum. */

    /* Allocate and clear the power spectrum: */
    double_vec_t pwr = double_vec_new(nf);
    int f;
    for (f = 0; f <= fMax; f++) { pwr.e[f] = 0; }
    
    int c;
    for (c = 0; c < nc; c++)
      { /* Copy channel {c} samples, scaled, to {in} vector, as complex numbers: */
        int i;
        for (i = 0; i < nd; i++)
          { dnae_sample_enc_t s = dnae_seq_get_sample_enc(seqp, i, c);
            in[i][0] = dnae_sample_decode(s, sfac[c]);
            in[i][1] = 0; 
          }

        /* Compute Fourier transform: */
        fftw_execute(plan);

        /* Compute the power spectrum and add it to {pwr[0..nf-1]}: */
        double norm = 1/((double)nd); /* Power normalization factor. */
        int f;
        for (f = 0; f <= fMax; f++)
          { /* Add square modulus of Fourier terms {f} and {-f}, if distinct: */
            double P = ot[f][0]*ot[f][0] + ot[f][1]*ot[f][1];
            int g = (nd - f) % nd;
            if (g != f) { P += ot[g][0]*ot[g][0] + ot[g][1]*ot[g][1]; }
            pwr.e[f] += norm*P;
          }
      }
    return pwr;
  }

void dnae_spectrum_postscript_plot(msm_ps_tools_t *dp, dnae_seq_t *seqp)
  { 
    /* Get sequence size: */
    int nd = dnae_seq_num_datums(seqp); /* Number of sample datums. */
    if (nd == 0) { return; }

    /* Get scaling factors: */
    double *sfac = seqp->sfac.f; /* Sample scaling factor per channel. */
    
    /* Estimate the maximum power {pwrMax} per specral frequency: */
    int c;
    double datPwrMax = 0.0;      /* Max power per datum. */
    for (c = 0; c < dnae_CHANNELS; c++)
      { double smpMax = dnae_sample_decode(dnae_sample_enc_VALID_MAX, sfac[c]);  /* Max value of decoded sample. */
        datPwrMax += smpMax*smpMax;
      }
    /* Assume that power spectrum is more or less uniform: */
    double pwrMax = 2*datPwrMax;

    /* Choose the nominal Y range {[yMin_yMax]}: */
    double yMin = 0.0;
    double yMax = 1.50 * pwrMax;
    
    /* Compute the sequence's power spectrum: */
    double_vec_t pwr = dnae_spectrum_from_seq(seqp);
      
    /* Plot the graphs: */
    msm_ps_tools_draw_histogram(dp, pwr.ne, NULL, pwr.e, yMin, yMax); 
    
    free(pwr.e);
  }

void dnae_spectrum_postscript_plot_named
  ( dnae_seq_t *seqp, 
    double hSize, 
    double vSize,
    double fontSize,
    char *name,
    char *tag
  )
  { /* Create the {msm_ps_tools_t} object stream from it: */
    msm_ps_tools_t *dp = msm_ps_tools_new_graph
      ( NULL, name, tag,
        hSize, vSize, 
        /*scaleL*/  TRUE,  /*titleL*/ FALSE,
        /*scaleR*/  FALSE, /*titleR*/ FALSE,
        /*scaleB*/  TRUE,  /*titleB*/ FALSE,
        /*scaleT*/  FALSE, /*titleT*/ FALSE,
        fontSize,
        /*maxXLabChars*/ 7,
        /*maxYLabChars*/ 5,
        msm_EPS_MARGIN_MM
      );
    
    /* Plot the sequence: */
    dnae_spectrum_postscript_plot(dp, seqp);
    /* Close and cleanup: */
    msm_ps_tools_close(dp);
  }
