#ifndef dnae_spectrum_H
#define dnae_spectrum_H

/* Filtered DNA sequences */
/* Last edited on 2014-08-05 22:46:29 by stolfilocal */

#define dnae_spectrum_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#include <vec.h>
#include <pswr.h>

#include <msm_ps_tools.h>

#include <dnae_sample.h>
#include <dnae_datum.h>
#include <dnae_seq.h>

/* This interface provides tools to compute and plot
  the power spectrum of a sampled multi-channel signal
  (a {dm-seq_t}). */

double_vec_t dnae_spectrum_from_seq(dnae_seq_t *seqp);
  /* Computes the power spectrum {pwr[0..fMax]} of sequence {seqp},
    where {fMax = nd/2} is the max absolute frequency, and {nd} is the
    number of datums in the sequence. The the three channels are
    decoded and scaled by {seqp->sfac}, and their power spectra are
    added together. */

/* POSTSCRIPT PLOTTING */

void dnae_spectrum_postscript_plot(msm_ps_tools_t *dp, dnae_seq_t *seqp);
  /* Generates a Postscript plot of of the power spectrum of sequence
    {seqp}, written out to the Postscript plotting stream {dp}. The
    power spectra of all three channels are added together. The plot
    uses the whole plotting area of {dp} (with a bit of slack). */ 

void dnae_spectrum_postscript_plot_named
  ( dnae_seq_t *seqp, 
    double hSize, 
    double vSize,
    double fontSize,
    char *name,
    char *tag
  );
  /* Generates a Postscript plot of the power spectrum of sequence
    {seqp}, written out to the Encapsulated Postscript file called
    "{name}{tag}.eps". The power spectra of all three channels are added
    together. The figure will have {hSize} by {vSize} millimeters,
    including a blank margin {msm_EPS_MARGIN_MM} millimeters wide all
    around. The labels will use the specified {fontSize} (in pt). The
    plot uses the whole plotting area of {dp} (with a bit of slack). */

#endif
