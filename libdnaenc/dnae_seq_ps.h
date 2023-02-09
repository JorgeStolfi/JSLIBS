#ifndef dnae_seq_ps_H
#define dnae_seq_ps_H

/* Postsscript plotting of numerically encoded DNA sequences */
/* Last edited on 2023-02-07 20:36:53 by stolfi */

#define dnae_seq_ps_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#define _GNU_SOURCE
#include <stdio.h>

#include <vec.h>

#include <msm_ps_tools.h>

#include <dnae_sample.h>
#include <dnae_datum.h>
#include <dnae_seq.h>

/* This interface defines a discrete representation for 
  signals, as strings of /datum/ elements ({dnae_datum_t}s). */

/* POSTSCRIPT PLOTTING */

void dnae_seq_ps_plot(msm_ps_tools_t *dp, dnae_seq_t *seq);
  /* Generates a Postscript plot of the sequence {seq}, written out
    to the Postscript plotting stream {dp}.  Each channel is decoded and
    drawn with a different color.  The plot uses the whole
    plotting area of {dp} (with a bit of slack). */ 

void dnae_seq_ps_plot_named
  ( dnae_seq_t *seq, 
    double hSize, 
    double vSize,
    double fontSize,
    char *name,
    char *tag
  );
  /* Generates a Postscript plot of the sequence {seq}, written out
    to the Encapsulated Postscript file called "{name}{tag}.eps". Each
    channel is decoded and drawn with a different color. The figure
    will have {hSize} by {vSize} millimeters, including a blank margin
    {msm_EPS_MARGIN_MM} millimeters wide all around. The spectrum plot
    is scaled to fit in the stated area.  */

#endif
