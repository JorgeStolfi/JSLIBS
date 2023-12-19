#ifndef neuromat_filter_bandpass_log_gauss_H
#define neuromat_filter_bandpass_log_gauss_H

/* NeuroMat basic bandpass filter transfer function with log-Gaussian profile. */
/* Last edited on 2023-12-15 08:04:42 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

double neuromat_filter_bandpass_log_gauss_compute_sigma(double fm, double fc, double gc);
  /* Computes the deviation {sigma} of a log-Gaussian function that
    has unit height, average {log(fm)}, and value {gc} at frequency {fc}.
    Namely, {sigma} is such that exp(-0.5*((log(fc) - log(fm))/sigma)^2) = gc}.
    
    The frequencies {fc} and {fm} must be positive and distinct, and
    {gc} must be strictly between 0 and 1. The sign of the result will
    be negative iff {fc} is less than {fm}. */

double neuromat_filter_bandpass_log_gauss(double f, double fm, int32_t np, double sigma);
  /* A filter whose profile is the sum of {np} Gaussians in log freq space 
    Each Gaussian has unit height, deviation {sigma}, and mean {log(fm) + 2*k*sigma}
    for {k} in {0..np-1}.  Note that the sign of {sigma} affects the position
    of the successive pulses. */

#endif
