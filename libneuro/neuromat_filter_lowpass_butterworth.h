#ifndef neuromat_filter_lowpass_butterworth_H
#define neuromat_filter_lowpass_butterworth_H

/* NeuroMat basic Butterworth low-pass filter transfer function. */
/* Last edited on 2024-01-05 17:33:33 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

int32_t neuromat_filter_lowpass_butterworth_compute_order(double fs, double fc, double gc);
  /* Computes an order {ord} such that the filter with knee frequency {fs}
    has value {gc} at frequency {fc}.  The frequencies {fs,fc} must be finite,
    positive, and distinct, and {gc} must be in {(0 _ 1)}.
    May fail if {fc} is too close to {fs} or {gc} is too small. */

double neuromat_filter_lowpass_butterworth_compute_fsup(double fs, int32_t ord);
  /* Returns a frequency {fsup} such that {|f|>fsup} implies 
    {neuromat_filter_lowpass_butterworth_eval(f,fs,ord)} is less than {10^{-15}}. */

double neuromat_filter_lowpass_butterworth_eval(double f, double fs, int32_t ord);
  /* A Butterworth lowpass filter with knee frequency {fs} and order
    {ord} (which must be positive).  If{fs} is zero, returns 0; else, if {f} is zero, returns 1.
    If {f} is negative, returns the same value as for {-f}. */

#endif
