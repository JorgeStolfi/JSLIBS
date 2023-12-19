#ifndef neuromat_filter_lowpass_butterworth_H
#define neuromat_filter_lowpass_butterworth_H

/* NeuroMat basic Butterworth low-pass filter transfer function. */
/* Last edited on 2023-12-16 01:02:06 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

int32_t neuromat_filter_lowpass_butterworth_compute_order(double fs, double fc, double gc);
  /* Computes an order {ord} such that the filter with knee frequency {fs}
    has value {gc} at frequency {fc}. Requires {0 < fs < fc < +oo} and {0 < gc < 1}.
    May fail if {fc} is too close to {fs} or {gc} is too small. */

double neuromat_filter_lowpass_butterworth(double f, double fs, int32_t ord);
  /* A Butterworth (?) lowpass filter with knee frequency {fs} and order
    {ord}.  If {f} is zero, returns 1; else if {fs} is zero, returns 0. */

#endif
