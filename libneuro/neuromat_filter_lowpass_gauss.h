#ifndef neuromat_filter_lowpass_gauss_H
#define neuromat_filter_lowpass_gauss_H

/* NeuroMat basic gaussian profile filter response. */
/* Last edited on 2024-01-05 17:33:12 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

double neuromat_filter_lowpass_gauss_compute_sigma(double fc, double gc);
  /* Computes the deviation {sigma} of a Gaussian bell with zero mean
    and unit height that has value {gc} (which must be strictly between 0 and 1)
    at parameter {fc}. */

double neuromat_filter_lowpass_gauss_compute_fsup(double sigma);
  /* Returns a frequency {fsup} such that {|f|>fsup} implies 
    {neuromat_filter_lowpass_gauss_eval(f,sigma)} is less than {10^{-15}}. */

double neuromat_filter_lowpass_gauss_eval(double f, double sigma);
  /* A filter with unit-height, zero-mean, {sigma}-deviation
    Gaussian bell function {exp(-0.5*(f/sigma)^2)}. */

#endif
