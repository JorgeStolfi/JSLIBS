#ifndef neuromat_filter_lowpass_sgerf_H
#define neuromat_filter_lowpass_sgerf_H

/* NeuroMat low-pass filter gain function with erf-sigmoid shoulder in log freq. */
/* Last edited on 2024-01-05 17:32:50 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

double neuromat_filter_lowpass_sgerf_compute_sigma(double fm, double fc, double gc);
  /* Computes the deviation {sigma} of an {erf} sigmoid that has value 0.5 at {fm} and 
    value {gc} at parameter {fc}.  The parameters {fm,fc} must be positive and distinct,
    and {gc} must be different from 0.5. */

double neuromat_filter_lowpass_sgerf_compute_fsup(double fm, double sigma);
  /* Returns a frequency {fsup} such that {|f|>fsup} implies 
    {neuromat_filter_lowpass_sgerf_eval(f,fm,sigma)} is zero. */

double neuromat_filter_lowpass_sgerf_eval(double f, double fm, double sigma);
  /* A real lowpass filter {W} with a {erf}-like sigmoid response in log freq space, that
    is 0.5 when {f=fm} and drops from 1 to 0 with deviation parameter {sigma}. The frequency {fm} 
    must be finite and positive. The function is 1 for {f = 0},
    and ignores the sign of {f}. */

#endif
