#ifndef neuromat_filter_lowpass_sgerf_H
#define neuromat_filter_lowpass_sgerf_H

/* NeuroMat basic sigmoid-shoulder filter response. */
/* Last edited on 2023-12-16 03:15:53 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

double neuromat_filter_lowpass_sgerf_compute_sigma(double fm, double fc, double gc);
  /* Computes the deviation {sigma} of an {erf} sigmoid that has value 0.5 at {fm} and 
    value {gc} at parameter {fc}.  The parameters {fm,fc} must be positive and distinct,
    and {gc} must be different from 0.5. */

double neuromat_filter_lowpass_sgerf(double f, double fm, double sigma);
  /* A real lowpass filter with a {erf}-like sigmoid response in log freq space, that
    is 0.5 at frequency {fm} and has deviation parameter {sigma}.
    
    The frequency {fm} must be finite and positive. */

#endif
