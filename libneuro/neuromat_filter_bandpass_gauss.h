#ifndef neuromat_filter_bandpass_gauss_H
#define neuromat_filter_bandpass_gauss_H

/* NeuroMat bandpass filter gain function defined as difference of gaussians in linear freq space. */
/* Last edited on 2024-01-05 17:34:42 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

void neuromat_filter_bandpass_gauss_compute_parms
  ( double flo0,  double flo1,
    double fhi1,  double fhi0,
    double *sigma_lo_P, double *sigma_hi_P, double *mag_P,
    double *fsup_P,
    bool_t verbose
  );
  /* Computes the parameters {sigma_lo}, {sigma_hi},and {mag} so that
    the bandpass filter function {W} defined by
    {neuromat_filter_bandpass_gauss_eval} is about zero at {flo0} and {fhi0}
    and about 1 between {flo1} and {fhi1}. Also computes the upper
    complete cutoff frequency {fsup}. Returns the parameters in
    {*sigma_lo_P,*sigma_hi_P,*fsup_P}.
  
    Both {flo0} and {fhi0} may be zero, in which case {sigma_lo} will be
    zero. Also {fhi1} and {fhi0} may be both {+INF}, in which case
    {sigma_hi} will be {+INF}. Otherwise the key frequencies must
    satisfy {0<flo0<flo1<fhi1<fhi0<+INF}. */

double neuromat_filter_bandpass_gauss_eval(double f, double sigma_lo, double sigma_hi, double mag);
  /* Evaluates a bandpass filter function {W} at frequency {f}, defined
    by the difference of two Gaussian bell functions with standard deviations
    {sigma_0} and {sigma_1}, respectively.  
    
    Namely, the filter is {W(f) = mag*(1-W_lo(f))*W_hi(f)} where
    {W_lo(f)=neuromat_filter_lowpass_gauss_eval(f,sigma_lo)} and
    {W_hi(f)=neuromat_filter_lowpass_gauss_eval(f,sigma_hi)}. The sign of {f}
    is ignored, so that {W} is even (symmetric about zero).
    
    If {sigma_lo} is zero, {W_lo(f)} is always zero, so {W} is the same
    as the low-pass filter {W_hi}. Symmetrically, if {sigma_hi} is
    {+INF}, {W_hi(f)} is always 1, so that {W} is the same as the
    high-pass filter {1-W_lo}. */

#endif
