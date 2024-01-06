#ifndef neuromat_filter_bandpass_sgerf_H
#define neuromat_filter_bandpass_sgerf_H

/* NeuroMat bandpass filter gain function with erf-sigmoid shoulders in log freq. */
/* Last edited on 2024-01-05 17:34:20 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

void neuromat_filter_bandpass_sgerf_compute_parms
  ( double flo0,  double flo1,
    double fhi1,  double fhi0,
    double *fm_lo_P, double *sigma_lo_P,
    double *fm_hi_P, double *sigma_hi_P,
    double *fsup_P,
    bool_t verbose
  );
  /* Computes the parameters {fm_lo,sigma_lo} and {fm_hi,sigma_hi} so
    that the bandpass filter function {W} defined by
    {neuromat_filter_bandpass_sgerf_eval} is about zero at {flo0} and {fhi0}
    and about 1 between {flo1} and {fhi1}. Also computes the upper
    complete cutoff frequency {fsup}. Returns the parameters in
    {*fm_lo_P,*sigma_lo_P,*fm_hi_P,*sigma_hi_P,*fsup_P}.
  
    Both {flo0} and {fhi0} may be zero, in which case {fm_lo} and
    {sigma_lo} will be zero. Also {fhi1} and {fhi0} may be both {+INF},
    in which case {fm_hi} and {sigma_hi} will be both {+INF}. Otherwise
    the key frequencies must satisfy {0<flo0<flo1<fhi1<fhi0<+INF}. */

double neuromat_filter_bandpass_sgerf_eval(double f, double fm_lo, double sigma_lo, double fm_hi, double sigma_hi);
  /* Evaluates a bandpass filter function {W} at frequency {f}, whose
    shoulders are {erf}-like sigmoids in log freq space, that are about
    0.5 when {f=fm_lo} or {f=fm_hi}, and drop from 1 to 0 with deviation
    parameters {sigma_lo} and {sigma_hi}, respectively.  
    
    Namely, the filter is {W(f) = (1-W_lo(f))*W_hi(f)} where
    {W_lo(f)=neuromat_filter_lowpass_sgerf_eval(f,fm_lo,sigma_lo)} and
    {W_hi(f)=neuromat_filter_lowpass_sgerf_eval(f,fm_hi,sigma_hi)}. The sign
    of {f} is ignored, so that {W} is even (symmetric about zero).
    
    If {fm_lo} is zero, {W_lo(f)} is always zero so {W} is the same as
    the low-pass filter {W_hi}. Symmetrically, if {fm_hi} is {+INF},
    {W_hi(f)} is always 1 so that {W} is the same as the high-pass filter
    {1-W_lo}. */

#endif
