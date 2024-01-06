#ifndef neuromat_filter_bandpass_butterworth_H
#define neuromat_filter_bandpass_butterworth_H

/* NeuroMat  bandpass filter gain function with Butterworth shoulders in log freq. */
/* Last edited on 2024-01-05 17:35:15 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

void neuromat_filter_bandpass_butterworth_compute_parms
  ( double flo0,  double flo1,
    double fhi1,  double fhi0,
    double *fs_lo_P, double *fs_hi_P, int32_t *ord_P,
    double *fsup_P,
    bool_t verbose
  );
  /* Computes the elbow fequencies parameters {fs_lo} and {fs_hi} and
    the order {ord} so that the bandpass filter function {W} defined by
    {neuromat_filter_bandpass_butterworth_eval} is about zero at {flo0} and
    {fhi0} and about 1 between {flo1} and {fhi1}. Also computes the
    upper complete cutoff frequency {fsup}. Returns the parameters in
    {*fs_lo_P,*fs_hi_P,*ord_P,*fsup_P}.
  
    Both {flo0} and {fhi0} may be zero, in which case {fs_lo} and
    {ord_lo} will be zero. Also {fhi1} and {fhi0} may be both {+INF}, in
    which case {fs_hi} and {ord_hi} will be both {+INF}. Otherwise the
    key frequencies must satisfy {0<flo0<flo1<fhi1<fhi0<+INF}. */

double neuromat_filter_bandpass_butterworth_eval(double f, double fs_lo, double fs_hi, int32_t ord);
  /* Evaluates a bandpass filter function {W} at frequency {f}, whose
    shoulders are Butterworth quasi-linear ramps in log freq x log gain space, 
    with elbows at {fs_lo} and {fs_hi}, respectively, and order {ord}.  
    
    Namely, the filter is {W(f) = (1-W_lo(f))*W_hi(f)} where
    {W_lo(f)=neuromat_filter_lowpass_butterworth_eval(f,fs_lo,ord)} and
    {W_hi(f)=neuromat_filter_lowpass_butterworth_eval(f,fs_hi,ord)}. The sign
    of {f} is ignored, so that {W} is even (symmetric about zero).
    
    If {fs_lo} is zero, {W_lo(f)} is always zero so {W} is the same as
    the low-pass filter {W_hi}. Symmetrically, if {fs_hi} is {+INF},
    {W_hi(f)} is always 1 so that {W} is the same as the high-pass filter
    {1-W_lo}. */

#endif
