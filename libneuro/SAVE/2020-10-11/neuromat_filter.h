#ifndef neuromat_filter_H
#define neuromat_filter_H

/* NeuroMat filtering and spectral analysis tools. */
/* Last edited on 2013-11-21 02:53:06 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <complex.h>
#include <bool.h>

typedef complex neuromat_filter_t(int kf, int nf, double fsmp);
  /* Type of a procedure that returns the (complex) gain of a component
    with frequency index {kf} in a discrete Fourier transform with {nf}
    coefficients. Requires {kf} in {0..nf-1}. Assumes the original data
    was sampled with frequency {fsmp} (hertz). Clients will assume that
    {gain((nf-kf)%nf,nf,fsmp)} is the complex conjugate of
    {gain(kf,nf,fsmp)}. */

void neuromat_filter_apply
  ( int nt, 
    int ne, 
    double **val, 
    double fsmp, 
    int tdeg, 
    bool_t tkeep, 
    neuromat_filter_t *gain,
    bool_t verbose
  );
  /* Applies a bandpass filter to the discrete signal set
    {val[0..nt-1][0..ne-1]} for {ne} electrodes sampled at {nt} equally
    spaced times. Assumes that the sampling rate is {fsmp} hertz
    (samples per second). Ignores any data channels ({val[0..nt-1][i]}
    with {i>=ne}), if present. Each Fourier coefficient with frequency
    index {kf} is multiplied by the complex number {gain(kf,nf,fsmp)}.
    
    If {tdeg} is non-negative, the procedure fits a `trend' polynomial {P} of
    degree {tdeg} to the data, and subtracts {P} from the data before applying
    the frequency filter.  Then, if {tkeep} is true, the procedure adds
    the trend {P} back to the filtered signal.  The trend polynomial is
    found with {neuromat_poly_fit_robust}.  If {tdeg} is negative, 
    the trend fitting is skipped, and {tkeep} is ignored.
    
    The trend fitting feature with {tdeg=0} can be used to preserve the
    mean level independently of the Fourier filter parameters. If the
    procedure is used with finite segments of non-periodic signals, one
    can specify {tdeg=1} or {tdeg=2} in order to avoid artifacts around
    the ends due to the implicit jump from {val[nt-1][i]} to
    {val[0][i]}. Presumably there is no advantage to use {tdeg} larger
    than 2, since terms with degree 3 and higher should be well
    approximated by Fourier series.
    
    If {verbose} is true prints the gain table and other debugging info. */
    
/* LOWPASS FILTERS */

/* 
  Each of these procedures computes a lowpass filter transfer function
  that has gain approximately 1 for frequencies below {fa} and approximately 
  zero for frequencies above {fb}, with smooth transitions in between.

  For the real-valued functions, one should assume that the phase is zero,
  so that the filter is time-symmetric, zero-delay.
  
  The frequencies {fa,fb} must be non-negative with {fa<=fb}. If
  {fa=fb=+INF} the filter is a no-op (returns gain 1 for all
  frequencies). If {fa==fb} the filter returns gain 0 when {f==fa}. */

double neuromat_filter_lowpass_gauss(double f, double fc, double gc, double fsmp);
  /* A Gaussian-profile filter with unit gain at frequency {f==0} and
    gain {gc} at frequency {f==fc}. The fiter is actually folded over
    itself by assuming a sampling frequency {fsmp}. */
    
double neuromat_filter_lowpass_biquadratic(double f, double fmax);
  /* A biquadratic filter that rolls off from gain 1 at {f==0} to gain zero at {f==fmax}. */ 

double neuromat_filter_lowpass_sigmoid(double f, double fa, double fb);
  /* A real lowpass filter with a log-sigmoid response, that has gain
     {0.999} at {f=fa} and {0.001} at {f=fb}. */

double neuromat_filter_lowpass_butterworth(double f, double fc, int n);
  /* A Butterworth (?) lowpass filter with cutoff frequency {fc} and order {n}.
    If {fc} is zero, returns 1 if {f} is zero, 0 otherwise. */

complex neuromat_filter_lowpass_cbutterworth(double f, double fc, int n);
  /* The complex Butterworth filter with cutoff frequency {fc} and order {n}.
    If {fc} is zero, returns 1 if {f} is zero, 0 otherwise. */

/* SPECTRAL ANALYSIS */

double **neuromat_filter_compute_spectra(int nt, int ne, double **val, int kfmax, bool_t verbose);
  /* Returns an array {pwr[0..ne-1][0..kfmax]} such that {pwr[i][kf]} is the 
    power (rms value) of the Fourier components of electrode {i} with
    frequency index {kf} (i.e. with {kf} full cycles in the input sample
    sequence). The parameter {kfmax} must not exceed {nt/2}.
    
    Unlike {neuromat_filter_apply}, this procedure uses no mirroring or
    padding before computing the discrete Fourier/Hartley transform of
    the signals; just windowing with a Hann (raised-cosine) window.
    Therefore, the power spectra will be blurred by the window's spectrum,
    which is 3 samples wide.
    
    If {verbose} is true prints some debugging info. */
    
void neuromat_filter_write_spectra
  ( FILE *wr, 
    int nt, 
    int ne, 
    int kfmax, 
    double fsmp, 
    double **pwr
  );
  /* Assumes that {pwr[0..ne-1][0..kfmax]} are the power spectra of {ne}
    electrodes computed from sampling at {nt} moments at the sampling 
    frequency {fsmp} (in samples per second, i.e. hertz).
    
    Assumes that each frame {pwr[kf]} is the power of the Fourier
    Hartley transform components with absolute discrete frequency {kf},
    ranging from 0 to {kfmax} inclusive. If {kf=0} or {2*kf=ns} there is
    only one component, othwerwise the powers of the two components {kf}
    and {nt-kf} are added together. Requires {kfmax} to be at most
    {nt/2}.    
    
    Writes the spectra to disk as {kfmax+1} lines with {ne+4}
    columns
      
      "{kf} {f} {flo[kf]} {fhi[kf]}  {pwr[0][kf]} ... {pwr[ne-1][kf]}"
    
    where {kflo..kfhi} is the discrete frequency (an integer), {f = fsmp*kf/nt} is
    the corresponding actual frequency (in hertz), and {[flo[kf] _
    fhi[kf]]} is an interval of frequencies suitable for histogram-style
    plots. The intervals are usually {[kf-1/2 _ kf+1/2]*fsmp/nt}, except
    when {kf=0}, in which case it is only {[0, 1/2]*fsmp/nt}, and when
    {nt} is even and {kf = nt/2}, in which case it is only {[kf-1/2,
    kf]*fsmp/nt}. 
    
    Thus the width of the interval is proportional to the
    number of Fourier (or Hartley) components whose power was added to
    make {pwr[i][kf]}. */

#endif
