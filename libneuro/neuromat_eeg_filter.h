#ifndef neuromat_eeg_filter_H
#define neuromat_eeg_filter_H

/* NeuroMat EEG Fourier/Hartley filtering. */
/* Last edited on 2023-12-16 04:01:41 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

void neuromat_eeg_filter_apply
  ( int32_t nt, 
    int32_t ne, 
    double **val, 
    int32_t tdeg, 
    bool_t tkeep,
    double *H,
    bool_t verbose
  );
  /* Applies a bandpass filter to the discrete signal set
    {val[0..nt-1][0..ne-1]} for {ne} electrodes sampled at {nt} equally
    spaced times. Assumes that the sampling rate is {fsmp} hertz
    (samples per second). Ignores any data channels ({val[0..nt-1][i]}
    with {i>=ne}), if present.
    
    The parameter {H} must be a table of {nt} coefficients that will be
    converted to Fourier gain coefficients and multiplied into the
    Fourier transform coefficients of each channel of {val}. Thus, in
    particular, {H[0]} is the gain for the zero-frequency (mean)
    component, and {H[nt/2]} is the gain for the highest frequency
    present.
    
    If the filter preserves phase of all sinusoids, the table {H} should
    be symmetric about the middle (that is, {H[nt-k] = H[k]} for {k} in
    {1..nt-1}).
    
    If {tdeg} is non-negative, the procedure fits a `trend' polynomial
    {P} of degree {tdeg} to each channel, and subtracts {P[0..n-1]}
    from the channel before applying the frequency filter. Then, if
    {tkeep} is true, the procedure adds the trend {P} back to the
    filtered signal. The trend polynomial is found with
    {neuromat_poly_fit_robust}. If {tdeg} is negative, the trend fitting
    is skipped, and {tkeep} is ignored.
    
    The trend fitting feature with {tdeg=0} can be used to preserve the
    mean level independently of the Fourier filter parameters. If the
    procedure is used with finite segments of non-periodic signals, one
    can specify {tdeg=1} or {tdeg=2} in order to avoid artifacts around
    the ends due to the implicit jump from {val[n-1][i]} to
    {val[0][i]}. Presumably there is no advantage to use {tdeg} larger
    than 2, since terms with degree 3 and higher should be well
    approximated by Fourier series.
    
    If {verbose} is true, prints various debugging info. */

#endif
