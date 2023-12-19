#ifndef neuromat_eeg_spectrum_compute_H
#define neuromat_eeg_spectrum_compute_H

/* NeuroMat tools for computing power spectra of signals. */
/* Last edited on 2023-12-14 08:34:33 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

double **neuromat_eeg_spectrum_compute(int32_t n, int32_t ne, double **val, int32_t kfmax, bool_t verbose);
  /* Returns an array {pwr[0..ne-1][0..kfmax]} such that {pwr[i][k]} is the 
    power (rms value) of the Fourier components of electrode {i} with
    frequency index {k} (i.e. with {k} full cycles in the input sample
    sequence). The parameter {kfmax} must not exceed {n/2}.
    
    Unlike {neuromat_filter_apply}, this procedure uses no mirroring or
    padding before computing the discrete Fourier/Hartley transform of
    the signals; just windowing with a Hann (raised-cosine) window.
    Therefore, the power spectra will be blurred by the window's spectrum,
    which is 3 samples wide.
    
    If {verbose} is true prints some debugging info. */

#endif
