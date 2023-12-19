#ifndef neuromat_eeg_spectrum_write_H
#define neuromat_eeg_spectrum_write_H

/* NeuroMat tools for writing power spectra of signals. */
/* Last edited on 2023-12-14 08:55:49 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <bool.h>
    
void neuromat_eeg_spectrum_write
  ( FILE *wr, 
    int32_t n, 
    int32_t ne, 
    int32_t kfmax, 
    double fsmp, 
    double **pwr
  );
  /* Assumes that {pwr[0..ne-1][0..kfmax]} are the power spectra of {ne}
    electrodes computed from sampling at {n} moments at the sampling 
    frequency {fsmp} (in samples per second, i.e. hertz).
    
    Assumes that each frame {pwr[k]} is the power of the Fourier
    Hartley transform components with absolute discrete frequency {k},
    ranging from 0 to {kfmax} inclusive. If {k=0} or {2*k=ns} there is
    only one component, othwerwise the powers of the two components {k}
    and {n-k} are added together. Requires {kfmax} to be at most
    {n/2}.    
    
    Writes the spectra to disk as {kfmax+1} lines with {ne+4}
    columns
      
      "{k} {f} {flo[k]} {fhi[k]}  {pwr[0][k]} ... {pwr[ne-1][k]}"
    
    where {kflo..kfhi} is the discrete frequency (an integer), {f = fsmp*k/n} is
    the corresponding actual frequency (in hertz), and {[flo[k] _
    fhi[k]]} is an interval of frequencies suitable for histogram-style
    plots. The intervals are usually {[k-1/2 _ k+1/2]*fsmp/n}, except
    when {k=0}, in which case it is only {[0, 1/2]*fsmp/n}, and when
    {n} is even and {k = n/2}, in which case it is only {[k-1/2,
    k]*fsmp/n}. 
    
    Thus the width of the interval is proportional to the
    number of Fourier (or Hartley) components whose power was added to
    make {pwr[i][k]}. */

#endif
