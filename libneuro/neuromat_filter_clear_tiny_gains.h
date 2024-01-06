#ifndef neuromat_filter_clear_tiny_gains_H
#define neuromat_filter_clear_tiny_gains_H

/* NeuroMat tools for cleaning up Hartley filter coeff tables */
/* Last edited on 2024-01-05 15:51:25 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>
 
int32_t neuromat_filter_clear_tiny_gains(int32_t nf, double H[], double eps, double fsmp, bool_t verbose);
  /* Assumes that {H[0..nf-1]} are the Hartley filter gains for a signal
    with {nf} samples, sampled at {fsmp} samples per second. Resets to
    zero all conjugate pairs {H[kf],H[kf']} (where {kf'=(nf-kf)%nf})
    which are both smaller than {eps} in modulus (presumably due
    to roundoff errors in gain formula).
    
    Gain pairs that are not cleared but are not much greater than {eps}
    are scaled down so as to eliminate discontinutites in the frequency
    response between cleared and non-cleared sections.
    
    In particular, clears {H[0]}, and also {H[nf/2]} if {nf} is
    even, if they are smaller than {eps} in absolute value.
    
    The procedure returns the largest frequency index {kmax} in {0..nf/2}
    such that {H[kmax]} and/or {H[kmax']} are still is non-zero; or
    {-1} if all gains {H[0..nf-1]} became became all zeros as a result of
    this clearing. */

#endif
