#ifndef neuromat_filter_tabulate_hartley_gains_H
#define neuromat_filter_tabulate_hartley_gains_H

/* NeuroMat tools for building a Hartley filter coeff table. */
/* Last edited on 2024-01-03 04:19:35 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <complex.h>
#include <bool.h>

#include <neuromat_filter.h>

void neuromat_filter_tabulate_hartley_gains
  ( int32_t nf,
    double fsmp,
    neuromat_filter_gain_t *gain,
    double fsup,
    bool_t normalize,
    double H[],
    bool_t verbose
  );
  /* Fills {H[0..nf-1]} with the real-valued Hartley filter coefficients
    suitable for scaling a discrete Hartley transform with
    {nf} coefficients, for a discrete real-valued signal with period {nf}
    and sampling frequency {fsmp} (Hz). 
    
    The procedure assumes that, for all integer {kf} in {0..nf-1}, the complex gain
    {w(kf)} for the discrete Fourier component (complex sinusoidal) with
    frequency {f=kf*fsmp/nf} is {gain(f)}. It also assumes that {gain(f)} will
    be effectively zero for {|f|>=fsup}.  
    
    For each {kf} in {0..nf-1}, the real coefficients {H(kf)} and {H(kf')} together
    will encode the Fourier gains {w(kf)} and {w(kf')}, where {kf'= (n-kf)%nf}. 
    
    The procedure also accounts for the fact that frequencies greater than the Nyquist
    frequency {fsmp/2} will be aliased with lower frequencies. Namely,
    that {i*fsmp + f} will be effectively the same frequency as {f}.
    Therefore, the procedure will define {w(kf)} to be the sum of {gain(i*fsmp + f)} over all
    integer {i}, positive or negative, such that {|f + i*fsmp| < fsup}.
    That is, the function {gain}, whose support spans {fsmp} Hz or more,
    gets "folded over" onto the frequency range {[0_fsmp)}.
    
    In any case, the filter must produce a real-value filtered signal
    for any real-valued discrete signal input. This condition means that
    the folded gains {w(kf)} and {w(kf')} must be complex conjugates. In
    particular, {w(0)} must be pure real; and, if {nf} is even,
    {w(nf/2)} must be pure real too. This condition is automatically
    satisfied if {gain(f)} itself is pure real for any frequency {f}.

    This folding is suppressed if {fsup} is zero or negative.  In this case,
    the filter {w(kf)} will be simply {gain(f) + cconj(gain(fsmp-f))}.
    
    In any case, if {normalize} is true, will scale all resulting coeffs
    {H[0..nf-1]} by {1/wmax} where {wmax} is the maximum of {|w(kf)|}.
    Then the maximum absolute gain for any discrete sinusouldal
    component will be 1. However, this normalizarion is skipped if
    {wmax} is zero (that is, if {w(kf)=0} for all {kf}).
    
    If {verbose} is true, prints the tabulated filter gains. */
 

#endif
