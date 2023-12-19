#ifndef neuromat_filter_tabulate_hartley_gains_H
#define neuromat_filter_tabulate_hartley_gains_H

/* NeuroMat tools for building a Hartley filter coeff table. */
/* Last edited on 2023-12-16 04:10:35 by stolfi */

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
    {nf} coefficients, for a discrete real-valued signal with period {nf} and 
    sampling frequency {fsmp} (Hz).
    
    If {fsup} is negative or zero, the procedure assumes that, for all
    integer {kf} in {0..nf-1}, the complex gain {w(kf)} for the discrete
    Fourier component (complex sinusoidal) with frequency {f=kf*fsmp/nf}
    is {gain(f)}. It implicitly assumes that {gain(f+fsmp) = gain(f)}
    for all real {f}.
    
    if {fsup} is positive, the procedure assumes that {gain(f)} is zero
    for {|f| > fsup}. In this case, if {fsup >= fsmp/2}, the procedure
    assumes that the gain {w(kf)} of the discrete Fourier component with
    frequency {f=kf*fsmp/nf} is the sum of {gain(f + i*fsmp)} over all
    integer {i}, positive or negative, such that {|f + i*fsmp| < fsup}.
    That is, the function {gain}, whose support spans {fsmp} Hz or more,
    gets "folded over" into the frequency range {[0_fsmp)}.
    
    In any case, the resuting complex gains {w(kf)}, for {kf} in
    {0..nf-1} must imply a real-value filtered signal for any
    real-valued discrete signal input. This condition means that the
    folded gains {w(kf)} and {w(kf')} must be complex conjugates, for
    all {kf} in {0..nf-1}, where {kf'} is {(nf-kf) % nf}. This condition
    is automatically satisfied if {gain(f)} itself is pure real for any
    frequency {f}. In any case, particular, {w(0)} must be pure real;
    and, if {nf} is even, {w(nf/2)} must be pure real too. In any case,
    for other {kf} in {1..nf/2}, the two conjugate gains {w(kf),w(kf')}
    are converted to the corresponding pair of real-valued Hartley gains
    {H[kf],H[kf']}.
    
    In any case, if {normalize} is true, will scale all resulting 
    coeffs {H[0..nf-1]} by {1/wmax} where {wmax} is the maximum of {|w(kf)|}.
    Then the maximum absolute gain for any discrete 
    sinusouldal component will be 1.  However, this normalizarion is
    skipped if {wmax} is zero (meaning {w(kf)=0} and {H[kf]=0} for all {kf}).
    
    If {verbose} is true, prints the tabulated filter gains. */
 

#endif
