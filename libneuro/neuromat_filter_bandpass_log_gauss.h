#ifndef neuromat_filter_bandpass_log_gauss_H
#define neuromat_filter_bandpass_log_gauss_H

/* NeuroMat basic bandpass filter transfer function with log-Gaussian profile. */
/* Last edited on 2024-01-06 08:04:01 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

void neuromat_filter_bandpass_log_gauss_compute_parms
  ( double flo0,  double flo1,
    double fhi1,  double fhi0,
    double *fm_P, int32_t *np_P, double *sigma_P, double *mag_P,
    double *fsup_P,
    bool_t verbose
  );
  /* Computes the parameters {fm,np,sigma,mag} so that the bandpass
    function {W} defined by {neuromat_filter_bandpass_log_gauss} is
    about zero at {flo0} and {fhi0} and about constant between {flo1}
    and {fhi1}, with max value 1. Also computes the upper complete
    cutoff frequency {fsup}. Returns the parameters in
    {*fm_P,*np_P,*sigma_P,*fsup_P,*mag_P}.
  
    The parameters must satisfy {0<flo0<flo1<fhi0<fhi1<+INF}. */

double neuromat_filter_bandpass_log_gauss_eval(double f, double fm, int32_t np, double sigma, double mag);
  /* A filter whose profile is {mag} times the sum of {np} Gaussians in
    log freq space Each Gaussian has unit height, deviation {sigma}, and
    mean {log(fm) + 2*k*sigma} for {k} in {0..np-1}. Note that the sign
    of {sigma} affects the position of the successive pulses relative to
    the first one.
    
    The parameter {fm} must be positive. If {f} is zero, or {np} is
    zero, returns zero. If {f} is negative, returns the same result as
    for {-f}. */

#endif
