/* See {neuromat_filter.h}. */
/* Last edited on 2023-12-16 17:56:53 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include <affirm.h>
#include <bool.h>

#include <neuromat_filter.h>

double neuromat_filter_hartley_basis_eval(int32_t n, int32_t f, int32_t t)
  { f = f % n; if (f < 0) { f = f + n; }
    t = t % n; if (t < 0) { t = t + n; }
    if (f == 0)
      { return 1.0/sqrt(n); }
    else
      { return cos(M_PI*((2.0*f*t)/n + 0.25))*2*M_SQRT1_2/sqrt(n); }
  }

complex neuromat_filter_fourier_basis_eval(int32_t n, int32_t f, int32_t t)
  { f = f % n; if (f < 0) { f = f + n; }
    t = t % n; if (t < 0) { t = t + n; }
    return cexp((2*f*t*M_PI*I)/n)/sqrt(n);
  }

void neuromat_filter_hartley_to_fourier(double Ha, double Hb, complex *Fa_P, complex *Fb_P)
  { double Wr = (Ha + Hb)/2;
    double Wi = (Hb - Ha)/2;
    (*Fa_P) = Wr - I*Wi;
    (*Fb_P) = Wr + I*Wi;
  }
  
void neuromat_filter_fourier_to_hartley(complex Fa, complex Fb, double *Ha_P, double *Hb_P)
  { demand(fabs(creal(Fa) - creal(Fb)) < 1.0e-12, "Fourier coeffs should have same real part");
    demand(fabs(cimag(Fa) + cimag(Fb)) < 1.0e-12, "Fourier coeffs should be conjugate");
    /* Force conjugation: */
    double Wr = 0.5*(creal(Fa) + creal(Fb));
    double Wi = 0.5*(cimag(Fa) - cimag(Fb));
    /* Convert Fourier coeff pair to Hartley coeff pair: */
    (*Ha_P)= Wr + Wi;
    (*Hb_P)= Wr - Wi;
  }
