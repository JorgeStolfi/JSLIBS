#ifndef neuromat_filter_H
#define neuromat_filter_H

/* NeuroMat basic defs for freq filter. */
/* Last edited on 2023-12-16 16:26:08 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <complex.h>
#include <bool.h>
    
/* HARTLEY AND FOURIER BASES

  The Hartley basis for vectors of {\RR^n} is {\eta[0..n-1]} where 
  {\eta[0][t] = 1} and {\eta[f][t] = sqrt(2/n)*cos(2*pi*f*t/n + pi/4)} 
  for all {f} in {1..n-1} and all {t} in {0..n-1}.  It is
  an orthonormal basis under the usual dot product 
  {<x|y> = SUM{ x[t]*y[t] : t \in {0..n-1}}}.
  
  The Fourier basis for vectors of {\RC^n} is {\phi[0..n-1]} where
  {\phi[f][t] = sqrt(1/n)*exp(i*2*pi*f*t/n)} for all {f}  and {t} in {0..n-1}.  
  It can also be written as {\phi[f][t] = sqrt(1/n)*\omg^{f*t}} where 
  {\omg = exp(i*2*pi/n)} is a primitive {n}th root of unity.  It is
  an orthonormal basis under the dot product 
  {<x|y> = SUM{ x[t]*conj(y[t]) : t \in {0..n-1}}}.
  
  By convention, the frequency index {[f]} and a vector element
  (time, sample) index {[t]} are always taken modulo {n}.
  
  The real frequency of a basis element {\eta[f]} or {\ph[f]} is 
  {f} if {f} is in {0..n/2}, or {f-n} if {f} is in {n/2..n-1}.
  Thus the real frequency ranges in {-n/2..n/2}. If {n} is even,
  the components {\eta[n/2]} and {\phi[n/2]} have ambiguour real
  frequency, either {-n/2} or {+n/2}. */
  
double neuromat_filter_hartley_basis_eval(int32_t n, int32_t f, int32_t t);
  /* The value {\eta[f][t]} at sample index {t} of the element {\eta[f]} of the Hartley basis 
    for {\RR^n}.  It is periodic on both {f} and {t} with period {n}. */
   
complex neuromat_filter_fourier_basis_eval(int32_t n, int32_t f, int32_t t);
  /* The value {\phi[f][t]} at sample index {t} of the element {\phi[f]} of the Fourier basis 
    for {\CR^n}.  It is periodic on both {f} and {t} with period {n}. */
   
/* HARTLEY TO FOURIER TRANSFORM CONVERSION

  The procedures in this section refer to the Hartley and Fourier transforms.
  
  The Hartley transform of a vector {x[0..n-1]} of {\RR^n}
  is the coefficients of {x} on the Hartley basis {\eta},
  namely the real vector {H[0..n-1] = (\CH x)} such that
  {x[t] = SUM{ H[f]*\eta[f][t] : f \in {0..n-1}}} for all {t}.
  
  The Hartley transform {H = \CH x} can be computed by brute 
  force as {H[f] = <x|\eta[f]> = <\eta[f]|x>.

  The Fourier transform of a vector {x[0..n-1]} of {\CR^n} is the
  coefficients of {x} on the Fourier basis {\eta}, namely the complex
  vector {F[0..n-1] = (\CF x)} such that {x[t] = SUM{ F[f]*\phi[f][t] :
  f \in {0..n-1}}} for all {t}.
  
  The Fourier transform {F = \CF x} can be computed by brute force as
  {F[f] = <x|\phi[f]>}. Note that, for generic vectors {x,y\in\RC^n},
  {<y|x>} is the complex conjugate of {<x|y>}; so {\phi[f]|x>} is not
  {F[f]} but its conjugate.
  
  If {x} is a vector of {\RR^n}, its Fourier transform {F = (\CF x)} has
  the property that {F[f]} and {F[f']} are complex conjugates, for all
  {f} in {0..n-1}; where {f' = (n-f) % n}. With the convention that
  indices are taken modulo {n}, this is simply saying that {F[f] =
  conj(F[-f])}. for all integer {f}. */

void neuromat_filter_hartley_to_fourier(double Ha, double Hb, complex *Fa_P, complex *Fb_P);
  /* Given elements elements {Ha=(\CH x)[fa]} and {Hb=(\CH
    x)[fb]} of a Hartley transform {\CH x} with complementary indices {fa} and {fb = (nf-fa) %
    nf}, computes the equivalent coefficients {Fa=(\CF x)[fa]} and
    {Fb=(\CF x)[fb]} in the Fourier transform, which are returned in
    {*Fa_P} and {*Fb_P}. These coefficients will be complex
    conjugates. */
    
void neuromat_filter_fourier_to_hartley(complex Fa, complex Fb, double *Ha_P, double *Hb_P);
  /* Given elements of a Fourier transform {Fa=(\CF x)[fa]} and {Fb=(\CF
    x)[fb]} {\CH x} with complementary indices {fa} and {fb = (nf-fa) %
    nf}, computes the equivalent coefficients {Fa=(\CF x)[fa]} and
    {Fb=(\CF x)[fb]} in the Fourier transform, which are returned in
    {*Fa_P} and {*Fb_P}. These coefficients will be complex
    conjugates. */


typedef complex neuromat_filter_gain_t(double f);
  /* Type of a procedure that returns the (complex) gain of a complex sinusoidal
    with frequency {f}; that is, the gain for the complex signal component
    {x(t)=exp(2*pi*I*f*t)} for all real {t}. */

#endif
