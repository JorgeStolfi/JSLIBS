#ifndef float_image_hartley_spectrum_H
#define float_image_hartley_spectrum_H

/* Tools for Hartley transform (real-valued Fourier-like transform). */
/* Last edited on 2024-12-05 10:29:42 by stolfi */ 

#include <stdint.h>

#include <bool.h>
#include <float_image.h>

void float_image_hartley_spectrum(float_image_t *H, float_image_t *P, bool_t center);
  /* Assumes that {H} is the Hartley transform of some image, and stores into {P}
    the corresponding power spectrum. See {float_image_hartley_spectrum_INFO} 
    below for details. */
    
#define float_image_hartley_spectrum_INFO \
  "  The power spectrum {P} of a single-channel image is derived from its Hartley" \
  " transform {H} by replacing each entry by the sum of squared" \
  " amplitudes of Hartley components with the same frequency" \
  " vector, ignoring the phase.\n" \
  "\n" \
  "  Specifically, let {NX} and {NY} be the column and row count" \
  " of {H}.  Element {P[fx,fy]} is generally {(H[fx,fy]^2 + H[gx,gy]^2)/2}, where" \
  " {gx,gy)} is the frequency pair congruent to {(-fx,-fy)} modulo {(NX,NY)};" \
  " unless {(gx,gy)} is the same as {(fx,fy)}, in which case" \
  " {P[fx,xy]} is just {H[fx,fy]^2}.  In any case, {P[gx,gy]} will be" \
  " equal to {P[fx,fy]}.\n" \
  "\n" \
  "  If {center} is true, the power spectrum image is shifted so that" \
  " the constant term, with frequency vector {(fx,fy)=(0,0)}, is placed" \
  " at the centermost pixel {(NX/2,NY/2)}.\n" \
  "\n" \
  "  The power spectrum of a multi-channel image is obtained by computing" \
  " the power spectrum of each channel separately."

#endif
