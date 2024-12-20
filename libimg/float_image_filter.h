#ifndef float_image_filter_H
#define float_image_filter_H

/* Tools for image filtering. */
/* Last edited on 2024-12-04 23:27:19 by stolfi */ 

#include <bool.h>
#include <r2.h>
#include <gauss_table.h>
#include <float_image.h>

void float_image_filter_gaussian_band
  ( float_image_t *F, 
    r2_t *wMin, 
    r2_t *wMax, 
    bool_t complement,
    bool_t verbose
  );
  /* Assumes that {F} is the Hartley transform of an image, and
    multiplies each component (F(fx,fy)} of {F} by a Gaussian band-pass
    filter mask {W(fx,fy)} defined by wavelength ranges {wMin,wMax};
    or, if {complement} is TRUE, by the complementary (band-kill)
    filter mask {1-W(fx,fy)}.
    
    The weight mask {W(fx,fy)} is designed to preserve all frequency
    components which, along each axis {ax}, have wavelengths between
    {wMin.c[ax]} and {wMax.c[ax]}. In particular, if {wMin.c[ax] ==
    0}, the mask {W} preserves the Fourier components with wavelengths
    substantially greater than {wMax.c[ax]}, including the image's
    mean level. Conversely, if {wMax.c[ax] == +INF}, the mask {W} will
    preserve only the Fourier components with wavelengths smaller than
    {wMin.c[ax]}. In both axes, and at both ends, the mask has a
    Gaussian-smooth (rather than sharp) cut-off.
    
    More precisely, the mask {W(fx,fy)} is {gHi(fx,fy)-gLo(fX,fY)}
    where {gHi,gLo} are two bidimensional Gaussian bells centered at
    frequency {0,0}, whose standard deviation, along each axis {ax},
    are {fHi[ax] = n[ax]/wMin.c[ax]} and {fLo[ax] = n[ax]/wMax.c[ax]},
    respectively; where {n[0]} and {n[1]} are the width and height of
    the image in pixels.
    
    If the wavelengths {wMin.c[ax]} and {wMax.c[ax]} are large
    compared to the image dimensions, the Gaussians {gHi,gLo} will
    fold over the edges of the transform image, as explained under
    {gauss_table_make}. */

double *float_image_filter_gaussian_freq_weights(int32_t n, double wRef);
  /* Creates a table of attenuation weights {w[0..n-1]} for the 
    Hartley transform terms with frequencies {0..n-1}
    in a signal with {n} samples.
    
    The weights correspond to a low-pass Gaussian filter with
    characteristic wavelength {wRef} and hence characteristic
    frequency {fRef = n/wRef}. The weight {w[i]} will be
    a Gaussian bell function of {i} with mean 0, max value 1, and 
    deviation {fRef}.  
    
    In particular, if {wRef} is 0, all weights will be 1. If {wRef} is
    finite but much larger than {n} (specifically, larger than
    {gauss_table_BIG_ARG*n}), then {w[0]} will be 1 and all other
    weights will be 0. As a special case, if {wRef} is {+INF}, (that
    is, {fRef = 0}), all weights will be zero. */

#endif
