#ifndef float_image_hartley_H
#define float_image_hartley_H

/* Tools for Hartley transform (real-valued Fourier-like transform). */
/* Last edited on 2008-10-05 10:01:15 by stolfi */ 

#include <bool.h>
#include <float_image.h>

void float_image_hartley_transform(float_image_t *A, float_image_t *T);
  /* Stores into {T} the Hartley transform of the image {A}. The
    result has the same channel, column, and row counts as the input.
    Each channel is transformed independently.
    See {float_image_hartley_INFO} below for details. */
    
#define float_image_hartley_INFO \
  "  The /Hartley transform/ is based on the fact that any monochromatic" \
  " image {I} can be written as a linear combination of elements of the" \
  " corresponding /Hartley basis/.  The transform is another" \
  " monochromatic image, with the same column and row counts as {I}, where" \
  " each pixel value is the coefficient of {I} in that linear combination.\n" \
  "\n" \
  "  Each element of the Hartley basis is characterized by two integer" \
  " frequencies {FX,FY}, respectively /horizontal/ and /vertical/.  In" \
  " an image with {NX} columns and {NY} rows, the element is a sinusoidal wave" \
  " with amplitude {sqrt(2/(NX*NY))} that begins with phase {PI/4} at" \
  " pixel {(0,0)} (the bottom left pixel) and appears to propagate in the" \
  " direction {(NX/FX,NY/FY)}, with {FX} full cycles along each row, and" \
  " {FY} full cycles along each column.  In particular, the component with" \
  " frequencies {0,0} is a uniform image with value {sqrt(1/(NX*NY))};" \
  " and, when {NX} and {NY} are both even, the component with" \
  " frequencies {NX/2,NY/2} is a checkerboard with values {± sqrt(1/(NX*NY))}.  Thus," \
  " the frequencies {FX} and {FY} are measures in waves per image.\n" \
  "\n" \
  "  Since a Hartley basis element is defined only at the pixel locations, its" \
  " frequencies {FX,FY} are determined only modulo {NX,NY}; that is, the" \
  " frequency pair {FX+r*NX,FY+s*NY} denotes the same component" \
  " as {FX,FY}, for any integers {r,s}.  For indexing into the transform" \
  " image, one uses the /modulo-reduced/ frequencies {FX = FX mod NX} and" \
  " {FY = FY mod NY} (with mathematical \"mod\"), which range in {0..NX-1} and" \
  " {0..NY-1}, respectively.  Thus, pixel {[0,0]} (the bottom left corner) of" \
  " the transform is the coefficient of the constant term.\n" \
  "\n" \
  "  However, one must observe that a wave with horizontal frequency" \
  " {FX} in the range {floor(NX/2)+1..NX-1} is more naturally" \
  " interpreted as a wave with horizontal frequency {FX-NX}, in" \
  " the range {-ceil(NX/2)+1..-1}.  Therefore, for image analysis" \
  " purposes one should consider the /signed frequencies/" \
  " {FX,FY} which have been reduced modulo {NX,NY} to the range" \
  " {-floor(NX/2) .. +floor(NX/2)} and {-floor(NY/2) .. +floor(NY/2)}," \
  " respectively.  Note that, when {NX} is even, the" \
  " signed frequency of the wave with {NX/2} cycles per row is ambiguous" \
  " (either {-NX/2} or {+NX/2}); and ditto for {NY}.  In particular, if" \
  " both {NX} and {NY} are even, the component with frequencies {NX/2,NY/2} is" \
  " a checkerboard with samples {±1/sqrt(NX*NY)}\n" \
  "\n" \
  "  The /pixel wavelength/ of a Hartley component with" \
  " apparent frequencies {FX,FY} is {W = 1/sqrt((FX/NX)^2,(FY/NY)^2)}.  This" \
  " number, measured in pixels, ranges between {sqrt(2)} (for a checkerboard)" \
  " and {+oo} (for the constant component).  The /pixel frequency/" \
  " (in waves per pixel) is therefore {F = 1/W}, and varies" \
  " from 0 (for the constant component) to {sqrt(0.5)} for" \
  " the checkerboard.  Apart from the constant component," \
  " the longest-period waves along each axis have pixel wavelengths {NX} and" \
  " {NY}, and pixel frequencies {1/NX} and {1/NY}.\n" \
  "\n" \
  "  The Hartley basis is orthonormal with respect to the image dot product\n" \
  "    {<I|J> = SUM{I[x,y]*J[x,y] : x in {0..NX-1}, y in {0..NY-1}}}\n" \
  " and the corresponding image norm\n" \
  "    {|I| = sqrt(<I|I>) = sqrt(SUM{I[x,y]^2 : x in {0..NX-1}, y in {0..NY-1}})}\n" \
  " It follows that the sum of the squares of the pixels of {I} is equal to the" \
  " sum of the squares of its transform (Parseval's theorem).  Moreover, the Hartley" \
  " transform is its own inverse."

void float_image_hartley_wave(float_image_t *A, int fx, int fy, double amp);
  /* Fills {A} with an image that is a Hartley component of
    frequencies {fx} and {fy} --- namely, a wave with {fx} cycles
    along the X axis and {fy} waves along the Y axis, with amplitude
    {amp} and initial phase {PI/4}. The frequencies {fx,fy} are
    automatically reduced modulo the image dimensions to fit inside
    the image's domain.  The mean value is zero, and the total energy
    is {amp^2*nx*ny/2} */


#endif
