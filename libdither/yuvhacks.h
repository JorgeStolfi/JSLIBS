/* yuvhacks.h - tools for perceptual color distance  */
/* Last edited on 2024-12-25 12:44:11 by stolfi */

#include <stdint.h>

void set_yuv_matrix(void);
  /* Computes coefficients for luminance/chrominance transformation,
    for RGB coordinates in the range {[0..maxval]}.
    See rgb_to_yuv below. */

void rgb_to_yuv(int32_t R, int32_t G, int32_t B, uint32_t maxval, float *yp, float *up, float *vp);
  /*  Computes luminance {*yp} and chrominance {*up}, {*vp} for the
    given RGB color {(R,G,B)}. 
    
    The results are a projective function of {(R,G,B)}, chosen so that
    Euclidean distance in {(y,u,v)} space corresponds more closely to
    perceptual color difference.

    For {(R,G,B)} in the {[0..maxval]^3} cube, the luminance {*yp} 
    lies in {[0..1]}. */

float rgb_to_y(int32_t R, int32_t G, int32_t B, uint32_t maxval);
  /*  Computes luminance only of pixel (R,G,B), as explined in
    rgb_to_yuv. */
