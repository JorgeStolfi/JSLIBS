#ifndef float_image_hdyn_H
#define float_image_hdyn_H

/* Tools for high-dynamic-range mixing of PGM/PBM images. */
/* Last edited on 2017-06-20 20:47:08 by stolfilocal */ 

#include <r2.h>
#include <r2x2.h>
#include <r3x3.h>
#include <bool.h>
#include <uint16_image.h>
#include <float_image.h>

/* PROCEDURES FOR ZERO OFFSET AND UNIT SATURATION

  The following procedures take as input a list {img[0..NI-1]} of
  images. They look only at channel {c} of each image, which is
  supposed to be a photo of the same scene with the same lighting,
  taken with different camera gains (aperture/exposure/speed
  settings).   Each photo is assumed to be contaminated by independent
  additive noise and may contain overexposed areas.

  These procedures assume that the input images have been normalized
  so that samples with value {vmin[i]} correspond to zero luminosity, and
  samples with value {vmax[i]} are overexposed. In both cases,
  allowance is made for noise.  Between these two limits, the 
  sample value {V} of each image is assumed to be a linear function
  of the true luminosity {Y}, namely {V = gain*Y + offset}, where
  {gain} and {offset} are different for each image.

  The procedures estimate the unknown parameters of this model:
  namely, the gain {gain[i]} of each image, the value offset {offset[i]}, 
  the deviation {sigma[i]} of the noise component, and the true luminosity
  of each pixel.  The latter is stored in channel {c} of the 
  corresponding pixel of image {omg}.
  
  If the {verbose} parameter is TRUE, the procedures print out
  various diagnostic messages. 
*/

  
void float_image_hdyn_estimate_gains_offsets_sigmas_values
  ( int c,                /* Channel. */
    int NI,               /* Number of input images. */
    float_image_t *img[], /* Input images. */
    double vmin[],        /* Value of underexposed pixels. */
    double vmax[],        /* Value of overexposed pixels. */
    bool_t verbose,       /* TRUE for diagnostics. */
    double gain[],        /* (OUT) Estimated gain of each image. */
    double offset[],      /* (OUT) Estimated value offset of each image. */
    double sigma[],       /* (OUT) Estimated noise deviation. */
    float_image_t *omg    /* (OUT) Estimated true values. */
  );
  /* Iteratively estimates {gain[i],offset[i],sigma[i],omg[c,x,y]} for each image {i}
     and each pixel {x,y}. 
     
     The luminosities should approximatelt span the interval {[0_1]}.
     May store into {omg} either {-INF} (underexposed), {+INF}
     (overexposed) or {NAN} (indeterminate) if there is no valid
     information about the pixel. */
  
void float_image_hdyn_estimate_gains_offsets
  ( int c,                /* Channel. */
    int NI,               /* Number of input images. */
    float_image_t *img[], /* Input images. */
    double vmin[],        /* Value of underexposed pixels. */
    double vmax[],        /* Value of overexposed pixels. */
    double sigma[],       /* Assumed noise deviation of each image. */
    bool_t verbose,       /* TRUE for diagnostics. */
    double gain[],        /* (OUT) Estimated gain of each image. */
    double offset[]       /* (OUT) Estimated value offset of each image. */
  );
  /* Assumes that the deviations {sigma[0..NI-1]} are known, and estimates
    the gains {gain[0..NI-1]} and the offsets {offset[0..NI-1]}.
    Adjusts them so that the implied luminosities for pixels that are
    exposed in some image span the interval {[0_1]}. */

void float_image_hdyn_estimate_values
  ( int c,                /* Channel. */
    int NI,               /* Number of input images. */
    float_image_t *img[], /* Input images. */
    double vmin[],        /* Value of underexposed pixels. */
    double vmax[],        /* Value of overexposed pixels. */
    double gain[],        /* Assumed gain of each image. */
    double offset[],      /* Assumed value offset of each image. */
    double sigma[],       /* Assumed noise deviation of each image. */
    bool_t verbose,       /* TRUE for diagnostics. */
    float_image_t *omg    /* (OUT) Estimated luminances. */
  );
  /* Assumes that the gains {gain[0..NI-1]}, offsets {offset[0..NI-1]} and 
    deviations {sigma[0..NI-1]} are known, estimates the
    true luminosities {omg[c,x,y]} for all pixels.
    
    Pixels that are overexposed in all images are set to {+INF}.
    Pixels that are underexposed in all images are set to {-INF}.
    Indeterminate pixels are set to {NAN}. Other pixels should be
    strictly between 0 and 1. */

void float_image_hdyn_estimate_sigmas
  ( int c,                /* Channel. */
    int NI,               /* Number of input images. */
    float_image_t *img[], /* Input images. */
    double vmin[],        /* Value of underexposed pixels. */
    double vmax[],        /* Value of overexposed pixels. */
    double gain[],        /* Assumed gain of each image. */
    double offset[],      /* Assumed value offset of each image. */
    float_image_t *omg,   /* Assumed luminances. */
    bool_t verbose,       /* TRUE for diagnostics. */
    double sigma[]        /* (OUT) Estimated noise deviation. */
  );
  /* Assumes that the gains {gain[0..NI-1]}, offsets{offset[0..NI-1]} and the
    true luminosities {omg[c,x,y]} are known; estimates the
    deviations {sigma[0..NI-1]} for all images. */

/* FOR LATER */

typedef struct float_image_hdyn_params_t
  { /* Adjusted parameters: */
    double Ymin;      /* Underexposure threshold. */
    double Ymax;      /* Overexposure threshold. */
    double Vmin;      /* Min reliable image value. */
    double Vmax;      /* Max reliable image value. */
    double sigma;     /* Noise level. */
  } float_image_hdyn_params_t;
  /*  Transfer function parameters for one channel of one image.  The transfer function gives 
    the sample value {V(Y)} corresponding to a `true' scene radiance {Y}.  The value
    is assumed to be {Vmin} if {Y <= Ymin}, {Vmax} if {Y > Ymax}, and a linear ramp
    between those two points.  In addition, {V(Y)} is assumed to be contaminated with
    noise of standard deviation {sigma}. */

#endif
