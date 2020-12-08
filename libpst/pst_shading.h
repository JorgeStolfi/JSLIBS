#ifndef pst_shading_H
#define pst_shading_H

/* pst_shading.h -- tools for shading computations. */
/* Last edited on 2007-01-06 22:15:37 by stolfi */

#include <bool.h>
#include <r3.h>
#include <float_image.h>
#include <argparser.h>

#include <pst_basic.h>
#include <pst_light.h>
#include <pst_lamp.h>
     
/* COMPUTING SYNTHETIC PHOTOS OF OBJECTS */

void pst_shading_add_diffuse_single
  ( float_image_t *NRM, 
    float_image_t *CLR,
    pst_lamp_t *src,
    float_image_t *IMG
  );
  /* Adds to a photo {IMG} of some scene the diffusive (Lambertian)
    term that results from illuminating it with the lamp {src}.
    
    Assumes that {NRM} is the normal map of the scene, and that {CLR}
    is its intrinsic color map. The number of channels in the lamp's
    power {pwr} must match the number of channels of {CLR} and {IMG}.
    If {CLR == NULL}, the intrinsic color of the surface is assumed to
    be white (1.0 in each channel).
    
    Adds to each pixel {IMG[p]} the amount {src.pwr[c]*coef*clr[c]} to
    channel {c} of every pixel, where {clr=CLR[p]} is the pixel's
    intrinsic color, and {coef} is the geometric factor,
    {coef=pst_lamp_geom_factor(NRM[p],src.dir,src.crad)}. Does not
    affect pixels where {NRM[p]} is the null vector, or where {CLR[p]}
    is zero.  */

void pst_shading_add_diffuse
  ( float_image_t *NRM, 
    float_image_t *CLR,
    pst_light_t *lht,
    float_image_t *IMG
  );
  /* Adds to {IMG} the diffusion shading term produced by the light field
    {lht}. Assumes a Lambertian surface with normal {NRM[p]} and
    intrinsic color (diffusion coefficient) {CLR[p]} at each pixel {p}. */

#define pst_shading_Lambertian_INFO \
  "  The visible surface of the scene is assumed to be Lambertian"\
  " (purely diffusive, with no specular or glossy reflection) and"\
  " of uniform lightness (diffusion coefficient)."

/* DIFFERENCE IMAGE */

float_image_t *pst_shading_difference_image
  ( float_image_t *NRM, 
    float_image_t *AIMG, 
    float_image_t *BIMG
  );
  /* Computes the sample-by-sample difference image {AIMG-BIMG},
    except that the difference is set to zero where the normal map
    {NRM} is zero. The three images must have the same size, and
    {AIMG} and {BIMG} must have the same number of channels. */


#endif
