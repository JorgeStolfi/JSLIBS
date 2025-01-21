#ifndef pst_shading_H
#define pst_shading_H

/* pst_shading.h -- tools for shading computations. */
/* Last edited on 2025-01-17 04:07:39 by stolfi */

#include <bool.h>
#include <r3.h>
#include <frgb.h>
#include <float_image.h>
#include <argparser.h>

#include <pst_basic.h>
#include <pst_light.h>
#include <pst_lamp.h>
       
typedef frgb_t pst_shading_func_t(r3_t *nrm);
  /* Type of a function that, given the normal {nrm} of the scene's
    surface at some point, returns the shading of the surface that
    point; namely, what would be the apparent color of the surface at
    that point if the surface intrinsic color was white (albedo = 1). It
    should take into account the normal, the local illumination
    conditions, and the surface finish (Lambertian, glossy, etc) but not
    the the intrinsic color (albedo). */

/* COMPUTING SYNTHETIC PHOTOS OF OBJECTS */

float_image_t* pst_shading_make_image
  ( int32_t NC, int32_t NX, int32_t NY,
    pst_normal_func_t *normal,
    pst_albedo_func_t *albedo,
    pst_shading_func_t *shading
  );
  /* Creates an image of a scene, given
    procedures {normal} and {albedo} that compute, respectively, the
    normal direction and intrinsic color of the point of the surface
    visible at a given image point; and a procedure {shading} that
    computes the shading (color of a white surface) of a surface point
    given its normal.  
    
    The number of channels {NC} must be 3 or 4. The procedure creates an
    image of the specified dimensions and then calls {pst_shading_paint}
    (q.v.) to paint the scene into it. */

double pst_shading_paint
  ( float_image_t *img,
    int32_t xlo, int32_t xhi,
    int32_t ylo, int32_t yhi,
    pst_normal_func_t *normal,
    pst_albedo_func_t *albedo,
    pst_shading_func_t* shading
  );
  /* Paints into {img} an image of a scene given
    procedures {normal} and {albedo} that compute, respectively, the
    normal direction and intrinsic color of the point of the surface
    visible at a given image point; and a procedure {shading} that
    computes the shading (color of a white surface) of a surface point
    given its normal.
    
    The painting is limited to the sub-image consisting of columns
    {xlo..xhi} and rows {ylo..yhi} of the image, clipped to the actual
    image bounds.
    
    The number {NC} of channels in the image must be 3 or 4.
    The three components of the computed color are 
    pained into channels {0..2}.  If {NC} is 4, the procecdure stores
    into channel 3 the opacity mask {alpha} (see below).
    
    The color painted into each pixel will be a weighted average
    of the colors of several points inside and around the pixel, excluding
    those where the functions {normal}, {albedo}, or {shading} return
    invalid values.  The opacity {alpha} of each pixel is the fraction
    of sample points which were included in the average.  This opacity is
    used to mix the computed color with the original color of that
    pixel in {img}, and is also stored in channel 3 if {NC} is 4.
    
    In any case, the returned result is the sum of all the opacities
    {alpha} of the painted pixels. */

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
  ( float_image_t *AIMG, 
    float_image_t *BIMG, 
    float_image_t *MSK,
    float_image_t *NRM
  );
  /* Computes the sample-by-sample difference image {AIMG-BIMG}.
    The two images must have the same column and row counts {NX} and {NY}, which will be those of the 
    resulting image.  They may have different depths (channel counts) {NCA,NCB};
    the depth of the result will be the {NCC=min(NCA,NCB)}.
    
    Optionally multiplies the samples in the resulting diference image by
    pixel weights obtained in any of the following ways:
     
       * If the depth {NCA} and {NCB} are different, then samples of channel {NCC+1} of the deepest image 
         will be used as pixel weights.
         
       * If {MSK} is not {NULL}, it must be a single-channed image image with {NX} columns
         and {NY} rows.  The value of each pixel is the 
         
       * If {NRM] is not {NULL}, it must be a three-channel image with {NX} columns and {NY} rows
         containing the assumed surface normal in each pixel.  The pixel weight implied by this 
         argument is zero if the normal is not finite or is {(0,0,0)}, 1 otherwise.
         
    Weights that are negative or not finite are treated as zero, and
    weights greater than 1 are treated as 1. If two or more weights are
    available for each pixel, they are multiplied together. */

#endif
