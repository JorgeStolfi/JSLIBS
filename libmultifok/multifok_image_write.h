/* Test tools for writing various maps related to multifocus stereo. */
/* Last edited on 2025-04-11 08:58:16 by stolfi */

#ifndef multifok_image_write_H
#define multifok_image_write_H

#include <stdint.h>

#include <float_image.h>

#include <multifok_image.h>

/* Every procedure that writes a PNG image file will first implicitly add
  of the input float image {img} a brightness scale, consisting of a row of
  rectangles where the samples are equally spaced values spanning the 
  assumed range of sample values of {img} (that depends on the procedure). 
  
  In the PNG file, pixel {[0,0]} of the input image {img} will be at the LOWER
  left corner, and the brightness scale will be added at the TOP. */

/* FRAME-INDEPENDENT PNG FILES */

void multifok_image_write_sample_weights(float_image_t *wsimg, char *outFolder);
  /* Writes the square {NW} by {NW} grayscale image {wsimg} as a file
    "{outFolder}/weights.png" The samples are rescaled from {0 _ vMax]} to
    {{0 _ 1]}, where {vMax} is the max sample value. */

void multifok_image_write_basis_kernels(uint32_t NB, float_image_t *bKer[], char *outFolder);
  /* For {kb} in {0..NB-1}, writes the square {NW} by {NW} grayscale image {bKer[kb]},
    presumed to be the kernel of a basis element, as a file "{outFolder}/basis-{BBB}.png"
    where {BBB} is {kb} formatted as "%03d". The samples
    of {bKer} are rescaled from {[-vAbsMax _ +vAbsMax]} to {[0_1]},
    where {vAbsMax} is the max absolute sample value in all elements. */

void multifok_image_write_pixel_mask(float_image_t *pSel, float_image_t *bgrd, char *outFolder);
  /* Writes the pixel mask image {pSel} to file "{outFolder}/pSel.png",
    The sample values are assumed to be in {[0_1]}.  
    
    If {bgrd} is {NULL}, the {pSel} image is written as is. 
    
    If {bgrd} is not {NULL}, it must be a color or grayscale image
    with the same size as {pSel}, and the samples of {pSel} must 
    be either 0 or 1.  The image written will be {0.75*pSel + 0.25*bgrd}. */

void multifok_image_write_selected_pixels(uint32_t NQ, i2_t pix[], float_image_t *bgrd, char *outFolder);
  /* Writes to file "{outFolder}/selected-pixels.png" the image {bgrd} with a cross
    centered at each of the pixels whose indices are listed in {pix[0..NQ-1]}.
    The {bgrd} image may be color or grayscale, and the sample values are assumed 
    to range in {{0_1]}. */
 
/* FRAME-SPECIFC PNG FILES */

void multifok_image_write_scene_view(float_image_t *sVal, char *frameFolder);
  /* Writes the simulated or real scene view image {sVal} to file "{frameFolder}/sVal.png".
    It may be color (3 channels) or grayscale (1 channel).
    The sample values are assumed to be in {[0_1]}. */

void multifok_image_write_height_average(float_image_t *hAvg, char *frameFolder, double hMin, double hMax);
  /* Writes the scene average {Z} image {hAvg} to file "{frameFolder}/hAvg.png".
    Pixels are scaled from {[hMin _ hMax]} to {[0_1]}. */

void multifok_image_write_height_deviation(float_image_t *hDev, char *frameFolder, double dMax);
  /* Writes the scene {Z} deviation image {hDev} to file "{frameFolder}/hDev.png".
    The pixels are scaled from {[0 _ dMax]} to {[0_1]}. */

void multifok_image_write_normal_average(float_image_t *sNrm, char *frameFolder);
  /* Writes the scene average normal vector image {sNrm} to file "{frameFolder}/sNrm.png".
    Normal coordinates {X}, {Y}, and {Z} are scaled from {[-1 _ +1]} to {[0_1]}
    and used as red, green, and blue, respectively. */

void multifok_image_write_sharpness(float_image_t *shrp, char *frameFolder);
  /* Writes the sharpness image {shrp} to file "{frameFolder}/shrp.png".
    Samples are assumed to range in {[0_1]} and thus are not 
    rescaled. */

void multifok_image_write_window_average(float_image_t *sAvg, char *frameFolder);
  /* Writes the image {sAvg} with window value averages as "{frameFolder}/sAvg.png".
    Assumes that the samples of {sAvg} are in {[0_1]}. */

void multifok_image_write_window_gradient(float_image_t *sGrd, char *frameFolder);
  /* Writes the image {sGrd} with window gradient moduli as
    "{frameFolder}/sGrd.png".  The number of channels {NC=sGrd.sz[0]} must be 3.
    The samples of {sGrd}, are implicitly mapped from {[-vMax _ vMax]} to {[0_1]}
    where {vMax} is the max absolute sample value. */

void multifok_image_write_window_deviation(float_image_t *sDev, char *frameFolder);
  /* Writes the image {sDev} with window value deviations as
    "{frameFolder}/sDev.png" The samples of {sDev}, which must be
    non-negative, are implicitly mapped from {[0 _ vMax]} to {[0_1]},
    where {vMax} is the max sample value. */

void multifok_image_write_focus_score(float_image_t *fVal, char *frameFolder);
  /* Writes the estimated focus sharpness score image {fVal} to file "{frameFolder}/fVal.png".
    Samples are implicitly mapped linearly from {[0 _ vMax]} to {[0_1]}]
    where {vMax} is the max sample value. */

void multifok_image_write_focus_score_error(float_image_t *fErr, char *frameFolder);
  /* Writes the focus sharpness error image {fErr} (estimated score
    {fVal} minus actual score {fInd}) to file "{frameFolder}/fErr.png".
    The sample values are rescaled from {[-1 _ +1]} to {[0_1]}. */
 
void multifok_image_write_estimated_height(float_image_t *zEst, char *frameFolder, double hMin, double hMax);
  /* Writes the image {zEst} with computed surface {Z} values 
    to file "{frameFolder}/zEst.png". The sample values are rescaled from 
    {[hMin _ hMax]} to {[0_1]}. */
 
void multifok_image_write_height_error(float_image_t *zErr, char *frameFolder, double hMin, double hMax);
  /* Writes the {Z} error image {zErr} (computed {Z} minus actual {Z})
    to file "{frameFolder}/zErr.png". The sample values are rescaled from 
    {[-hMag _ +hMag]} to {[0_1]}, where {hMag = hMax-hMin}. */
      
void multifok_image_write_merged_scene_view(float_image_t *sMrg, char *frameFolder);
  /* Writes the composite sharp scene view image {sMrg} as "{frameFolder}/sMrg.png".
    The image may be color (3 channels) or grayscale (1 channel).
    Assimes that samples are in {[0_1]}. */
      
void multifok_image_write_merged_scene_view_error(float_image_t *sErr, char *frameFolder);
  /* Writes the reconstruction sErrror image {sErr} (msErrged image {sMrg} minus 
    actual sharp view {sVal}) as "{frameFolder}/sErr.png", mapping its samples
    from {[-1 _ +1]} to {[0_1]}. */

/* FRAME-RELATED FNI FILES WITH WEIGHTS */

void multifok_image_write_fni_height_average(float_image_t *hAvgWht, char *frameFolder);
  /* Writes to file "{frameFolder}/hAvg.fni" the scene average {Z} image
    {hAvgWht}, with no conversion and no brightness scale.
    It must have two channels, channel 1 being the reliability weight. */

void multifok_image_write_fni_normal_average(float_image_t *sNrmWht, char *frameFolder);
  /* Writes to file "{frameFolder}/sNrm.fni" the scene average {Z} image
    {sNrmWht}, wwith no conversion and no brightness scale.
    It must have four channels, channel 3 being the reliability weight. */

/* OPERATOR-RELATED PNG FILES */
    
void multifok_image_write_normalized(float_image_t *sNrm, char *frameFolder);
  /* Writes the locally normalized grayscale image {sNrm} as
    "{frameFolder}/sNrm.png", mapping its samples from {[-1 _ +1]} to
    {[0_1]}. */

void multifok_image_write_basis_coeffs(uint32_t NB, float_image_t *bVal[], char *frameFolder, double bMax);
  /* For each {kb} in {0..NB-1}, writes the basis coeff image {bVal[kb]} 
    to file "{frameFolder}/bVal-{BBB}.png" where {BBB} is {kb} formatted as "%03d".
    The sample values are rescaled from {[-bMax _ +bMax]} to {[0_1]}. */

void multifok_image_write_basis_coeffs_squared(uint32_t NB, float_image_t *bSqr[], char *frameFolder, double bMax);
  /* For each {kb} in {0..NB-1}, writes the squared basis coeff image {bSqr[kb]}
    to file "{frameFolder}/bSqr-{BBB}.png" where {BBB} is {kb} formatted as "%03d".
    The sample values are rescaled from {[0 _ bMax^2]} to {[0_1]}. */

void multifok_image_write_quadratic_terms(uint32_t NT, float_image_t *tVal[], char *frameFolder, double tMax);
  /* For each {kt} in {0..NT-1}, wWrites the quadratic term values image {tVal[kt]} 
    to file "{frameFolder}/tVal-{TTT}.png" where {TTT} is {kt} formatted as "%03d".
    The sample values are rescaled from {[-tMax _ +tMax]} to {[0_1]}. */
   
/* GENERIC IMAGE FILE WRITING */
    
void multifok_image_write_png(float_image_t *img, char *dir, char *name, float vMin, float vMax);
  /* Writes the image {img} to file "{dir}/{name}.png". If {img} has 3
    or more channels, the resulting PNG will be a color image, where
    channels {0..2} will be interpreted as R, G, and B. Otherwise the
    resulting PNG will be a grayscale image, with the intensity taken
    from channel 0. In any case, all samples considered are remapped
    from {[vMin _ vMax]} to {[0_1]}. */
  
void multifok_image_write_fni(float_image_t *img, char *dir, char *name);
  /* Writes the image {img} to file  "{dir}/{name}.fni" in FNI format, without any conversion. */

#endif
