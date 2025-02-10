/* Test tools for {multifok_focus_op} and related funcs. */
/* Last edited on 2025-02-01 16:45:38 by stolfi */

#ifndef multifok_image_H
#define multifok_image_H

#include <stdint.h>

#include <interval.h>
#include <i2.h>
#include <r2.h>
#include <r3.h>
#include <bool.h>
#include <frgb.h>
#include <float_image.h>

#include <multifok_scene.h>
#include <multifok_term.h>

/* TEST IMAGE I/O 
  
  All procedures in this section writte their image files so that pixel {0,0}
  is at the LOWER left corner. */

#define multifok_image_gamma 1.000
  /* Assumed encoding gamma of input and output images. */

#define multifok_image_bias 0.0327
  /* Assumed encoding bias of input and output images. */

/* FRAME-INDEPENDENT IMAGES */

void multifok_image_sample_weights_write(float_image_t *wsimg, char *outFolder);
  /* Writes the square {NW} by {NW} grayscale image {wsimg} as a file
    "{outFolder}/weights.png" The samples are rescaled from {0 _ vMax]} to
    {{0 _ 1]}, where {vMax} is the max sample value. */

void multifok_image_basis_kernels_write(uint32_t NB, float_image_t *bKer[], char *outFolder);
  /* For {kb} in {0..NB-1}, writes the square {NW} by {NW} grayscale image {bKer[kb]},
    presumed to be the kernel of a basis element, as a file "{outFolder}/basis-{BBB}.png"
    where {BBB} is {kb} formatted as "%03d". The samples
    of {bKer} are rescaled from {[-vAbsMax _ +vAbsMax]} to {[0 _ 1]},
    where {vAbsMax} is the max absolute sample value in all elements. */

void multifok_image_pixel_mask_write(float_image_t *pSel, float_image_t *bgrd, char *outFolder);
  /* Writes the pixel mask image {pSel} to file "{outFolder}/pSel.png",
    The sample values are assumed to be in {[0_1]}.  
    
    If {bgrd} is {NULL}, the {pSel} image is written as is. 
    
    If {bgrd} is not {NULL}, it must be a color or grayscale image
    with the same size as {pSel}, and the samples of {pSel} must 
    be either 0 or 1.  The image written will be {0.75*pSel + 0.25*bgrd}. */

void multifok_image_selected_pixels_write(uint32_t NQ, i2_t pix[], float_image_t *bgrd, char *outFolder);
  /* Writes to file "{outFolder}/selected-pixels.png" the image {bgrd} with a cross
    centered at each of the pixels whose indices are listed in {pix[0..NQ-1]}.
    The {bgrd} image may be color or grayscale. */
 
/* FRAME-SPECIFC IMAGES */

void multifok_image_scene_view_write(float_image_t *sVal, char *frameFolder);
  /* Writes the simulated or real scene view image {sVal} to file "{frameFolder}/sVal.png".
    It may be color (3 channels) or grayscale (1 channel).
    The sample values are assumed to be in {[0 _ 1]}. */

float_image_t *multifok_image_scene_view_read(char *frameFolder);
  /* Reads the simulated or real scene view image {sVal} from file "{frameFolder}/sVal.png".
    It may be color (3 channels) or grayscale (1 channel).
    The sample values will be to be in {[0 _ 1]}.   */ 


void multifok_image_height_average_write(float_image_t *hAvg, char *frameFolder, double hMin, double hMax);
  /* Writes the scene average {Z} image {hAvg} to file "{frameFolder}/hAvg.png".
    Pixels are scaled from {[hMin _ hMax]} to {[0 _ 1]}. */

float_image_t *multifok_image_height_average_read(char *frameFolder, double hMin, double hMax);
  /* Reads the scene average {Z} image {hAvg} from file "{frameFolder}/hAvg.png"
    and scales the samples from {[0_1]} to {[hMin _ hMax]} */ 


void multifok_image_height_deviation_write(float_image_t *hDev, char *frameFolder, double hMin, double hMax);
  /* Writes the scene {Z} deviation image {hDev} to file "{frameFolder}/hDev.png".
    The pixels are scaled from {[0 _ hMag]} to {[0_1]} where {hMag = (hMax-hMin)/2}. */

float_image_t *multifok_image_height_deviation_read(char *frameFolder, double hMin, double hMax);
  /* Reads the scene Z deviation image {hDev} from file "{frameFolder}/hDev.png"
    and scales the samples from {[0_1]} to {[0 _ hMag]} {hMag = (hMax-hMin)/2}. */ 


void multifok_image_normal_average_write(float_image_t *sNrm, char *frameFolder);
  /* Writes the scene average {Z} image {sNrm} to file "{frameFolder}/sNrm.png".
    Pixels are scaled from {[-1 _ +1]} to {[0 _ 1]}. */

float_image_t *multifok_image_normal_average_read(char *frameFolder);
  /* Reads the scene average {Z} image {sNrm} from file "{frameFolder}/sNrm.png"
    and scales the samples from {[0_1]} to {[-1 _ +1]} */ 


void multifok_image_sharpness_write(float_image_t *shrp, char *frameFolder);
  /* Writes the sharpness image {shrp} to file "{frameFolder}/shrp.png".
    Samples are assumed to range in {[0 _ 1]} and thus are not 
    rescaled. */

float_image_t *multifok_image_sharpness_read(char *frameFolder);
  /* Reads the actual sharpness image {shrp} from file "{frameFolder}/shrp.png".
    Sample values will be in {[0 _ 1]}. */ 


void multifok_image_window_average_write(float_image_t *sAvg, char *frameFolder);
  /* Writes the image {sAvg} with window value averages as "{frameFolder}/sAvg.png".
    Assumes that the samples of {sAvg} are in {[0 _ 1]}. */

void multifok_image_window_gradient_write(float_image_t *sGrd, char *frameFolder);
  /* Writes the image {sGrd} with window gradient moduli as
    "{frameFolder}/sGrd.png".  The number of channels {NC=sGrd.sz[0]} must be 3.
    The samples of {sGrd}, are implicitly mapped from {[-vMax _ vMax]} to {[0 _ 1]}
    where {vMax} is the max absolute sample value. */

void multifok_image_window_deviation_write(float_image_t *sDev, char *frameFolder);
  /* Writes the image {sDev} with window value deviations as
    "{frameFolder}/sDev.png" The samples of {sDev}, which must be
    non-negative, are implicitly mapped from {[0 _ vMax]} to {[0_1]},
    where {vMax} is the max sample value. */

void multifok_image_normalized_write(float_image_t *sNrm, char *frameFolder);
  /* Writes the locally normalized grayscale image {sNrm} as
    "{frameFolder}/sNrm.png", mapping its samples from {[-1 _ +1]} to
    {[0_1]}. */


void multifok_image_focus_score_write(float_image_t *fVal, char *frameFolder);
  /* Writes the estimated focus sharpness score image {fVal} to file "{frameFolder}/fVal.png".
    Samples are implicitly mapped linearly from {[0 _ vMax]} to {[0 _ 1]}]
    where {vMax} is the max sample value. */

void multifok_image_focus_score_error_write(float_image_t *fErr, char *frameFolder);
  /* Writes the focus sharpness error image {fErr} (estimated score
    {fVal} minus actual score {fInd}) to file "{frameFolder}/fErr.png".
    The sample values are rescaled from {[-1 _ +1]} to {[0_1]}. */

 
void multifok_image_estimated_height_write(float_image_t *zEst, char *frameFolder, double hMin, double hMax);
  /* Writes the image {zEst} with computed surface {Z} values 
    to file "{frameFolder}/zEst.png". The sample values are rescaled from 
    {[hMin _ hMax]} to {[0_1]}. */
 
void multifok_image_height_error_write(float_image_t *zErr, char *frameFolder, double hMin, double hMax);
  /* Writes the {Z} error image {zErr} (computed {Z} minus actual {Z})
    to file "{frameFolder}/zErr.png". The sample values are rescaled from 
    {[-hMag _ +hMag]} to {[0_1]}, where {hMag = hMax-hMin}. */
    
      
void multifok_image_merged_scene_view_write(float_image_t *sMrg, char *frameFolder);
  /* Writes the composite sharp scene view image {sMrg} as "{frameFolder}/sMrg.png".
    The image may be color (3 channels) or grayscale (1 channel).
    Assimes that samples are in {[0 _ 1]}. */
      
void multifok_image_merged_scene_view_error_write(float_image_t *sErr, char *frameFolder);
  /* Writes the reconstruction sErrror image {sErr} (msErrged image {sMrg} minus 
    actual sharp view {sVal}) as "{frameFolder}/sErr.png", mapping its samples
    from {[-1 _ +1]} to {[0_1]}. */

    
void multifok_image_basis_coeffs_write(uint32_t NB, float_image_t *bVal[], char *frameFolder, double bMax);
  /* For each {kb} in {0..NB-1}, writes the basis coeff image {bVal[kb]} 
    to file "{frameFolder}/bVal-{BBB}.png" where {BBB} is {kb} formatted as "%03d".
    The sample values are rescaled from {[-bMax _ +bMax]} to {[0_1]}. */

void multifok_image_basis_coeffs_squared_write(uint32_t NB, float_image_t *bSqr[], char *frameFolder, double bMax);
  /* For each {kb} in {0..NB-1}, writes the squared basis coeff image {bSqr[kb]}
    to file "{frameFolder}/bSqr-{BBB}.png" where {BBB} is {kb} formatted as "%03d".
    The sample values are rescaled from {[0 _ bMax^2]} to {[0_1]}. */

void multifok_image_quadratic_terms_write(uint32_t NT, float_image_t *tVal[], char *frameFolder, double tMax);
  /* For each {kt} in {0..NT-1}, wWrites the quadratic term values image {tVal[kt]} 
    to file "{frameFolder}/tVal-{TTT}.png" where {TTT} is {kt} formatted as "%03d".
    The sample values are rescaled from {[-tMax _ +tMax]} to {[0_1]}. */
   
/* GENERIC IMAGE I/O */
    
void multifok_image_write(float_image_t *img, char *dir, char *name, float vMin, float vMax);
  /* Writes the color image {img} to file "{dir}/{name}.png".
    It may have any number of channels. Before writing, all samples are remapped
    from {[vMin _ vMax]} to {[0 _ 1]}. */

float_image_t *multifok_color_read(char *dir, char *name, float vMin, float vMax);
  /* Reads a color image from file "{dir}/{name}.png". It may have any number of  channels. 
    After reading, all samples are remapped from {[0 _ 1]} to {[vMin _ vMax]}. */
     
void multifok_image_write(float_image_t *img, char *dir, char *name, float vMin, float vMax);
  /* Writes the grayscale image {img} to file
    "{dir}-{name}.png". It should have 1 channel. BValfore
    writing, all samples are remapped from {[vMin _ vMax]} to {[0 _ 1]}. */
 
float_image_t *multifok_image_read(char *dir, char *name, float vMin, float vMax);
  /* Reads a grayscale image from file "{dir}-{name}.png". It
    should have 1 channel. AftsErr reading, all samples are remapped from
    {[0 _ 1]} to {[vMin _ vMax]}. */

/* OTHER TOOLS */

void multifok_image_draw_crosses(float_image_t *img, int32_t ch, uint32_t NQ, i2_t pix[], float val);
  /* Draws open crosses into image {img} at the positions listed in {pix[0..NQ-1]}
    The crosses are drawn with value {val} into channel {ch}, and with value 0 in every other
    channels. */

#endif
