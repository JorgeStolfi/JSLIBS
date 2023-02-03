/* Test tools for {multifok_focus_op} and related funcs. */
/* Last edited on 2023-01-30 18:04:49 by stolfi */

#ifndef multifok_test_H
#define multifok_test_H

#define _GNU_SOURCE
#include <stdint.h>

#include <interval.h>
#include <r3.h>
#include <bool.h>
#include <frgb.h>
#include <float_image.h>

#include <multifok_scene.h>
#include <multifok_term.h>

void multifok_test_images_make
  ( int32_t NX, 
    int32_t NY, 
    multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    double zFoc, 
    double zDep, 
    int32_t NP,
    int32_t NR_min,
    float_image_t **cimg_P,
    float_image_t **simg_P,
    float_image_t **zimg_P,
    float_image_t **dimg_P
  );
  /* Creates images {csimg}, {shimg}, {azimg}, and {dzimg} of the given {scene},
    with simulated depth-of-focus blur. All four images will have {NX}
    columns of {NY} pixels.
    
    The {csimg} image, returned in {*cimg_P}, will have 3 channels, and
    shows the color of the scene as seen through a camera or muscope
    focused at {Z=zFoc} with focal depth {zDep}.
    
    The {shimg} image, returned in {*simg_P}, will have 1 channel, and
    shows the estimate sharpness of the image at each pixel, a number
    between 0 (totally fuzzy) and 1 (as sharp as possible given the
    pixel averaging).
    
    The {azimg} image, returned in {*zimg_P}, will have 1 channel, and
    shows the average {Z} coordinate of the scene at each pixel. Its
    values will range in {[zMin _ zMax]}, the range of {Z} coordinates
    in the scene, that is, {scene->box[2]}.

    The {dzimg} image, returned in {*zimg_P}, will have 1 channel, and
    shows the standard deviation of the {Z} coordinate of the scene
    inside each pixel. It should be in the range {[0 _ zdMax]} where
    {zdMax} is {max(|zMin-zFoc|, |zMax-zFoc|)}.

    The value of each pixel of the four images is computed by taking a
    regular grid of {NP} by {NP} sample points around the center of the
    pixel on the plane {Z = zFoc}. The procedure computes the apparent
    color {clr(p)}, sharpness {shr(p), average {Z} height {zht(p)}, and
    deviation of height {zdv(p)} at each sample point {p}, then averages
    those values with a 2D Hann window weight function.
    
    The values {clr(p)}, {shr(p)}, {zht(p)} at each sample point {p} are
    obtained by tracing at least {NR_min} rays through {p}, computing
    color {clr(R)}, the nominal sharpness {shr(R)}, and the {Z}
    coordinate {zht(R)} of each ray {R}, and averaging those values over
    all rays, with proper weights to simulate an apodizing mask. The
    value of {zdv(p)} is the deviation of {zht(R)} over those rays.
    
    The rays will have random directions that deviate from the vertical
    by maximum distance of about 2 pixels over a vertical travel of {zDep}. If
    {zDep=+INF}, then {NR_min} must be 1, and there will be only one
    strictly vertical ray per sampling piint {p}.
    
    The color {clr(R)} is determined by computing {r = pattern(x,y,z,kd,3,fs)}
    where {x,y,z} are the coordinates of the hit point, and {kd} is the 
    index of the object hit by the ray, or {-1} if the ray hit the background surface. */

/* BASIS AND TERM DATA I/O */

void multifok_test_read_term_indices_names_and_weights
  ( char *fname, 
    int32_t NB, 
    char *belName[],
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P, 
    int32_t *NT_P, 
    double **wt_P, 
    char ***termName_P,
    bool_t verbose
  );
  /* Reads from file "{fname}" a list of a certain number {NT} of
    quadratic local operators that are products of basis
    coefficients, and the weights of each term in a sharpness score.
    
    See {multifok_score_read_term_names_and_weights} for the file format.
    
    The information is returned in {*NP_P,*prix_P,*NT_P,*wt_P,*termName_P}. */
  
void multifok_test_write_term_names_and_weights(char *outPrefix, int32_t NT, char *termName[], double wt[]); 
  /* Writes the names {termName[0..NT-1]} and weights {wt[0..NT-1]} of the terms of a quadratic score
    formula to a file file "{outPrefix}-termNames.txt", one per line. */

void multifok_test_write_basis_elem_names(char *outPrefix, int32_t NB, char *belName[]);
  /* Writes the basis element names {belName[0..NB-1]} to a file file "{outPrefix}-bnames.txt"
    one per line. */

void multifok_test_write_term_names(char *outPrefix, int32_t NT, char *termName[])  ;
  /* Writes the basis element names {termName[0..NT-1]} to a file file "{outPrefix}-tnames.txt"
    one per line. */

/* TEST IMAGE I/O */

#define multifok_test_image_gamma 1.000
  /* Assumed encoding gamma of input and output images. */

#define multifok_test_image_bias 0.0327
  /* Assumed encoding bias of input and output images. */

float_image_t *multifok_test_read_scene_color_image(char *inPrefix, char *tag);
  /* Reads the simulated color image {csimg} from file "{inPrefix}{tag}-cs.ppm". */ 

float_image_t *multifok_test_read_sharpness_image(char *inPrefix, char *tag);
  /* Reads the actual sharpness image {shimg} from file "{inPrefix}{tag}-sh.pgm". */ 

float_image_t *multifok_test_read_zave_image(char *inPrefix, char *tag);
  /* Reads the actual relative scene {Z} coordinate image {azimg} from file "{inPrefix}{tag}-az.pgm"
    and scales the samples from {[0_1]} to {[0 _ ZMAX]} where {ZMAX} is 
    {multifok_scene_ZMAX}. */ 

float_image_t *multifok_test_read_zdev_image(char *inPrefix, char *tag);
  /* Reads the scene Z deviation image {dzimg} from file "{inPrefix}{tag}-dz.pgm"
    and scales the samples from {[0_1]} to {[0 _ ZMAX]} where {ZMAX} is 
    {multifok_scene_ZMAX}.. */ 

void multifok_test_write_scene_color_image(float_image_t *csimg, char *outPrefix, char *tag);
  /* Writes the color image {csimg} to file "{outPrefix}{tag}-cs.ppm". */

void multifok_test_write_sharpness_image(float_image_t *shimg, char *outPrefix, char *tag);
  /* Writes the sharpness image {shimg} to file "{outPrefix}{tag}-sh.pgm".
    Assumes that samples are in {[0 _ 1]}. */

void multifok_test_write_zavg_image(float_image_t *azimg, char *outPrefix, char *tag);
  /* Writes the scene {Z} average image {azimg} to file "{outPrefix}{tag}-az.pgm".
    Pixels are scaled from {[0 _ ZMAX]} to {[0 _ 1]} where {ZMAX}
    is {multifok_scene_ZMAX}. */

void multifok_test_write_zdev_image(float_image_t *dzimg, char *outPrefix, char *tag);
  /* Writes the scene {Z} deviation image {dzimg} to file "{outPrefix}{tag}-dz.pgm".
    The pixels are scaled from {[0 _ ZMAX]} to {[0_1]} where {ZMAX}
    is {multifok_scene_ZMAX}. */

void multifok_test_write_score_image(float_image_t *scimg, char *outPrefix, char *tag);
  /* Writes the estimated sharpness score image {scimg} to file "{outPrefix}{tag}-sc.pgm".
    Assumes that the sample values are in {[0 _ 1]}. */

void multifok_test_write_score_error_image(float_image_t *esimg, char *outPrefix, char *tag);
  /* Writes the sharpness error image {esimg} (computed {score} minus actual sharpness 
    {sharp}) to file "{outPrefix}{tag}-es.pgm".
    The sample values are rescaled from {[-1 _ +1]} to {[0_1]}. */
 
void multifok_test_write_estimated_Z_image(float_image_t *czimg, char *outPrefix, char *tag);
  /* Writes the image {czimg} with computed surface {Z} values 
    to file "{outPrefix}{tag}-cz.pgm". The sample values are rescaled from 
    {[0 _ ZMAX]} to {[0_1]}, where where {ZMAX}
    is {multifok_scene_ZMAX}. */
 
void multifok_test_write_Z_error_image(float_image_t *ezimg, char *outPrefix, char *tag);
  /* Writes the {Z} error image {ezimg} (computed {Z} minus actual {Z})
    to file "{outPrefix}{tag}-ez.pgm". The sample values are rescaled from 
    {[-vAbsMax_+vAbsMax]} to {[0_1]}, where {vAbsMax} is the max 
   absolute sample value. */
      
void multifok_test_write_reconstructed_color_image(float_image_t *crimg, char *outPrefix, char *tag);
  /* Writes the color image {crimg} as "{outPrefix}{tag}-cr.ppm" */
      
void multifok_test_write_color_error_image(float_image_t *erimg, char *outPrefix, char *tag);
  /* Writes the color image {erimg} as "{outPrefix}{tag}-er.ppm", mapping its samples
    from {[-1 _ +1]} to {[0_1]}. */

void multifok_test_write_window_average_image(float_image_t *avimg, char *outPrefix, char *tag);
  /* Writes the average sample image {avimg} as "{outPrefix}{tag}-av.pgm". */

void multifok_test_write_window_deviation_image(float_image_t *dvimg, char *outPrefix, char *tag);
  /* Writes the sample deviation image {dvimg} as "{outPrefix}{tag}-dv.pgm", mapping its samples
    from {[0 _ vMax]} to {[0_1]} where {vMax} is the max sample value. */

void multifok_test_write_normalized_image(float_image_t *nrimg, char *outPrefix, char *tag);
  /* Writes the locally normalized grayscale image {nrimg} as "{outPrefix}{tag}-nr.pgm", mapping its samples
    from {[-1 _ +1]} to {[0_1]}. */

void multifok_test_write_sample_weights_image(float_image_t *wsimg, char *outPrefix, char *tag);
  /* Writes the square {NW} by {NW} grayscale image {wsimg} as a file "{outPrefix}{tag}-ws.pgm" 
    The samples are rescaled from {0 _ vMax]} to {{0 _ 1]}, where {vMax}
    is the max sample value. */

void multifok_test_write_basis_elem_image(float_image_t *beimg, char *outPrefix, char *tag);
  /* Writes the square {NW} by {NW} grayscale image {beimg}, presumed to be a basis element,
    as a file "{outPrefix}{tag}-be.pgm".
    The samples if {beimg} are rescaled from {[-vAbsMax _ +vAbsMax]} to {[0 _ 1]}, where {vAbsMax}
    is the max absolute sample value. */

void multifok_test_write_basis_coeff_image(float_image_t *bcimg, char *outPrefix, char *tag);
  /* Writes the basis coeff image {bcimg} to file "{outPrefix}{tag}-bc.pgm".
    The sample values are rescaled from {[-vAbsMax _ +vAbsMax]} to {[0_1]}, where {vAbsMax}
    is the max absolute sample value. */

void multifok_test_write_basis_coeff_squared_image(float_image_t *bqimg, char *outPrefix, char *tag);
  /* Writes the squared basis coeff image {bqimg} to file "{outPrefix}{tag}-bq.pgm".
    The sample values are rescaled from {[0 _ vMax]} to {[0_1]}, where {vMax}
    is the average sample value plus a few standard deviations. */
    
void multifok_test_write_quadratic_term_image(float_image_t *tmimg, char *outPrefix, char *tag);
  /* Writes the quadrratic term values image {tmimg} to file "{outPrefix}{tag}-tm.pgm".
    The sample values are rescaled from {[-vR _ +vR]} to {[0_1]}, where {vR}
    is the max absolute sample value. */
    
void multifok_test_write_pixel_mask_image(float_image_t *mkimg, char *outPrefix, char *tag);
  /* Writes the pixel mask image {mkimg} to file "{outPrefix}{tag}-mk.pgm".
    The sample values are assumed to be in {[0_1]}. */

/* GENERIC IMAGE I/O */

float_image_t *multifok_test_read_color_image(char *inPrefix, char *tag, char *code, float vMin, float vMax);
  /* Reads a color image from file "{inPrefix}{tag}-{code}.ppm". It should have 3 channels. 
    After reading, all samples are remapped from {[0 _ 1]} to {[vMin _ vMax]}. */
  
float_image_t *multifok_test_read_grayscale_image(char *inPrefix, char *tag, char *code, float vMin, float vMax);
  /* Reads a grayscale image from file "{inPrefix}{tag}-{code}.pgm". It should have 1 channel. 
    After reading, all samples are remapped from {[0 _ 1]} to {[vMin _ vMax]}. */
    
void multifok_test_write_color_image(float_image_t *img, char *outPrefix, char *tag, char *code, float vMin, float vMax);
  /* Writes the color image {img} to file "{outPrefix}{tag}-{code}.ppm". It should have 3 channels. 
    Before writing, all samples are remapped from {[vMin _ vMax]} to {[0 _ 1]}. */
    
void multifok_test_write_grayscale_image(float_image_t *img, char *outPrefix, char *tag, char *code, float vMin, float vMax);
  /* Writes the grayscale image {img} to file "{outPrefix}{tag}-{code}.pgm". It should have 1 channel. 
    Before writing, all samples are remapped from {[vMin _ vMax]} to {[0 _ 1]}. */

#endif
