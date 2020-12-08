#ifndef msm_image_H
#define msm_image_H

/* Image tools for {msm_match}. */
/* Last edited on 2013-10-22 01:32:50 by stolfilocal */

#define msm_image_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <stdint.h>
#include <float_image.h>

typedef struct msm_image_t
  { int NV[2];
    int scale;
    float_image_t *fim;
  } msm_image_t; 
   /* A virtual image with {NV[0]} by {NV[1]} pixels, implemented as a
     smaller float image {fim}, at a reduction factor {scale:1}. */

msm_image_t *msm_image_alloc(int NC, int NVX, int NVY, bool_t scale);
  /* Allocates a virtual image with {NC} channels, {NVX} columns and
    {NVY} rows. If {scale} is TRUE then the underlying {float_image}'s actual size is
    {ceil(NVX/ppp)} by {ceil(NVY/ppp)}, where {ppp} is chosen so
    that neither dimension exceeds some reasonable maximum. */

/* PAINTING DOTS 

  The procedures in this section use virtual coordinates; so the
  arguments {x,y} map to the pixel {((xr/scale,yr/scale)}, 
  where {xr,yr} are {x,y} reduced modulo {img->NV[0],img->NV[1]}. */

void msm_image_set(msm_image_t *img, int c, int x, int y, float v);
  /* Stores {v} into the sample in channel {c} of {img} corresponding to index pair
    {(x,y)} */

void msm_image_add(msm_image_t *img, int c, int x, int y, float v);
  /* Adds {v} to the sample in channel {c} of {img} corresponding to
    index pair {(x,y)}.  If that sample is {-oo}, sets it to {v}. */

void msm_image_max(msm_image_t *img, int c, int x, int y, float v);
  /* Stores {v} into the sample in channel {c} of {img} corresponding to
    index pair {(x,y)}, provided that it is greater than the current value
    of that sample. */

/* CHOOSING OUTPUT PIXEL RANGES

  The procedures in this section choose the low value {lo} and high value {hi} 
  that should be mapped to black and white (in either order) for better 
  visualization of a given image {img}.  In any case, the limits {lo} and {hi}
  are fudged as needed so as to ensure that {lo < hi}. */
  
void msm_image_compute_avg_dev_range(msm_image_t *img, int c, double *lo, double *hi, double ns);
  /* Computes the average {avg} and standard deviation {dev} of all
    finite sample values in channel {c} of {img}. Ignores samples that
    are {±INF} or {NAN}. The {lo} and {hi} values are set to
    {avg±ns*dev}, respectively.  
    
    If all valid samples have the same value, returns a tiny interval
    containing that value. If there are no valid sample values, returns
    a tiny interval centered at 0. */ 

void msm_image_compute_min_max_range(msm_image_t *img, int c, double *lo, double *hi);
  /* Computes the minimum {vmin} and maximum {vmax} of all finite
    sample values in channel {c} of {img}. Ignores samples that are
    {±INF} or {NAN}. The {lo} and {hi} values are set to {vmin} and
    {vmax}, respectively.  
    
    If all valid samples have the same value, returns a tiny interval
    containing that value. If there are no valid sample values, returns
    a tiny interval centered at 0. */ 

msm_image_t *msm_image_colorize(msm_image_t *gim, double maxf);
  /* Converts a single-channel image {gim} to a pseudocolor image {cim}.
    A value {v} of {gim} is mapped to blue if equal to {-maxf} or less,
    to red if  equal to {+maxf} or greater, and to white if close to zero.
    If {v} is {±INF}, it is treated as `undefined' and mapped to black. */

void msm_image_write_as_pnm
  ( msm_image_t *img, 
    double minv, 
    double maxv, 
    char *name,
    char *tag
  );
  /* Writes the image {img} to disk as a PGM or PPM image, with
    file name "{name}{tag}.pgm", maxval 255, converting affinely the
    samples from {[minv _ maxv]} to {[0..maxval]} with clipping.
    The output is PGM if {img} has 1 channel, and PPM if it has
    3 channels. */ 

void msm_image_normalize_and_write_as_pgm(msm_image_t *img, char *name, char *tag);
  /* Writes the image {img} to disk as a PGM file called
    "{name}{tag}.pgm". The entries of {img} are scaled so that the
    smallest finite pixel maps to gray, the largest pixel maps to
    black. In particular, pixels with value {-INF} map to white. */

#endif
