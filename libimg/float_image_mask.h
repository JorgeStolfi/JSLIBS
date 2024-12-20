#ifndef float_image_mask_H
#define float_image_mask_H

/* Procedures to create mask images for convolution, filtering, etc.. */
/* Last edited on 2024-12-04 23:29:35 by stolfi */ 

#include <bool.h>
#include <float_image.h>

/* BASIC CLIPPING MASKS */

void float_image_mask_window(float_image_t *msk, int32_t ic, int32_t ord, bool_t round);
  /* Fills channel {ic} of image {msk} with a basic rectangular or
    circular window mask, with specified continuity order at
    the boundary.
    
    If {round} is FALSE, the mask spans the image's domain rectangle
    {R}. If {round} is TRUE, the mask spans the right ellipse {E}
    inscribed in {R}.
    
    The {order} parameter specifies the general profile of the sample
    values inside the window region {S} ({R} or {E}), as a function of
    the distance {r} from the center relative to the window radius:
    
      0  The constant value 1 throughout {S}.
      
      1  A down-turned quadratic {1-r^2}: value 1 and slope 0 at the center,
         value 0 and negative slope at the boundary of {S}.
         
      2 A biquadratic {(1-r^2)^2}: value 1 and slope 0 at the center,
        value 0 and slope 0 at the boundary of {S}. 
        
    Other values of {order} are illegal. In all cases, the mask falls
    smoothly to 0, with slope 0, at the boundary of {S}. To achieve
    this condition, when {order} is 0, the basic constant profile
    above is replaced by a cubic sigmoid starting 2 pixels inwards
    from the boundary of {S}. When {order} is 1, the basic downwards
    parabola above is replaced by and upwards parabola, starting 1
    pixel away from the boundary of {S}. */
  
/* MASK MODIFIERS 

  These procedures multiply a selected channel in a given mask {msk}
  by some 2D function {F(p) = F(x,y)}. In these procedures,
      
      {xc,yc} are the coordinates of the center of the image's domain;
      
      {sx,sy} are given width and height parameters, in pixels;
      
      {tx,ty} are relative coordinates of some point. Namely, if the
      coordinates of the point are {x,y}, then {tx = (x-xc)/sx} and
      {ty = (y-yc)/sy}.
      
  The function {F} is sampled at pixel centers, that is,
  {F(x+0.5,y+0.5)} where {x,y} are the pixel indices.  */
  
void float_image_mask_mul_gauss
  ( float_image_t *msk, 
    int32_t ic, 
    double sx, 
    double sy
  );
  /* Multiplies channel {ic} of image {msk} by a 2D Gaussian bell
    function {F(x,y) = G(tx)*G(ty)} where {G(z) = exp(-z^2/2)}. Note
    that {F} is 1 at the center of the image's domain. */
    
void float_image_mask_mul_power
  ( float_image_t *msk, 
    int32_t ic, 
    double sx, 
    double sy, 
    double pwr
  );
  /* Multiplies channel {ic} of image {msk} by the inverse power function
    {1/r^pwr}, modified to avoid infinities.
    
    More precisely, multiplies {msk} by {F(x,y) = 1/r^{pwr}}, where 
    {r = hypot(tx,ty,1.0)}. Note that {r} and {F} are 1 at the center of
    the image's domain.*/

/* MASK PROPERTIES */

typedef struct float_image_mask_stats_t
  { int32_t NX;            /* Number of pixel columns. */
    int32_t NY;            /* Number of pixel rows. */
    int32_t nINF;          /* Number of infinite samples. */
    int32_t nNAN;          /* Number of NAN samples. */
    double min;        /* Minimum sample value. */
    double max;        /* Max sample value. */               
    double avg;        /* Average of sample values. */       
    double dev;        /* Deviation of sample values. */     
    double ctr[2];     /* Barycenter of pixel centers. */ 
    double rad;        /* Radius enclosing all non-zero pixels. */
    double ext[2];     /* Max X,Y extent of non-zero pixels. */     
    double mmt[2];     /* X,Y Moments of nonzero pixels. */     
  } float_image_mask_stats_t;
  /* Statistical parameters of a mask image, namely:

    {.nINF}    Number of samples that are {±INF}.

    {.nNAN}    Number of {NAN} samples.

    {.min}     Minimum sample value.

    {.max}     Maximum sample value.

    {.avg}     Average of all sample values 

    {.dev}     Standard deviation of all sample values.

    {.ctr[ax]} Average position of nonzero samples, weighted by the
               absolute sample values.

    {.rad}     Max distance from the center {(xc,yc}} of the image's
               domain to the farthest point of any nonzero pixel.

    {.ext[ax]} Max distance along axis {ax} from the center {(xc,yc}}
               of the image's domain to the farthest point of 
               any nonzero pixel.

    {.mmt[ax]} Root-mean-square distance along axis {ax} from
               {(xc,yc}} and any nonzero sample, weighted by the
               absolute sample values.

  Infinite and {NAN} samples are counted in {.nINF} and {.nNAN},
  but ignored in  all the other statistics.

  Note that the {rad}, {ext}, and {mmt} are always computed relative
  to the image domain's ctr {(xc,yc)} and not to the barycenter
  {ctr[0..1]}. The sample with indices {ix,iy} is assumed to be a
  unit-side square centered at coordinates {(ix+0.5,iy+0.5)}. */

float_image_mask_stats_t float_image_mask_stats_get(float_image_t *msk, int32_t ic);
  /* Computes some general properties of channel {ic} of image {msk},
    viewed as a weight mask. */
    
void float_image_mask_stats_print(FILE *wr, float_image_mask_stats_t *S);
  /* Writes {S} to {wr}, in human-readable format. */

#endif
