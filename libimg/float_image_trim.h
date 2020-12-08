#ifndef float_image_trim_H
#define float_image_trim_H

/* Tools for trimming uniform background and margins from an image. */
/* Last edited on 2007-05-09 08:38:00 by stolfi */

#include <bool.h>
#include <irange.h>
#include <frgb.h>
#include <float_image.h>

void float_image_trim_background(float_image_t *I, double noise, double except, irange_t OBox[], bool_t paddable[], frgb_t padcolor[]);
  /* Identifies the margins of image {I} that look like background.

    Considers only the pixels of {I} that lie inside the rectangle
    {IBox}. Sets {padcolor[e]}, for {e} in {0..3}, to the dominant
    color along edge {e}.  Then repeatedly strips rows or columns of
    pixels along the edges of {IBox}, while they are mostly uniform
    and equal to {padcolor[e]}. The parameters {noise} and {except}
    are used to determine the dominant color and to check for
    uniformity. Returns in {OBox} the (possibly empty) sub-rectangle
    of {IBox} that remains at the end of this process.

    Also sets {paddable[e]} to TRUE iff at least one row or column of
    pixels was stripped along edge {e}.  Edges are numbered as 
    in {edge_t}. */ 

frgb_t float_image_dominant_color(float_image_t *I, irange_t *MBox, double noise);
  /* Obtains the dominant color {C} among the pixels of {I} contained
    in the rectangle {MBox}.  
    
    The result is such that {C} is the average of all pixels in {MBox}
    that lie within radius {noise} of {C}. Hopefully it also maximizes
    the number of such pixels. */

bool_t float_image_is_uniform(float_image_t *I, irange_t MBox[], frgb_t *clr, double noise, double except);
  /* Let {M} be the set of pixels in {MBox}, and {U} the subset
    whose colors lie more than {noise} away from {clr}.
    The procedure returns TRUE iff {|U|/|M| <= except}. */
    
#endif
