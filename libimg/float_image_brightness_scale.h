#ifndef float_image_brightness_scale_H
#define float_image_brightness_scale_H

/* Shrink an image by one row and one column. */
/* Last edited on 2025-02-02 02:14:10 by stolfi */ 

#include <bool.h>
#include <float_image.h>

void float_image_brightness_scale_paint
  ( float_image_t *A,
    int32_t xMin, int32_t xMax,  
    int32_t yMin, int32_t yMax,  
    float vMin, float vMax
  );
  /* Paints a brightness scale onto image {A}, consisting of a
    horizontal row of rectangular chips with values increasing linearly
    from {vMin} to {vMax}.
    
    The scale will span columns {xMin..xMax} and rows {yMin..yMax} of
    {A}. If the range {xMin..xMax} is not too small, there will be a
    one-pixel margin with value {vMid=(vMin+vMax)/2} between the chips
    and around the scale.
    
    The scale will normally have 11 chips, but the number of chips may
    be reduced if {A} has too few columns. At worst there may be just
    three chips with values {vMin}, {vMid}, {vMax}; or just one chip
    with value {vMid}. */

float_image_t *float_image_brightness_scale_add(float_image_t *A, bool_t loY, float vMin, float vMax);
  /* Makes a copy {R} of image {A} appending to it a brightness scale.
    
    The column and channel count of {A} will be preserved. The scale
    will span the whole with of {R}. The row count will be increased by
    the height of the brightness scale, that will be approximately
    proportional to the row count of {A}. The scale will be inserted
    before the first row of {A} if {loY} is true, or after the last row
    if {loY} is false. The scale is painted with
    {float_image_brightness_scale_paint} */

#endif
