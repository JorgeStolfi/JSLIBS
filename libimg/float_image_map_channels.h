#ifndef float_image_map_channels_H
#define float_image_map_channels_H

/* Applying color transformations to images. */
/* Last edited on 2023-01-07 13:17:03 by stolfi */ 

#include <bool.h>
#include <ix.h>
#include <r2.h>
#include <r2x2.h>
#include <r3x3.h>
#include <interval.h>
#include <float_image.h>

typedef void float_image_map_channels_proc_t (int32_t NCA, float vA[], int32_t NCB, float vB[]);
  /* Type of a procedure that maps sample values {vA[0..NCA-1]} to sample values {vB[0..NCB-1]}. */
  
void float_image_map_channels
  ( float_image_t *imgA, 
    float_image_t *imgB, 
    float_image_map_channels_proc_t map 
  );
  /* Sets each pixel of image {imgB} to the corresponding pixel of {imgA},
    mapped by the function {map}.  The two images must have the same number
    of rows and columns.  */

void float_image_map_channels_RGB_to_YUV
  ( float_image_t *imgA, 
    float_image_t *imgB
  );
  /* The two images must have the same column and row count, {imgA} must have at least 3 channels.
    Maps the first 3 channels of {imgA} from RGB space to  European TV YUV coordinates
    (see {frgb_to_YUV}) and stores in the first 3 channels of {imgB}. The remaining channels
    of {imgA} are copied unchanged.
    
    If {imgB} has fewer channels than {imgA}, excess channels of the
    result are not stored. If it has more channels, excess chanels of
    {imgB} are filled with zeros. */

#endif
