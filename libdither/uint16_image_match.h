/* uint16_image_match.h - generic color matching. */
/* Last edited on 2024-12-26 19:34:21 by stolfi */ 

#ifndef uint16_image_match_H
#define uint16_image_match_H

#include <stdint.h>
#include <assert.h>

#include <jspnm.h>
#include <bool.h>
#include <uint16_color_table.h>

#define MAXCOLORS 32767

typedef int32_t uint16_image_match_proc_t
  ( int32_t R,
    int32_t G,
    int32_t B,
    uint16_t maxval,
    uint16_color_table_t *chv,
    uint16_t mapmaxval
  );
  /* A function that finds the index into the color histogram {chv} that
    best matches {(R,G,B)}. Assumes that {R,G,B} are intensities
    linearly mapped from {[0_1]} to {0..maxval}, while the colors in
    {chv} are mapped from {[0 _ 1]} to {0..mapmaxval}. 
    
    Does a simple linear search.
    
    To be useful in algorithms like Floyd-Steingerg, the procedure must work even
    when {R,G,B} are somewhat outside the range {0..maxval}; at least in
    the range {-maxval .. 2*maxval}. */

int32_t uint16_image_match_rgb 
  ( int32_t r, int32_t g, int32_t b, 
    uint16_t maxval, 
    uint16_color_table_t *chv, 
    uint16_t mapmaxval
  );
  /* Search colormap for closest match in RGB metric. */
  
int32_t uint16_image_match_yuv
  ( int32_t r, int32_t g, int32_t b, 
    uint16_t maxval, 
    uint16_color_table_t *chv, 
    uint16_t mapmaxval
  );
  /* Search colormap for closest match in YUV metric. */

#endif
