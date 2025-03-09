#ifndef pst_map_clear_H
#define pst_map_clear_H

/* pst_map_clear.h -- procedures for comparing maps (height, slopes, etc). */
/* Last edited on 2025-03-06 12:59:32 by stolfi */

#include <stdint.h>

#include <bool.h>
#include <float_image.h>

void pst_map_clear
  ( float_image_t *A,
    int32_t wch,
    int32_t xmin,
    int32_t xmax,
    int32_t ymin,
    int32_t ymax
  );
  /* For each channel {c!=wch}, sets to {NAN} every sample {A[c,x,y]}
    with {x} in {xmin..xmax} and {y} in {ymin..ymax}. If {wch} is a
    valid channel index, also sets {A[wch,x,y]} to zero. */
  
#endif
