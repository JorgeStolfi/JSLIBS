/* uint16_image_quantize_floyd.h - image quantization by the Floyd-Steinberg algorithm. */
/* Last edited on 2024-12-26 19:51:22 by stolfi */ 

#ifndef uint16_image_quantize_floyd_H
#define uint16_image_quantize_floyd_H

#include <stdint.h>

#include <jspnm.h>
#include <uint16_image.h>
#include <uint16_image_match.h>
#include <uint16_color_table.h>
#include <uint16_color_tree.h>

void uint16_image_quantize_floyd
  ( uint16_image_t *img, 
    uint16_color_table_t *chv, 
    uint16_t mapmaxval,
    uint32_t maxColors,
    uint16_image_match_proc_t *match
  );
  /* Replaces each pixel in the image {img} by a pixel from
    {chv[0..NH-1]} where {NH = chv.ne}. The difference between the two
    pixels is propagated to the adjacent pixels that have not been
    processed yet.
    
    Assumes that the {R}, {G}, and {B} samples in {chv} are in the range
    {0..mapmaxval}, represening the linear scale range {[0 _ 1]}. Note
    that {mapmaxval} may be distinct from the initial {img.maxval}, so
    the search for similar colors implies a scaling from {0..img.maxval}
    to {0..mapmaxval}. On exit, {img.maxval} will be {mapmaxval}.
    
    While {NH} is less than {maxColors}, may append to {chv} more
    (scaled) colors that are encoutered in the scan, provided that their
    coordinates are in the range {0..mapmaxval}. On exit, {chv.ne} will
    be the number of colors actually used. */

#endif
