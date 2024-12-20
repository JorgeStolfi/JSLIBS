/* uint16_image_RGB_table.h - hash table of RGB colors */
/* Last edited on 2024-12-04 23:38:06 by stolfi */ 

#ifndef uint16_image_RGB_table_H
#define uint16_image_RGB_table_H

#include <jspnm.h>
#include <uint16_image_RGB_hist.h>

typedef struct uint16_image_RGB_bucket_item *uint16_image_RGB_bucket;
typedef struct uint16_image_RGB_bucket_item { struct uint16_image_RGB_hist_item ch; uint16_image_RGB_bucket next; } uint16_image_RGB_bucket_item;

typedef uint16_image_RGB_bucket *uint16_image_RGB_table;

uint16_image_RGB_table uint16_image_RGB_table_build
  ( uint16_t **samples,
    int32_t chns, 
    int32_t cols, 
    int32_t rows, 
    int32_t maxcolors, 
    int32_t *colorsP );

int32_t uint16_image_RGB_table_lookup (uint16_image_RGB_table cht, ppm_pixel_t *colorP);

uint16_image_RGB_hist_vector uint16_image_RGB_table_to_hist (uint16_image_RGB_table cht, int32_t maxcolors);
uint16_image_RGB_table uint16_image_RGB_hist_to_table (uint16_image_RGB_hist_vector chv, int32_t colors);

int32_t uint16_image_RGB_table_add (uint16_image_RGB_table cht, ppm_pixel_t *colorP, int32_t value);
  /* Returns -1 on failure. */

uint16_image_RGB_table uint16_image_RGB_table_alloc (void);

void uint16_image_RGB_table_free (uint16_image_RGB_table cht);

#endif
