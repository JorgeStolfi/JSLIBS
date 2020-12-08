/* uint16_image_RGB_table.h - hash table of RGB colors */
/* Last edited on 2009-01-07 01:40:57 by stolfi */ 

#ifndef uint16_image_RGB_table_H
#define uint16_image_RGB_table_H

#include <jspnm.h>
#include <uint16_image_RGB_hist.h>

typedef struct uint16_image_RGB_bucket_item *uint16_image_RGB_bucket;
typedef struct uint16_image_RGB_bucket_item { struct uint16_image_RGB_hist_item ch; uint16_image_RGB_bucket next; } uint16_image_RGB_bucket_item;

typedef uint16_image_RGB_bucket *uint16_image_RGB_table;

uint16_image_RGB_table uint16_image_RGB_table_build
  ( uint16_t **samples,
    int chns, 
    int cols, 
    int rows, 
    int maxcolors, 
    int *colorsP );

int uint16_image_RGB_table_lookup (uint16_image_RGB_table cht, ppm_pixel_t *colorP);

uint16_image_RGB_hist_vector uint16_image_RGB_table_to_hist (uint16_image_RGB_table cht, int maxcolors);
uint16_image_RGB_table uint16_image_RGB_hist_to_table (uint16_image_RGB_hist_vector chv, int colors);

int uint16_image_RGB_table_add (uint16_image_RGB_table cht, ppm_pixel_t *colorP, int value);
  /* Returns -1 on failure. */

uint16_image_RGB_table uint16_image_RGB_table_alloc (void);

void uint16_image_RGB_table_free (uint16_image_RGB_table cht);

#endif
