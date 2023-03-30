/* uint16_image.h - simple in-memory format for image-like arrays of {uint16_t}. */
/* Last edited on 2023-03-17 20:45:23 by stolfi */

#ifndef uint16_image_H
#define uint16_image_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <jspnm.h>

#define uint16_image_MAX_SAMPLE (65535u)
  /* The maximum sample value that can be stored in a {uint16_image_t}. */

/* Beware that the order of parameters is inconsitent: sometimes {chn,col,row},
   sometimes {col,row,chn}. */

typedef struct uint16_image_t
  { uint16_t maxval;
    int32_t cols;
    int32_t rows;
    int32_t chns;            /* Normally either 1 or 3. */
    uint16_t **smp;  /* Rows of samples. */
  } uint16_image_t;
  /* A PGM ({chns=1}) or PPM ({chns=3}) image in memory.
    May also be used for other image formats with other
    values of {chns}. The sample in channel {c} of the
    pixel in row {y} and column {x} is {smp[y][x*chns + c]}. */

uint16_image_t *uint16_image_new(int32_t cols, int32_t rows, int32_t chns);
  /* Alocates a new PPM image, including the header and 
    (if {rows*cols != 0}) the pixel arrays.  Sets the 
    {maxval} field to 0 (an invalid value); clients
    must set it explicitly. */
    
void uint16_image_free(uint16_image_t *img);
  /* Discards the pixel array and header of a PNM image. */

uint16_t** uint16_image_alloc_pixel_array(int32_t cols, int32_t rows, int32_t chns);
void uint16_image_free_pixel_array(uint16_t **smp, int32_t cols, int32_t rows, int32_t chns);
  /* Alocates/discards the sample array for a PNM image of the
    specified size and number of channels. */
    
uint16_t* uint16_image_alloc_pixel_row(int32_t cols, int32_t chns);
void uint16_image_free_pixel_row(uint16_t *smp, int32_t cols, int32_t chns);
  /* Alocates/discards a row of pixels for a PNM image of the
    specified width and number of channels. */
    
uint16_image_t *uint16_image_new_header(int32_t cols, int32_t rows, int32_t chns);
  /* Creates a new header record for a image with the specified 
    size and number of channels, but leaves the pixel arrays as NULL
    and the {maxval} to 0 (an invalid value); clients
    must set it explicitly. */

uint16_t uint16_image_get_sample(uint16_image_t *img, int32_t c, int32_t x, int32_t y);
  /* Returns the sample value in channel {c}, column {x}, row {y} 
    of the image {img}. Non-existent channels
    are assumed to be all zeros. */

void uint16_image_set_sample(uint16_image_t *img, int32_t c, int32_t x, int32_t y, uint16_t pv);
  /* Stores the sample value {pv} in channel {c}, column {x}, row {y}
    of the image {img}. Fails if {x,y,c} are not a valid. */

uint16_image_t *uint16_image_crop(uint16_image_t *img, int32_t ic, int32_t nc, int32_t ix, int32_t nx, int32_t iy, int32_t ny);
  /* Returns the sub-image of {img} consisting of channels {ic..ic+nc-1},
     columns {ix..ix+nx-1}, and rows {iy..iy+ny-1}. */

void uint16_image_describe(FILE *wr, char *name, uint16_image_t *img);
  /* Prints the {name} and other parameters of image {img} 
    to stderr. */

#endif
