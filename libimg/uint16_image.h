/* uint16_image.h - simple in-memory format for image-like arrays of {uint16_t}. */
/* Last edited on 2024-12-26 12:44:28 by stolfi */

#ifndef uint16_image_H
#define uint16_image_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

#define uint16_image_MAX_SAMPLE (65535u)
  /* The maximum sample value that can be stored in a {uint16_image_t}. */

/* Beware that the order of parameters is inconsitent: sometimes {chn,col,row},
   sometimes {col,row,chn}. */

typedef struct uint16_image_t
  { uint16_t maxval;  /* Maximum sample value. */
    uint32_t cols;    /* Comput of pixel columns. */
    uint32_t rows;    /* Count of pixel rows. */
    uint32_t chns;    /* Count of samples per pixel. */
    uint16_t **smp;   /* The rows of samples. */
  } uint16_image_t;
  /* A discrete 2D multi-channel image in memory.  Could be
    the result of reading a PGM or PBM ({chns=1}) or PPM ({chns=3}) image file.
    May also be used for other image formats with other
    values of {chns}. The sample in channel {c} of the
    pixel in row {y} and column {x} is {smp[y][x*chns + c]}. */

uint16_image_t *uint16_image_new(uint32_t cols, uint32_t rows, uint32_t chns);
  /* Alocates a new image, including the header and 
    (if {rows*cols != 0}) the pixel arrays.  Sets the 
    {maxval} field to 0 (an invalid value); clients
    must set it explicitly. */
    
void uint16_image_free(uint16_image_t *img);
  /* Discards the pixel arrays and header of the image {img}. */

uint16_t** uint16_image_alloc_pixel_array(uint32_t cols, uint32_t rows, uint32_t chns);
void uint16_image_free_pixel_array(uint16_t **smp, uint32_t cols, uint32_t rows, uint32_t chns);
  /* Alocates/discards the sample array for an image of the
    specified size and number of channels. */
    
uint16_t* uint16_image_alloc_pixel_row(uint32_t cols, uint32_t chns);
void uint16_image_free_pixel_row(uint16_t *smp, uint32_t cols, uint32_t chns);
  /* Alocates/discards a row of pixels for an image of the
    specified width and number of channels. */
    
uint16_image_t *uint16_image_new_header(uint32_t cols, uint32_t rows, uint32_t chns);
  /* Creates a new header record for a image with the specified 
    size and number of channels, but leaves the pixel arrays as NULL
    and the {maxval} as 0 (an invalid value); clients
    must set both explicitly. */

uint16_t uint16_image_get_sample(uint16_image_t *img, int32_t c, int32_t x, int32_t y);
  /* Returns the sample value in channel {c}, column {x}, row {y} 
    of the image {img}. Non-existent channels
    are assumed to be all zeros. */

void uint16_image_set_sample(uint16_image_t *img, int32_t c, int32_t x, int32_t y, uint16_t smp);
  /* Stores the sample value {smp} in channel {c}, column {x}, row {y}
    of the image {img}. Fails if {x,y,c} are not a valid index combination
    of {smp} is not in {0..img->maxval}. */

bool_t uint16_image_same_color(uint32_t chns, uint16_t clra[], uint16_t clrb[]);
  /* Returns true iff {clra[c]==clrb[c]} for all {c} in {0..chns-1}. */

uint16_image_t *uint16_image_crop(uint16_image_t *img, uint32_t ic, uint32_t nc, uint32_t ix, uint32_t nx, uint32_t iy, uint32_t ny);
  /* Returns the sub-image of {img} consisting of channels {ic..ic+nc-1},
     columns {ix..ix+nx-1}, and rows {iy..iy+ny-1}. */

void uint16_image_describe(FILE *wr, char *name, uint16_image_t *img);
  /* Prints the {name} and other parameters of image {img} 
    to stderr. */

#endif
