/* See uint16_image.h */
/* Last edited on 2023-03-17 20:46:42 by stolfi */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include <affirm.h>
#include <jsfile.h>
#include <jspnm.h>

#include <uint16_image.h>

uint16_image_t *uint16_image_new(int32_t cols, int32_t rows, int32_t chns)
  { /* Allocate header: */
    uint16_image_t *img = uint16_image_new_header(cols, rows, chns);
    /* Allocate pixels: */
    if ((cols != 0) && (rows != 0))
      { img->smp = uint16_image_alloc_pixel_array(cols, rows, chns); }
    return(img);
  }
  
void uint16_image_free(uint16_image_t *img)
  { if (img == NULL) return;
    if (img->smp != NULL) 
      { uint16_image_free_pixel_array(img->smp, img->cols, img->rows, img->chns); }
    free(img);
  }

uint16_image_t *uint16_image_new_header(int32_t cols, int32_t rows, int32_t chns)
  { /* Allocate header: */
    uint16_image_t *img = (uint16_image_t *)pnm_malloc(sizeof(uint16_image_t));
    assert(img != NULL);
    /* Initialize fields: */
    img->cols = cols;
    img->rows = rows;
    img->chns = chns;
    /* Leave sample arrays as NULL: */
    img->smp = (uint16_t**)NULL;
    /* Initialize the {maxval} field as 0 (invalid): */
    img->maxval = 0;
    return(img);
  }
  
uint16_t** uint16_image_alloc_pixel_array(int32_t cols, int32_t rows, int32_t chns)
  { int32_t row;
    uint16_t **smp = (uint16_t **)pnm_malloc(rows*sizeof(uint16_t *));
    for (row = 0; row < rows; row++)
      { smp[row] = uint16_image_alloc_pixel_row(cols, chns); }
    return smp;
  }
  
void uint16_image_free_pixel_array(uint16_t **smp, int32_t cols, int32_t rows, int32_t chns)
  { if (smp == NULL) { return; }
    for (int32_t row = 0; row < rows; row++) 
      { uint16_image_free_pixel_row(smp[row], cols, chns); }
    free(smp);
  }

uint16_t* uint16_image_alloc_pixel_row(int32_t cols, int32_t chns)
  { if (cols == 0)
      { return (uint16_t *)NULL; }
    else
      { return (uint16_t *)pnm_malloc(cols*chns*sizeof(uint16_t)); }
  }
  
void uint16_image_free_pixel_row(uint16_t *row, int32_t cols, int32_t chns)
  { 
    if (row != NULL) { free(row); }
  }

uint16_t uint16_image_get_sample(uint16_image_t *img, int32_t c, int32_t x, int32_t y)
  { if ((c < 0) || (c >= img->chns)) { return 0; }
    if ((x < 0) || (x >= img->cols)) { return 0; }
    if ((y < 0) || (y >= img->rows)) { return 0; }
    return img->smp[y][x*img->chns + c];
  }

void uint16_image_set_sample(uint16_image_t *img, int32_t c, int32_t x, int32_t y, uint16_t pv)
  { if ((c < 0) || (c >= img->chns)) { pnm_error("invalid channel"); }
    if ((x < 0) || (x >= img->cols)) { pnm_error("invalid column"); }
    if ((y < 0) || (y >= img->rows)) { pnm_error("invalid row"); }
    if (pv > img->maxval) { pnm_error("invalid sample value"); }
    img->smp[y][x*img->chns + c] = pv;
  }

uint16_image_t *uint16_image_crop(uint16_image_t *img, int32_t ic, int32_t nc, int32_t ix, int32_t nx, int32_t iy, int32_t ny)
  { demand((ic >= 0) && (ic+nc <= img->chns), "invalid channel range");
    demand((ix >= 0) && (ix+nx <= img->cols), "invalid column range");
    demand((iy >= 0) && (iy+ny <= img->rows), "invalid row range");
    uint16_image_t *omg = uint16_image_new(nx, ny, nc);
    omg->maxval = img->maxval;
    int32_t chn, col, row;
    for (row = 0; row < ny; row++)
      { for (col = 0; col < nx; col++)
          { uint16_t *pi = &(img->smp[iy + row][(ix + col)*img->chns + ic]);
            uint16_t *po = &(omg->smp[row][col*nc]);
            for (chn = 0; chn < nc; chn++) { (*po) = (*pi); po++; pi++; }
          }
      }
    return omg;
  }        

void uint16_image_describe(FILE *wr, char *name, uint16_image_t *img)
  {
    fprintf(stderr, "image %s\n", name);
    fprintf(stderr, "\n");
    fprintf(stderr, "cols =   %5d\n", img->cols);
    fprintf(stderr, "rows =   %5d\n", img->rows);
    fprintf(stderr, "chns =   %5d\n", img->chns);
    fprintf(stderr, "maxval = %5d\n", img->maxval);
    /* Determine min and max pixel value: */
    uint16_t minsmp = img->maxval;
    uint16_t maxsmp = 0;
    for (int32_t y = 0; y < img->rows; y++)
      { uint16_t *ip = img->smp[y];
        for (int32_t x = 0; x < img->cols; x++)
          { for (int32_t c = 0; c < img->chns; c++)
              { uint16_t smp = (*ip);
                assert(smp <= img->maxval);
                if (smp < minsmp) { minsmp = smp; }
                if (smp > maxsmp) { maxsmp = smp; }
                ip++;
              }
          }
      }
    fprintf(wr, "minsmp = %5d\n", minsmp);
    fprintf(wr, "maxsmp = %5d\n", maxsmp);
  }
