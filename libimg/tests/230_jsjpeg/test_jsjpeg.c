/* Test of jsjpeg.h, uint16_image_io_jpeg.h */
/* Last edited on 2024-12-26 16:46:48 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <jspnm.h>
#include <jsprintf.h>
#include <uint16_image.h>
#include <uint16_image_read_jpeg.h>
#include <uint16_image_write_jpeg.h>


int32_t main (int32_t argc, char **argv);

void do_uint16_image_read_write_jpeg_tests(char *prefix);
  /* Performs various image I/O tests, reading files
    with names "{prefix}-rd-{VAR}.jpg" and writing files
    with names "{prefix}-wr-{QQQ}-{VAR}.jpg", where 
    {VAR} is "gry" or "rgb" and {QQQ} is the 3-digit 
    JPEG quality parameter. */

void do_uint16_image_read_write_jpeg_test(char *prefix, int32_t ik, int32_t iq);
  /* Reads an image from file "{prefix}-rd-{VAR}.jpg" and writes the
    complemented image as files  with names "{prefix}-wr-{QQQ}-{VAR}.jpg".

    The {ik} paramter indicates the colorspace. Currently the only 
    valid choices are 0 (the space {JCS_GRAYSCALE} of {jpeglib.h},
    in which case {VAR} is "GRY") or 1 ({JCS_RGB}in which case {VAR}
    is "RGB").
    
    The {iq} parameter determines the JPEG {quality} setting. When
    {iq} is 0 the quality is 100% (no compression). As {iq} grows the
    quality decreases (non-linearly). The {QQQ} field in utput
    filenames is the 3-digit decimal value of the JPEG {quality}. */
    
void frobnicate_image(uint16_image_t *img, uint16_image_t *omg);  
  /* Copies {img} ino {omg}, complementing some pixels
    relative to {img.maxval}, flipping left-right and
    top-bottom. Also sets {omg.maxval = img.maxval}. */

int32_t main (int32_t argc, char **argv)
  { do_uint16_image_read_write_jpeg_tests("out/test");
    return 0;
  }
  
void do_uint16_image_read_write_jpeg_tests(char *prefix)
  { int32_t iq, ik;
    for (iq = 0; iq <= 6; iq ++)
      for (ik = 0; ik < 2; ik++)
        { do_uint16_image_read_write_jpeg_test(prefix, ik, iq); }
  }
      
void do_uint16_image_read_write_jpeg_test(char *prefix, int32_t ik, int32_t iq)
  { 
    /* Tag indicating format variant: */
    /* Filename extension: */
    char *kindtbl[2] = {"GRY", "RGB"};
    char *xkind = kindtbl[ik]; 
    
    /* JPEG quality parameter: */
    int32_t quality = 100 - iq*iq;
    
    fprintf(stderr, "-----------------------------------------\n");
    fprintf(stderr, "testing prefix = %s kind = %s quality = %d\n", prefix, xkind, quality);
    fprintf(stderr, "\n");

    char *fname = jsprintf("%s-rd-%s.jpg", prefix, xkind);
    int32_t kind;
    uint16_image_t *img = uint16_image_read_jpeg_named(fname, TRUE, &kind);
    uint16_image_describe(stderr, fname, img);
    fprintf(stderr, "JPEG colorspace = %d\n", kind);
    free(fname);

    uint16_image_t *omg = uint16_image_new(img->cols, img->rows, img->chns);
    frobnicate_image(img, omg);
    fprintf(stderr, "\n");

    fname = jsprintf("%s-wr-%03d-%s.jpg", prefix, quality, xkind);
    uint16_image_describe(stderr, fname, omg);
    uint16_image_write_jpeg_named(fname, omg, quality, TRUE);
    free(fname);
    fprintf(stderr, "-----------------------------------------\n");
  }

void frobnicate_image(uint16_image_t *img, uint16_image_t *omg)
  {
    uint32_t cols = img->cols; assert(img->cols == omg->cols);
    uint32_t rows = img->rows; assert(img->rows == omg->rows);
    uint32_t chns = img->chns; assert(img->chns == omg->chns);

    omg->maxval = img->maxval;
    
    int32_t skip = 20;

    int32_t x, y, c;
    for (y = 0; y < rows; y++)
      { bool_t yrev = ((y >= skip) && (y < (int32_t)rows-skip));
        uint16_t *ip = img->smp[y];
        uint16_t *op = omg->smp[(int32_t)rows - 1 - y];
        int32_t ik = 0;
        int32_t ok = (int32_t)(cols*chns - 1);
        for (x = 0; x < cols; x++)
          { bool_t xrev = ((x >= skip) && (x < (int32_t)cols-skip));
            bool_t rev = (xrev && yrev);
            for (c = 0; c < chns; c++)
              { op[ok] = ( rev ? (uint16_t)(img->maxval - ip[ik]) : ip[ik]);
                ik++; ok--;
              }
          }
      }
  }
