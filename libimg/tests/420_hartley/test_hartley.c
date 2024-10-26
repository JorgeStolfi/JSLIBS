#define PROG_NAME "test_hartley"
#define PROG_DESC "test of {float_image_transform.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-10-12 17:38:17 by stolfi */ 
/* Created on 2008-09-21 by J. Stolfi, UNICAMP */

#define test_hartley_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <jsfile.h>
#include <bool.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <float_image_hartley.h>
#include <float_image_to_uint16_image.h>
#include <float_image.h>

int main(int argn, char **argv);

void do_test_transform_set
  ( float_image_t *iimg,
    int kx, 
    int ky,
    float_image_t *timg,
    float_image_t *oimg
  );
  /* Tests the direct and inverse transform on blip-image inputs and
    wave-image inputs with indices {±kx,±ky}, {±kx,0}, {±0,ky}, skipping
    duplicated tests. */

void do_test_transform
  ( float_image_t *iimg,
    int kx, 
    int ky,
    bool_t wave,
    float_image_t *timg,
    float_image_t *oimg
  );
  /* Tests the direct and inverse transform on a blip image with
    position {kx,ky} (if {wave=FALSE}) or a wave image with
    frequencies {kx,ky} (if {wave=TRUE}).  The parameters {kx,ky}
    will be automatically reduced to the image's domain. */
  
void do_test_basis(float_image_t *oimg, int nx, int ny);
  /* Generates images of all Hartley wave basis elements with frequency vectors {(kx,ky)} in 
    {{-nx..+nx} × {-ny..+ny}}.  The images are called "out/elem_{kx}_{ky}.png". */
    
void do_test_basis_element
  ( float_image_t *oimg,
    int kx, 
    int ky
  );
  /* Generates an image the Hartley wave basis element with frequency vector {(kx,ky)}.
    The image will be called "out/elem_{kx}_{ky}.png". */
    
void write_image(int kx, int ky, bool_t wave, char *prefix, char *suffix, float_image_t *img);

int main (int argn, char **argv)
  {
    /* Choose image size: */
    int chns = 3;
    int cols = atoi(argv[1]);
    int rows = atoi(argv[2]);
    
    /* Input image for testing: */
    float_image_t *iimg = float_image_new(chns, cols, rows);
    
    /* Transform image: */
    float_image_t *timg = float_image_new(chns, cols, rows);

    /* Output image: */
    float_image_t *oimg = float_image_new(chns, cols, rows);

    /* Tests with constant images: */
    do_test_transform_set(iimg, 0, 0, timg, oimg);

    /* Tests with lowest-freq waves: */
    do_test_transform_set(iimg, 1, 1, timg, oimg);

    /* Tests with mid-freq waves: */
    do_test_transform_set(iimg, cols/3, rows/3, timg, oimg);

    /* Highest possible frequencies: */
    do_test_transform_set(iimg, cols/2, rows/2, timg, oimg);
    
    /* Generate a full wave basis: */
    int nx = 3;
    int ny = 3;
    do_test_basis(oimg, nx, ny);

    return 0;
  }

void do_test_transform_set
  ( float_image_t *iimg,
    int kx, 
    int ky,
    float_image_t *timg,
    float_image_t *oimg
  )
  {
    do_test_transform(iimg, +kx,  00, FALSE, timg, oimg);
    do_test_transform(iimg, +kx,  00, TRUE,  timg, oimg);
    
    if (kx != 0)
      { do_test_transform(iimg, -kx,  00, FALSE, timg, oimg);
        do_test_transform(iimg, -kx,  00, TRUE,  timg, oimg);
      }
    
    if (ky != 0)
      { do_test_transform(iimg, 00,  +ky, FALSE, timg, oimg);
        do_test_transform(iimg, 00,  +ky, TRUE,  timg, oimg);

        do_test_transform(iimg, 00,  -ky, FALSE, timg, oimg);
        do_test_transform(iimg, 00,  -ky, TRUE,  timg, oimg);
      }
    
    if ((kx != 0) && (ky != 0))
      { 
        do_test_transform(iimg, +kx, +ky, FALSE, timg, oimg);
        do_test_transform(iimg, +kx, +ky, TRUE,  timg, oimg);

        do_test_transform(iimg, +kx, -ky, FALSE, timg, oimg);
        do_test_transform(iimg, +kx, -ky, TRUE,  timg, oimg);

        do_test_transform(iimg, -kx, +ky, FALSE, timg, oimg);
        do_test_transform(iimg, -kx, +ky, TRUE,  timg, oimg);

        do_test_transform(iimg, -kx, -ky, FALSE, timg, oimg);
        do_test_transform(iimg, -kx, -ky, TRUE,  timg, oimg);
      }
  }

void do_test_transform
  ( float_image_t *iimg,
    int kx, 
    int ky,
    bool_t wave,
    float_image_t *timg,
    float_image_t *oimg
  )
  {
    int chns = (int)iimg->sz[0];
    int cols = (int)iimg->sz[1];
    int rows = (int)iimg->sz[2];
    int c;
    
    fprintf(stderr, "------------------------------------------------------------\n");
    fprintf(stderr, "testing with %s image \n", (wave ? "wave" : "blip")); 
    fprintf(stderr, "sizes = %d %d  indices = %d %d\n", cols, rows, kx, ky); 
    
    /* Reduce frequencies to the natural range: */
    int rkx = kx % cols; if (rkx < 0) { rkx += cols; }
    int rky = ky % rows; if (rky < 0) { rky += rows; }
    
    fprintf(stderr, "creating input image %d %d %c...\n", kx, ky, "FT"[wave]);
    if (wave)
      { double amp = sqrt(2.0/(cols*rows)); /* So that the channel energy is 1. */
        float_image_hartley_wave(iimg, rkx, rky, amp);
      }
    else
      { float_image_fill(iimg, 0.0);
        double amp = 1.0; /* so that the total channel energy is 1. */
        float_image_fill_pixel(iimg, rkx, rky, (float)amp);
      }
      
    write_image(kx, ky, wave, "out/tr", "0-orig", iimg);
    
    double ierg = 0.0;
    for (c = 0; c < chns; c++) 
      { ierg += float_image_compute_squared_sample_sum(iimg, c, 0.0, NULL); }
    fprintf(stderr, "input image energy = %24.16e\n", ierg);
    demand(fabs(ierg/chns - 1.0) < 1.0e-6, "input image does not have unit power per channel!"); 
    
    float foo = (float)(1/M_PI); /* An arbitrary value. */
    fprintf(stderr, "filling transform image with %8.6f ...\n", foo);
    float_image_fill(timg, foo);
    
    fprintf(stderr, "applying transform...\n");
    float_image_hartley_transform(iimg, timg);
    write_image(kx, ky, wave, "out/tr", "1-hart", timg);
    
    double terg = 0.0;
    for (c = 0; c < chns; c++)
      { terg += float_image_compute_squared_sample_sum(timg, c, 0.0, NULL); }
    double tmag = sqrt(terg/ierg);
    fprintf(stderr, "transform energy = %24.16e  magnification = %24.16e\n", terg, tmag);
    demand(fabs(tmag - 1.0) < 1.0e-6, "transform did not preserve power!"); 
    
    float bar = (float)(1/M_SQRT2); /* An arbitrary value. */
    fprintf(stderr, "filling output image with %8.6f ...\n", bar);
    float_image_fill(oimg, bar);
    
    fprintf(stderr, "applying inverse transform...\n");
    float_image_hartley_transform(timg, oimg);
    write_image(kx, ky, wave, "out/tr", "2-hinv", oimg);
    
    double oerg = 0.0;
    for (c = 0; c < chns; c++) 
      { oerg += float_image_compute_squared_sample_sum(oimg, c, 0.0, NULL); }
    double omag = sqrt(oerg/terg);
    fprintf(stderr, "transform energy = %24.16e  magnification = %24.16e\n", oerg, omag);
    demand(fabs(omag - 1.0) < 1.0e-6, "inverse transform did not preserve power!"); 
    
    fprintf(stderr, "------------------------------------------------------------\n");
  }  

void do_test_basis(float_image_t *oimg, int nx, int ny)
  {
    
    for (int ky = -ny; ky <= ny; ky++)
      { for (int kx = -nx; kx <= nx; kx++)
          { do_test_basis_element(oimg, kx, ky); }
      }
  }

void do_test_basis_element
  ( float_image_t *oimg,
    int kx, 
    int ky
  )
  {
    int cols = (int)oimg->sz[1];
    int rows = (int)oimg->sz[2];
    
    fprintf(stderr, "------------------------------------------------------------\n");
    fprintf(stderr, "geenrating Hartley basis element image\n"); 
    fprintf(stderr, "sizes = %d %d  indices = %d %d\n", cols, rows, kx, ky); 
    
    /* Reduce frequencies to the natural range: */
    int rkx = kx % cols; if (rkx < 0) { rkx += cols; }
    int rky = ky % rows; if (rky < 0) { rky += rows; }
    
    fprintf(stderr, "creating element %+d %+d...\n", kx, ky);
    double amp = 1; /* For display. */
    float_image_hartley_wave(oimg, rkx, rky, amp);
      
    write_image(kx, ky, TRUE, "out/elem", "wave", oimg);
    fprintf(stderr, "------------------------------------------------------------\n");
  }  

void write_image(int kx, int ky, bool_t wave, char *prefix, char *suffix, float_image_t *img)
  {
    char *fname = NULL;
    asprintf(&fname, "%s-%+04d-%+04d-%c-%s.ppm", prefix, kx, ky, "FT"[wave], suffix);
    FILE *wr = open_write(fname, TRUE);
    int chns = (int)img->sz[0];
    /* Map true min to 0, true max to {fim->maxval}: */
    double vLo[chns];
    double vHi[chns];
    int c;
    for (c = 0; c < chns; c++) 
      { float vMin = +INF, vMax = -INF;
        float_image_update_sample_range(img, c, &vMin, &vMax);
        vLo[c] = vMin; vHi[c] = vMax;
      }
    bool_t yup = TRUE, verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    uint16_image_t *pimg = float_image_to_uint16_image(img, isMask, chns, vLo, vHi, NULL, 255, yup, verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pimg, forceplain, verbose);
    uint16_image_free(pimg);
    fclose(wr);
    free(fname);
  }
