#define PROG_NAME "test_geostereo"
#define PROG_DESC "test of {float_image_geostereo.h} and {float_image_geostereo_uniscale.h}"
#define PROG_VERS "1.0"

/* Last edited on 2017-06-26 04:20:01 by stolfilocal */ 
/* Created on 2009-06-02 by J. Stolfi, UNICAMP */

#define test_geostereo_COPYRIGHT \
  "Copyright © 2017  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <jsfile.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_read_pnm.h>
#include <float_image_write_pnm.h>
#include <float_image_geostereo.h>
#include <float_image_geostereo_uniscale.h>

int32_t main(int32_t argn, char **argv);

/* In all file names, the {ext} part. */ 

void do_geostereo_tests(char *name0, char *name1, int32_t NC, int32_t ncands);
  /* Applies {float_image_geostereo_uniscale} to images
    "in/{name1}.{ext}" and "in/{name2}.{ext}"; where {ext} is
    "pgm" if {NC} is 1, or "ppm" if it is 3. 
    
    Keeps {ncands} best matches for each pixel. Writes results as
    "out/disp-{ncands}.{ext}" and "out/score-{ncands}.{ext} where {ext}
    is "fni" always, "pgm" if {ncands} is 1, "ppm" if {ncands} is 3. */

float_image_t *get_test_image(char *name, int32_t NC);
  /* Reads a test image from "in/{name}.{ext}" where {ext} is "pgm"
    if {NC} is 1, "ppm" if it is 3. */ 

void write_image(float_image_t *img, char *name, double vmin, double vmax);
  /* Writes the image {img} to file "out/disp-{NC}.fni" where {NC} is
    the number of channels in the image. IF {NC} is 1, also writes
    "out/disp-{NC}.pgm". If {NC} is 3, also writes "out/disp-{NC}.ppm".
    In the PGM/PPM files, samples are mapped affinely from {vmin} to 0
    and {vmax} to 65535.  May modify {img}.  */
 
char *makefname(char *name, int32_t NC, char *ext);
  /* Returns the string "out/{name}-{NC}.{ext}". */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    demand(argc == 3, "wrong num of parameters");
    char *name0 = argv[1];
    char *name1 = argv[2];
    
    do_geostereo_tests(name0, name1, 1, 3);
    do_geostereo_tests(name0, name1, 3, 3);

    return 0;
  }

void do_geostereo_tests(char *name0, char *name1, int32_t NC, int32_t ncands)
  {
    fprintf(stderr, "testing with images %s %s\n", name0, name1);
    
    float_image_t *img0 = get_test_image(name0, NC);
    float_image_t *img1 = get_test_image(name1, NC);

    int32_t nwx = 3;  /* Window width. */
    int32_t nwy = 3;  /* Window height. */
    double dmin = -10.0; /* Min displacement to consider (pixels). */
    double dmax = +10.0; /* Max displacement to consider (pixels). */
    float_image_t *imgd; /* Displacement map. */
    float_image_t *imgs; /* Score map. */
    float_image_geostereo_uniscale
      ( img0, img1, 
        nwx, nwy, 
        dmin, dmax, 
        ncands,
        &imgd, &imgs
      );
    
    write_image(imgd, "disp", dmin, dmax);
    double smax = 1.0; /* Guessing. */
    write_image(imgs, "score", 0.0, smax);
    
    float_image_free(imgd);
    float_image_free(imgs);
    
    float_image_free(img0);
    float_image_free(img1);
  }   

float_image_t *get_test_image(char *name, int32_t NC)
  {
    demand((NC == 1) || (NC == 3), "bad num of channels");
    char *fname = NULL;
    asprintf(&fname, "in/%s.%s", name, (NC == 3 ? "ppm" : "pgm"));
    bool_t isMask = FALSE; /* Assume pixels have a smooth distribution. */
    bool_t yup = FALSE;
    bool_t warn = TRUE;
    bool_t verbose = FALSE;
    float_image_t *img = float_image_read_pnm_named(fname, isMask, 1.0000, 0.0327, yup, warn, verbose);
    free(fname);
    return img;
  }
  
void write_image(float_image_t *img, char *name, double vmin, double vmax)
  {
    int32_t NC = (int32_t)img->sz[0];
    bool_t isMask = FALSE; /* Assume pixels have a smooth distribution. */
    bool_t yup = FALSE; /* Reverse order of lines in file. */
    bool_t warn = TRUE; /* Print message on open. */
    bool_t verbose = FALSE;
    
    /* Write as FNI: */
    { char *fname = makefname(name, NC, "fni");
      FILE *wr = open_write(fname, warn);
      float_image_write(wr, img);
      fclose(wr);
      free(fname);
    }
          
    if ((NC == 1) || (NC == 3))
      { /* Write as PGM/PPM: */
        char *fname = makefname(name, NC, (NC == 3 ? "ppm" : "pgm"));
        for (int32_t c = 0; c < NC; c++)
          { float_image_rescale_samples(img, c, (float)vmin, (float)vmax, 0.0, 1.0); }
        float_image_write_pnm_named(fname, img, isMask, 1.0000, 0.0327, yup, warn, verbose);
        free(fname);
      }
  }

char *makefname(char *name, int32_t NC, char *ext)
  { 
    char *fname = NULL;
    asprintf(&fname, "out/%s-%02d.%s", name, NC, ext);
    return fname;
  }

