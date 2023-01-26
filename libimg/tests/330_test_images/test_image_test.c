#define PROG_NAME "test_image_test"
#define PROG_DESC "test of {float_image_test.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-01-14 20:29:17 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_image_test_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <ix.h>
#include <jsfile.h>
#include <bool.h>
#include <r2x2.h>
#include <r2.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <float_image_test.h>
#include <float_image.h>
#include <float_image_waves.h>
#include <float_image_write_pnm.h>

int main(int argn, char **argv);

void do_test_gen(float_image_test_generator_t *gen_proc, char *gen_name, char *outPrefix);

void write_color_image(float_image_t *img, char *outPrefix, char *tag);

int main (int argn, char **argv)
  {
    char *outPrefix = "out/test";
    
    do_test_gen(&float_image_test_gen_stripes, "stripes", outPrefix);
    do_test_gen(&float_image_test_gen_ripples, "ripples", outPrefix);
    do_test_gen(&float_image_test_gen_checker, "checker", outPrefix);
    do_test_gen(&float_image_test_gen_chopsea, "chopsea", outPrefix);
    
    return 0;
  }

void do_test_gen(float_image_test_generator_t *gen_proc, char *gen_name, char *outPrefix)
  {
    fprintf(stderr, "=== test gen=%s============\n", gen_name);

    fprintf(stderr, "creating image \"%s\" ...\n", gen_name);
    int NC = 3;
    int NX = 1024;
    int NY = 768;
    float_image_t *img = float_image_new(NC, NX, NY);
    
    fprintf(stderr, "filling image \"%s\" ...\n", gen_name);
    float_image_test_paint(img, gen_proc, 1);
    
    /* Write the input and output images: */
    write_color_image(img, outPrefix, gen_name);
    fprintf(stderr, "===================================================\n");
  }  

void write_color_image(float_image_t *img, char *outPrefix, char *tag)
  {
    assert(img->sz[0] == 3);

    char *fname = NULL;
    asprintf(&fname, "%s-%s.ppm", outPrefix, tag);
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    double gamma = 1.000;
    double bias = 0.000;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = TRUE;
    float_image_write_pnm_named(fname, img, isMask, gamma, bias, yup, warn, verbose);
    free(fname);
  }
