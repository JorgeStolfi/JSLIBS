#define PROG_NAME "test_test_image"
#define PROG_DESC "test of {float_image_test.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-02-04 19:35:10 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_test_image_COPYRIGHT \
  "Copyright � 2007  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <ix.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <bool.h>
#include <r2x2.h>
#include <r2.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <float_image_test.h>
#include <float_image.h>
#include <float_image_waves.h>
#include <float_image_write_gen.h>

int main(int argn, char **argv);

void do_all_tests(int32_t NC, int32_t NX, int32_t NY);

void do_test_gen(int32_t NC, int32_t NX, int32_t NY, float_image_test_generator_t *gen_proc, char *gen_name, char *outPrefix);

void write_color_image(float_image_t *img, char *outPrefix, int32_t NC, int32_t NX, int32_t NY, char *funcName);

int main (int argn, char **argv)
  {
    do_all_tests(3,  320,  240);
    do_all_tests(1, 1024, 1024);
    
    return 0;
  }

void do_all_tests(int32_t NC, int32_t NX, int32_t NY)
  {
    char *outPrefix = "out/test";
    
    do_test_gen(NC, NX, NY, &float_image_test_gen_bullsex, "bullsex", outPrefix);
    do_test_gen(NC, NX, NY, &float_image_test_gen_bullsqr, "bullsqr", outPrefix);
    do_test_gen(NC, NX, NY, &float_image_test_gen_checker, "checker", outPrefix);
    do_test_gen(NC, NX, NY, &float_image_test_gen_chopsea, "chopsea", outPrefix);
    do_test_gen(NC, NX, NY, &float_image_test_gen_grittie, "grittie", outPrefix);
    do_test_gen(NC, NX, NY, &float_image_test_gen_ripples, "ripples", outPrefix);
    do_test_gen(NC, NX, NY, &float_image_test_gen_stripes, "stripes", outPrefix);
  }

void do_test_gen(int32_t NC, int32_t NX, int32_t NY, float_image_test_generator_t *gen_proc, char *gen_name, char *outPrefix)
  {
    fprintf(stderr, "=== test gen = %s============\n", gen_name);

    fprintf(stderr, "creating image \"%s\" ...\n", gen_name);
    float_image_t *img = float_image_new(NC, NX, NY);
    
    fprintf(stderr, "filling image \"%s\" ...\n", gen_name);
    float_image_test_paint(img, gen_proc, 3);
    
    /* Write the input and output images: */
    write_color_image(img, outPrefix, NC, NX, NY, gen_name);
    fprintf(stderr, "===================================================\n");
  }  

void write_color_image(float_image_t *img, char *outPrefix, int32_t NC, int32_t NX, int32_t NY, char *funcName)
  {
    char *fname = jsprintf("%s-%04dx%04d-%d-%s.png", outPrefix, NX, NY, NC, funcName);
    double expoEnc = 1.000;
    double bias = 0.000;
    bool_t verbose = TRUE;
    image_file_format_t ffmt = image_file_format_PNG;
    bool_t yUp = FALSE; /* !!! Should test with {TRUE} too !!! */
    float_image_write_gen_named(fname, img, ffmt, yUp, 0.0, 1.0, expoEnc, bias, verbose);
    free(fname);
  }
