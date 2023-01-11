#define PROG_NAME "test_image_test"
#define PROG_DESC "test of {float_image_test.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-01-10 16:03:39 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_image_test_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <affirm.h>
#include <ix.h>
#include <jsfile.h>
#include <bool.h>
#include <r2x2.h>
#include <r2.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <float_image_test.h>
#include <float_image_to_uint16_image.h>
#include <float_image.h>

int main(int argn, char **argv);

void do_test_gen(float_image_test_generator_t *gen_proc, char *gen_name);

void do_test_comb_waves(int32_t NF0, int32_t NF1);

void write_image(char *gen_name, float_image_t *img);

int main (int argn, char **argv)
  {
    do_test_gen(&float_image_test_gen_stripes, "stripes");
    do_test_gen(&float_image_test_gen_ripples, "ripples");
    do_test_gen(&float_image_test_gen_checker, "checker");
    do_test_gen(&float_image_test_gen_chopsea, "chopsea");
    
    do_test_comb_waves(0, 0);
    do_test_comb_waves(1, 1);
    do_test_comb_waves(2, 2);
    do_test_comb_waves(0, 2);
    
    return 0;
  }

void do_test_gen(float_image_test_generator_t *gen_proc, char *gen_name)
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
    write_image(gen_name, img);
    fprintf(stderr, "===================================================\n");
  }  

void do_test_comb_waves(int32_t NF, int32_t NS)
  {
    fprintf(stderr, "=== test comb_waves = %d..%d ============\n", NF, NS);

    fprintf(stderr, "creating image...\n");
    int NC = 3;
    int NX = 1024;
    int NY = 768;
    float_image_t *img = float_image_new(NC, NX, NY);
    
    fprintf(stderr, "choosing frequencies...\n");
    double amp[NF];
    double fx[NF];
    double fy[NF]; 
    double phase[NF];
    
    fprintf(stderr, "filling image ...\n");
    float_image_test_paint(img, gen_proc, 1);
    
    /* Write the input and output images: */
    write_image(fname, img);
    fprintf(stderr, "===================================================\n");
  }  

void write_image(char *gen_name, float_image_t *img)
  {
    char *fname = NULL;
    asprintf(&fname, "out/test-%s.ppm", gen_name);
    FILE *wr = open_write(fname, TRUE);
    int chns = (int)img->sz[0];
    bool_t yup = TRUE, verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    uint16_image_t *pimg = float_image_to_uint16_image(img, isMask, chns, NULL, NULL, NULL, 255, yup, verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pimg, forceplain, verbose);
    uint16_image_free(pimg);
    fclose(wr);
    free(fname);
  }
