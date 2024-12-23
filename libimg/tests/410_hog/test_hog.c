#define PROG_NAME "test_hog"
#define PROG_DESC "test of {float_image_hog.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-21 11:24:28 by stolfi */ 
/* Created on 2008-10-05 by J. Stolfi, UNICAMP */

#define test_hog_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <uint16_image.h>
#include <uint16_image_read_pnm.h>
#include <float_image.h>
#include <float_image_gradient.h>
#include <float_image_hog.h>
#include <float_image_from_uint16_image.h>
#include <sample_conv.h>
#include <sample_conv_gamma.h>
#include <bool.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <jsrandom.h>

#define BT_ENC_EXPO sample_conv_gamma_BT709_ENC_EXPO
#define BT_ENC_BIAS sample_conv_gamma_BT709_BIAS 
  /* Values of {expo} and {bias} parameters for {sample_conv_gamma}
    that approximate the BT.709 encoding.  */

typedef enum { kind_PEAK, kind_WAVE, kind_BUMP, kind_REAL } kind_t;
#define kind_LAST (kind_REAL)

int main(int argn, char **argv);

void do_test
  ( char *prefix,
    bool_t masked,
    bool_t oriented,
    double noise,
    int nh
  );
  /* Tests {float_image_hog_collect} with the gradient of image {I}
    and mask {M}, with parameters {orietned} and {noise}. */
  
float_image_t *read_image(char *prefix, char *suffix, char *ext);
  /* Reads an image file called "in/{prefix}-{suffix}.{ext}", where {ext} is "ppm" or "pgm" */
  
void write_image(char *prefix, float_image_t *A);
  /* Writes the image {A} to a file called "{prefix}-G.fni",
    in the FNI format, without gamma correction. */

void write_hog(char *prefix, int nh, double h[], bool_t masked, bool_t oriented, double noise);
  /* Writes the histogram {h[0..nh-1]} as file
    "out/{prefix}-{nh}-{masked}-{oriented}-{noise*1000}.txt". */ 

int main (int argn, char **argv)
  {
    srandom(4615*417);
    
    /* Parse arguments: */
    char *rest;

    int na = 1;                                      /* Current argument index. */
    
    char *prefix = argv[na++];                       /* Prefix for inut and output files. */

    double noise;
    noise = strtod(argv[na++],&rest);                /* Standard deviation of gradient noise. */
    demand((*rest) == 0, "invalid {noise}");

    int nh = (int)strtol(argv[na++],&rest,10);            /* Number of bins in histogram. */
    demand((*rest) == 0, "invalid {nh}");

    /* Run the test: */
    do_test(prefix, FALSE, FALSE, noise, nh);
    do_test(prefix, TRUE,  FALSE, noise, nh);
    do_test(prefix, TRUE,  TRUE,  noise, nh);

    return 0;
  }

void do_test
  ( char *prefix,
    bool_t masked,
    bool_t oriented,
    double noise,
    int nh
  )
  {
    fprintf(stderr, "------------------------------------------------------------\n");
    fprintf(stderr, "testing with image %s.ppm\n", prefix); 
    
    /* Get the input image: */
    float_image_t *I = read_image(prefix, "img", "ppm");
    int cI = 0; /* Channel of image {I} to process. */
    int NX = (int)I->sz[1]; 
    int NY = (int)I->sz[2];
    
    /* Get the mask, if any: */
    float_image_t *M = (masked ? read_image(prefix, "msk", "pgm") : NULL);
    
    /* Compute the gradient of {I} channel {cI}, save in {G} channels 0 (DX) and 1 (DY): */
    float_image_t *G = float_image_new(3,NX,NY);
    float_image_gradient_sobel(I, cI, G, 0, G, 1);

    /* Save the original image in {G}, channel 2: */
    float_image_set_channel(G, 2, I, cI);
    
    /* Write the original+gradient image: */
    write_image(prefix, G);
    
    /* Collect the HOG: */
    double h[nh];
    float_image_hog_collect(G, 0, G, 1, M, noise, oriented, nh, h);
    
    /* Write out the HOG: */
    write_hog(prefix, nh, h, masked, oriented, noise);

    float_image_free(I);
    if (M != NULL) { float_image_free(M); }
    float_image_free(G);

    fprintf(stderr, "------------------------------------------------------------\n");
  }  

float_image_t *read_image(char *prefix, char *suffix, char *ext)
  { char *fname = jsprintf("in/%s-%s.%s", prefix, suffix, ext);
    FILE *rd = open_read(fname, TRUE);
    uint16_image_t *pim = uint16_image_read_pnm_file(rd);
    fclose(rd);
    free(fname);
    bool_t yup = TRUE, verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, NULL, NULL, yup, verbose);
    uint16_image_free(pim);
    int c;
    for (c = 0; c < fim->sz[0]; c++) { float_image_apply_gamma(fim, c, 1/BT_ENC_EXPO, BT_ENC_BIAS); }
    return fim;
  }

void write_image(char *prefix, float_image_t *A)
  { char *fname = jsprintf("out/%s-G.fni", prefix);
    FILE *wr = open_write(fname, TRUE);
    float_image_write(wr, A);
    fclose(wr);
    free(fname);
  }

void write_hog(char *prefix, int nh, double h[], bool_t masked, bool_t oriented, double noise)
  {
    /* Write to PPM file: */
    char *fname = jsprintf("out/%s-%04d-%d-%d-%04d.txt", prefix, nh, masked, oriented, (int)floor(noise*1000+0.5));
    FILE *wr = open_write(fname, TRUE);
    double span = (oriented ? 2*M_PI : M_PI);
    double step = span/nh;
    int ih;
    for (ih =0; ih < nh; ih++)
      { 
        double amd = ih*step;
        double alo = amd - 0.5*step;
        double ahi = amd + 0.5*step;
        fprintf(wr, "%4d %8.5f %8.5f %8.5f  %12.9f\n", ih, amd, alo, ahi, h[ih]);
      }
    fclose(wr);
    
    /* Cleanup: */
    free(fname);
  }
