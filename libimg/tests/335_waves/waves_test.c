#define PROG_NAME "waves_test"
#define PROG_DESC "test of {float_image_test.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-01-10 20:16:33 by stolfi */ 
/* Created on 2023-01-10 by J. Stolfi, UNICAMP */

#define waves_test_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <jsfile.h>
#include <bool.h>
#include <r2.h>
#include <float_image_write_pnm.h>
#include <float_image_waves.h>
#include <float_image_test.h>
#include <float_image.h>

int main(int argn, char **argv);

void do_test_waves(int32_t NF, double amp[], double fx[], double fy[], double phase[], int32_t kf0, int32_t kf1);
  /* Generates an image using the waves with indices {kf0..kf1} which must be in {0...NF-1}. */

void write_image(char *fname, float_image_t *img);

int main (int argn, char **argv)
  {
    fprintf(stderr, "choosing frequencies...\n");
    int32_t NF = 5;
    double amp[NF];
    double fx[NF];
    double fy[NF]; 
    double phase[NF];
    bool_t verbose = TRUE;
    float_image_waves_pick(NF, amp, fx, fy, phase, verbose);

    for(int32_t kf = 0; kf < NF; kf++)
      { do_test_waves(NF, amp, fx, fy, phase, kf, kf); }
    do_test_waves(NF, amp, fx, fy, phase, 0, 3);
    do_test_waves(NF, amp, fx, fy, phase, 0, NF-1);
   
    return 0;
  }

void do_test_waves(int32_t NF, double amp[], double fx[], double fy[], double phase[], int32_t kf0, int32_t kf1)
  {
    fprintf(stderr, "=== test float_image_waves = %d..%d ============\n", kf0, kf1);
    demand((0 <= kf0) && (kf0 <= kf1) && (kf1 < NF), "invalid wave indices");
    
    auto void wave_gen(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[]);
    
    fprintf(stderr, "creating image...\n");
    int NC = 1;
    int NX = 256;
    int NY = 192;
    float_image_t *img = float_image_new(NC, NX, NY);
    
    /* Choosing a squash amplitude: */
    double sum_a2 = 0;
    for (int32_t kf = kf0; kf <= kf1; kf++) { sum_a2 += amp[kf]*amp[kf]; }
    double amp_rms = sqrt(sum_a2); /* Root-mean-sum amplitude. */
    double squash = 0.5*amp_rms;

    fprintf(stderr, "filling image ...\n");
    float_image_test_paint(img, wave_gen, 3);
    
    /* Write the input and output images: */
    char *fname = NULL;
    asprintf(&fname, "out/test-kf%03d-%03d.pgm", kf0, kf1);
    write_image(fname, img);
    free(fname);
    fprintf(stderr, "===================================================\n");
    
    return;
    
    void wave_gen(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[])
      {
        /* Eval the waves: */
        int32_t KF = kf1 - kf0 + 1;
        double r = float_image_waves_eval(p->c[0], p->c[1], KF, &(amp[kf0]), &(fx[kf0]), &(fy[kf0]), &(phase[kf0]));
        
        /* Map {r} to {[-1 _ +1]}: */
        r = r/squash;
        r = r/hypot(1, r);
        
        /* Map {r} to {[0_1]}: */
        r = (r + 1.0)/2;
        assert(NC == 1);
        fs[0] = (float)r;
      }
  }  

void write_image(char *fname, float_image_t *img)
  {
    assert(img->sz[0] == 1);

    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    double gamma = 1.000;
    double bias = 0.000;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = TRUE;
    float_image_write_pnm_named(fname, img, isMask, gamma, bias, yup, warn, verbose);
  }
