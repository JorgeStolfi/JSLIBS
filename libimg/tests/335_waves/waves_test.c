#define PROG_NAME "waves_test"
#define PROG_DESC "test of {float_image_test.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-01-22 03:50:36 by stolfi */ 
/* Created on 2023-01-10 by J. Stolfi, UNICAMP */

#define waves_test_COPYRIGHT \
  "Copyright � 2023  by the State University of Campinas (UNICAMP)"
  
#define PROG_INFO \
  "For various combinations of {KF0} and {KF1}, creates" \
  " images \"{outPrefix}-{KF0}{KF1}.pgm\" and  \"{outPrefix}-{KF0}-{KF1}.ppm\" showing" \
  " the waves patterns generated by {float_image_waves} including only the terms with" \
  " frequencies {KF0..KF1}.  The \".ppm\" images are colorized versions of" \
  " the \".pgm\" ones, with asimple linear interpolation color scale."

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
#include <frgb.h>
#include <float_image_write_pnm.h>
#include <float_image_waves.h>
#include <float_image_test.h>
#include <float_image.h>

int main(int argn, char **argv);

void wt_create_images_from_wave_range
  ( int32_t NF, 
    double amp[],
    double fx[],
    double fy[],
    double phase[], 
    int32_t kf0, 
    int32_t kf1
  );
  /* Generates images using the waves with indices {kf0..kf1} which must be in {0...NF-1},
    of various channel counts and sizes. */

void wt_create_images_from_wave_range_and_size
  ( int32_t NF, 
    double amp[],
    double fx[],
    double fy[],
    double phase[], 
    int32_t kf0, 
    int32_t kf1, 
    double squash,
    int32_t NC,
    int32_t NX,
    int32_t NY
  );
  /* Generates images using the waves with indices {kf0..kf1} which must be in {0...NF-1},
    with {NC} channels (1 or 3) and size {NX} by {NY}. */

void wt_write_image(int32_t NX, int32_t NY, int32_t kf0, int32_t kf1, float_image_t *img);

int main (int argn, char **argv)
  {
    fprintf(stderr, "choosing wave parameters...\n");
    int32_t NF = 10;
    double amp[NF];
    double fx[NF];
    double fy[NF]; 
    double phase[NF];
    bool_t verbose = TRUE;
    float_image_waves_pick(NF, amp, fx, fy, phase, verbose);
    
    for(int32_t kf = 0; kf < NF; kf++)
      { wt_create_images_from_wave_range(NF, amp, fx, fy, phase, kf, kf); }
    wt_create_images_from_wave_range(NF, amp, fx, fy, phase, 0, 3);
    wt_create_images_from_wave_range(NF, amp, fx, fy, phase, 0, NF-1);
   
    return 0;
  }

void wt_create_images_from_wave_range
  ( int32_t NF, 
    double amp[],
    double fx[],
    double fy[],
    double phase[], 
    int32_t kf0, 
    int32_t kf1
  )
  {
    fprintf(stderr, "=== test float_image_waves = %d..%d ============\n", kf0, kf1);
    demand((0 <= kf0) && (kf0 <= kf1) && (kf1 < NF), "invalid wave indices");

    /* Choosing a squash amplitude: */
    double sum_a2 = 0;
    for (int32_t kf = kf0; kf <= kf1; kf++) { sum_a2 += amp[kf]*amp[kf]; }
    double amp_rms = sqrt(sum_a2); /* Root-mean-sum amplitude. */
    double squash = 0.5*amp_rms;

    wt_create_images_from_wave_range_and_size(NF, amp, fx, fy, phase, kf0, kf1, squash, 3, 400, 400);
    wt_create_images_from_wave_range_and_size(NF, amp, fx, fy, phase, kf0, kf1, squash, 1, 512, 512);
    wt_create_images_from_wave_range_and_size(NF, amp, fx, fy, phase, kf0, kf1, squash, 1, 256, 256);
    wt_create_images_from_wave_range_and_size(NF, amp, fx, fy, phase, kf0, kf1, squash, 1, 128, 128);

  }

void wt_create_images_from_wave_range_and_size
  ( int32_t NF, 
    double amp[],
    double fx[],
    double fy[],
    double phase[], 
    int32_t kf0, 
    int32_t kf1, 
    double squash,
    int32_t NC,
    int32_t NX,
    int32_t NY
  )
  {
    auto void wave_gen(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[]);
    
    fprintf(stderr, "creating image...\n");
    int NC_gray = 1;
    float_image_t *gimg = float_image_new(NC_gray, NX, NY); /* Grayscale image */

    fprintf(stderr, "filling grayscale image ...\n");
    float_image_test_paint(gimg, wave_gen, 3);
    
    if (NC == NC_gray)
      {
        wt_write_image(NX, NY, kf0, kf1, gimg);
      }
    else
      { int NC_color = 3;
        demand(NC == NC_color, "invalid {NC}");
        /* Use the grayscale values to interpolate between two contrasting colors: */
        float_image_t *cimg = float_image_new(NC_color, NX, NY); /* Color image. */
        fprintf(stderr, "Generating the color image ...\n");
        frgb_t bg = (frgb_t){{ 1.0f, 0.9f, 0.0f }};
        frgb_t fg = (frgb_t){{ 0.0f, 0.1f, 0.8f }};
        for (int32_t ix = 0; ix < NX; ix++)
          { for (int32_t iy = 0; iy < NY; iy++)
              { double r = float_image_get_sample(gimg, 0, ix, iy);
                for (int32_t ic = 0; ic < NC_color; ic++)
                  { float ci = (float)(r*fg.c[ic] + (1-r)*bg.c[ic]);
                    float_image_set_sample(cimg, ic, ix, iy, ci);
                  }
              }
          }
        /* Write the color image images: */
        wt_write_image(NX, NY, kf0, kf1, cimg);
      }
    
    fprintf(stderr, "===================================================\n");
    
    return;
    
    void wave_gen(r2_t *p, int32_t NC1, int32_t NX1, int32_t NY1, float fs[])
      {
        assert((NC1 == 1) && (NX1 == NX) && (NY1 == NY));
        
        /* Eval the waves: */
        int32_t KF = kf1 - kf0 + 1;
        double r = float_image_waves_eval(p->c[0], p->c[1], KF, &(amp[kf0]), &(fx[kf0]), &(fy[kf0]), &(phase[kf0]), squash);
        
        /* Map {r} to {[0_1]}: */
        r = (r + 1.0)/2;
        fs[0] = (float)r;
      }
  }  

void wt_write_image(int32_t NX, int32_t NY, int32_t kf0, int32_t kf1, float_image_t *img)
  {
    int32_t NC = (int32_t)img->sz[0];
    
    char *ext = (NC == 1 ? "pgm" : "ppm"); 
    char *fname = NULL;
    asprintf(&fname, "out/%s-%03dx%03d/wavys-%d%d.%s", ext, NX, NY, kf0, kf1, ext);

    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    double gamma = 1.000;
    double bias = 0.000;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = TRUE;
    float_image_write_pnm_named(fname, img, isMask, gamma, bias, yup, warn, verbose);
    free(fname);
  }
