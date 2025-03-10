/* See {multifok_image.h}. */
/* Last edited on 2025-02-02 02:23:18 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <wt_table.h>
#include <wt_table_hann.h>
#include <affirm.h>
#include <interval.h>
#include <fget.h>
#include <i2.h>
#include <r2.h>
#include <r3.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <jsqroots.h>
#include <rn.h>

#include <float_image.h>
#include <float_image_paint.h>
#include <float_image_read_gen.h>
#include <float_image_write_gen.h>
#include <float_image_brightness_scale.h>

#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_score.h>
#include <multifok_term.h>
#include <multifok_scene.h>

#include <multifok_image.h>

/* FRAME-INDEPENDENT IMAGES */

void multifok_image_sample_weights_write(float_image_t *wsimg, char *outFolder)
  { float vMin = 0.0f; 
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(wsimg, 0, &vMin, &vMax);
    demand(vMin == 0.0f, "image has negative sample weights");

    multifok_image_write(wsimg, outFolder, "weights", vMin, vMax);
  }
       
void multifok_image_basis_kernels_write(uint32_t NB, float_image_t *bKer[], char *outFolder)
  {
    float vMin = -1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    for (uint32_t kb = 0;  kb < NB; kb++)
      { float_image_update_sample_range(bKer[kb], 0, &vMin, &vMax); }

    double vR = fmax(fabs(vMin), fabs(vMax));
    vMin = (float)-vR;
    vMax = (float)+vR;

    for (uint32_t kb = 0;  kb < NB; kb++)
      { char *imgName = jsprintf("basis-%03d", kb);
        multifok_image_write(bKer[kb], outFolder, imgName, vMin, vMax);
        free(imgName);
      }
  }

void multifok_image_pixel_mask_write(float_image_t *pSel, float_image_t *bgrd, char *outFolder)
  {
    if (bgrd == NULL)
      { multifok_image_write(pSel, outFolder, "pSel", 0.0, 1.0); }
    else
      { int32_t NC, NX, NY;
        float_image_get_size(bgrd, &NC, &NX, &NY);
        float_image_check_size(pSel, 1, NX, NY, "bad pixel pick image");
        float_image_t *comp = float_image_new(NC, NX, NY);
        for (int32_t c = 0;  c < NC; c++)
          { float_image_mix_channels(0.75, pSel, 0, 0.25, bgrd, c, comp, c); }
        multifok_image_write(comp, outFolder, "pSel", 0.0, 1.0);
      }
  }

void multifok_image_selected_pixels_write(uint32_t NQ, i2_t pix[], float_image_t *sVal, char *outFolder)
  { int32_t NC = (int32_t)sVal->sz[0];
    demand((NC == 3) || (NC == 1), "background must be a color image");
    float_image_t *fimg = float_image_copy(sVal);
    multifok_image_draw_crosses(fimg, 0, NQ, pix, 1.0f);
    multifok_image_write(fimg, outFolder, "selected-pixels", 0.0f, 1.0f);
    float_image_free(fimg);
  }

/* FRAME-SPECIFC IMAGES */

void multifok_image_scene_view_write(float_image_t *sVal, char *frameFolder)
  {
    multifok_image_write(sVal, frameFolder, "sVal", 0.0, 1.0);
  }
      
float_image_t *multifok_image_scene_view_read(char *frameFolder)
  {
    return multifok_image_read(frameFolder, "sVal", 0.0, 1.0);
  } 


void multifok_image_height_average_write(float_image_t *hAvg, char *frameFolder, double hMin, double hMax)
  {
    float vMin = +INF;
    float vMax = -INF;
    float_image_update_sample_range(hAvg, 0, &vMin, &vMax);
    fprintf(stderr, "nominal Z range = [ %.6f _ %.6f ]", hMin, hMax);
    fprintf(stderr, " actual in hAvg = [ %.6f _ %.6f ]\n", vMin, vMax);

    multifok_image_write(hAvg, frameFolder, "hAvg", (float)hMin, (float)hMax);
  }
  
float_image_t *multifok_image_height_average_read(char *frameFolder, double hMin, double hMax)
  {
    return multifok_image_read(frameFolder, "hAvg", (float)hMin, (float)hMax);
  }


void multifok_image_height_deviation_write(float_image_t *hDev, char *frameFolder, double hMin, double hMax)
  {
    double hMag = (hMax-hMin)/2;
    multifok_image_write(hDev, frameFolder, "hDev", 0.0, (float)hMag);
  }

float_image_t *multifok_image_height_deviation_read(char *frameFolder, double hMin, double hMax)
  {
    double hMag = (hMax-hMin)/2;
    return multifok_image_read(frameFolder, "hDev", 0.0, (float)hMag);
  }


void multifok_image_normal_average_write(float_image_t *sNrm, char *frameFolder)
  {
    multifok_image_write(sNrm, frameFolder, "sNrm", -1.0, +1.0);
  }

float_image_t *multifok_image_normal_average_read(char *frameFolder)
  {
    return multifok_image_read(frameFolder, "sNrm", -1.0, +1.0);
  }


void multifok_image_sharpness_write(float_image_t *shrp, char *frameFolder)
  {
    multifok_image_write(shrp, frameFolder, "shrp", 0.0, 1.0);
  }

float_image_t *multifok_image_sharpness_read(char *frameFolder)
  {
    return multifok_image_read(frameFolder, "shrp", 0.0, 1.0);
  }

void multifok_image_window_average_write(float_image_t *sAvg, char *frameFolder)
  { 
    multifok_image_write(sAvg, frameFolder, "sAvg", 0.0f, 1.0f);
  }

void multifok_image_window_gradient_write(float_image_t *sGrd, char *frameFolder)
  { int32_t NC = (int32_t)sGrd->sz[0];
    demand(NC == 3, "gradient image must have {NC=3}"); 
    float vMin = +INF; 
    float vMax = -INF;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(sGrd, 0, &vMin, &vMax);
    float_image_update_sample_range(sGrd, 1, &vMin, &vMax);
    vMax = (float)fmax(1.0e-38, fmax(fabs(vMin), fabs(vMax)));
    
    multifok_image_write(sGrd, frameFolder, "sGrd", -vMax, +vMax);
  };

void multifok_image_window_deviation_write(float_image_t *sDev, char *frameFolder)
  { float vMin = 0.0f; 
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(sDev, 0, &vMin, &vMax);
    demand(vMin == 0.0f, "negative samples in deviation image");

    multifok_image_write(sDev, frameFolder, "sDev", vMin, vMax);
  };

void multifok_image_normalized_write(float_image_t *sNrm, char *frameFolder)
  { float vMin = -1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(sNrm, 0, &vMin, &vMax);
    double vR = fmax(fabs(vMin), fabs(vMax));
    vMin = (float)-vR;
    vMax = (float)+vR;

    multifok_image_write(sNrm, frameFolder, "sNrm", vMin, vMax);
  };


void multifok_image_focus_score_write(float_image_t *fVal, char *frameFolder)
  {
    float vMin = 0.0f; 
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(fVal, 0, &vMin, &vMax);
    demand(vMin == 0.0f, "negative samples in score image");

    multifok_image_write(fVal, frameFolder, "fVal", 0.0, 1.0);
  }

void multifok_image_focus_score_error_write(float_image_t *fErr, char *frameFolder)
  {
    multifok_image_write(fErr, frameFolder, "fErr", -1.0, +1.0);
  }


    
void multifok_image_estimated_height_write(float_image_t *zEst, char *frameFolder, double hMin, double hMax)
  {
    multifok_image_write(zEst, frameFolder, "zEst", (float)hMin, (float)hMax);
  }
 
void multifok_image_height_error_write(float_image_t *zErr, char *frameFolder, double hMin, double hMax)
  {
    double hMag = hMax - hMin;
    assert(hMag >= 1.0e-8);
    multifok_image_write(zErr, frameFolder, "zErr", -(float)hMag, +(float)hMag);
  }

      
void multifok_image_merged_scene_view_write(float_image_t *sMrg, char *frameFolder)
  {
    multifok_image_write(sMrg, frameFolder, "sMrg", 0.0, 1.0);
  }

void multifok_image_merged_scene_view_error_write(float_image_t *sErr, char *frameFolder)
  { 
    multifok_image_write(sErr, frameFolder, "sErr", -1.0f, +1.0f);
  }


void multifok_image_basis_coeffs_write(uint32_t NB, float_image_t *bVal[], char *frameFolder, double bMax)
  {
    float vMax = (float)bMax;
    for (uint32_t kb = 0;  kb < NB; kb++)
      { char *imgName = jsprintf("bVal-%03d", kb);
        multifok_image_write(bVal[kb], frameFolder, imgName, -vMax, +vMax);
        free(imgName);
      }
  }

void multifok_image_basis_coeffs_squared_write(uint32_t NB, float_image_t *bSqr[], char *frameFolder, double bMax)
  {
    float vMin = 0.0; 
    float vMax = (float)(bMax*bMax);
    
    for (uint32_t kb = 0;  kb < NB; kb++)
      { char *imgName = jsprintf("bSqr-%03d", kb);
        multifok_image_write(bSqr[kb], frameFolder, imgName, vMin, vMax);
        free(imgName);
      }
  }

void multifok_image_quadratic_terms_write(uint32_t NT, float_image_t *tVal[], char *frameFolder, double tMax)
  {
    float vMin = 0.0; /* Assumes we only work with positive quadratic terms. */
    float vMax = (float)tMax;
    for (uint32_t kt = 0;  kt < NT; kt++)
      { char *imgName = jsprintf("tVal-%03d", kt);
        multifok_image_write(tVal[kt], frameFolder, imgName, vMin, vMax);
        free(imgName);
      }
  }

/* GENERIC IMAGE I/O */

float_image_t *multifok_image_read(char *dir, char *name, float vMin, float vMax)
  {
    char *fname = jsprintf("%s/%s.png", dir, name);
    image_file_format_t ffmt = image_file_format_PNG;
    bool_t yUp = FALSE;
    double expo_dec, bias;
    uint16_t *maxval = NULL;
    bool_t verbose = FALSE;
    float_image_t *img = float_image_read_gen_named(fname, ffmt, yUp, vMin, vMax, &maxval, &expo_dec, &bias, verbose);
    
    free(fname);
    if (maxval != NULL) { free(maxval); }
    return img;
  }

void multifok_image_write(float_image_t *img, char *dir, char *name, float vMin, float vMax)
  { 
    /* Add the linear brightness scale: */
    bool_t yLo = FALSE; /* Add scale after last row. */
    float_image_t *simg = float_image_brightness_scale_add(img, yLo, vMin, vMax);
  
    char *fname = jsprintf("%s/%s.png", dir, name);
    image_file_format_t ffmt = image_file_format_PNG;
    bool_t yUp = TRUE;
    double gammaEnc = multifok_image_gamma;
    double bias = multifok_image_bias;
    bool_t verbose = FALSE;
    if (! verbose) { fprintf(stderr, "writing %s ...\n", fname); }
    float_image_write_gen_named(fname, simg, ffmt, yUp, vMin, vMax, gammaEnc, bias, verbose);
    free(fname);
    float_image_free(simg);
  }

void multifok_image_draw_crosses(float_image_t *img, int32_t ch, uint32_t NQ, i2_t pix[], float val)
  { 
    int32_t NC = (int32_t)img->sz[0];
    
    /* Draw crosses over the image: */
    double rad = 6.0;
    double hwd = 0.5;
    bool_t empty = TRUE;
    bool_t diagonal = TRUE;
    for (uint32_t kq = 0;  kq < NQ; kq++)
      { double xctr = pix[kq].c[0] + 0.5;
        double yctr = pix[kq].c[1] + 0.5;
        for (int32_t ic = 0;  ic < NC; ic++)
          { float_image_paint_cross(img, ic, xctr, yctr, rad, empty, 2.0*hwd, diagonal, 0.0f, 3);
            float valk = (ic == ch ? val : 0.0f);
            float_image_paint_cross(img, ic, xctr, yctr, rad, empty, hwd, diagonal, valk, 3);
          }
      }
  }

/* SPECIALIZED IMAGE OUTPUT */

#define multifok_image_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

