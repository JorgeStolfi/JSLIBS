/* See {multifok_image_write.h}. */
/* Last edited on 2025-04-11 09:00:40 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <jsprintf.h>

#include <float_image.h>
#include <float_image_write_gen.h>
#include <float_image_brightness_scale.h>

#include <multifok_image.h>

#include <multifok_image_write.h>

/* FRAME-INDEPENDENT PNG FILES */

void multifok_image_write_sample_weights(float_image_t *wsimg, char *outFolder)
  { float vMin = 0.0f; 
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(wsimg, 0, &vMin, &vMax);
    demand(vMin == 0.0f, "image has negative sample weights");

    multifok_image_write_png(wsimg, outFolder, "weights", vMin, vMax);
  }
       
void multifok_image_write_basis_kernels(uint32_t NB, float_image_t *bKer[], char *outFolder)
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
        multifok_image_write_png(bKer[kb], outFolder, imgName, vMin, vMax);
        free(imgName);
      }
  }

void multifok_image_write_pixel_mask(float_image_t *pSel, float_image_t *bgrd, char *outFolder)
  {
    if (bgrd == NULL)
      { multifok_image_write_png(pSel, outFolder, "pSel", 0.0, 1.0); }
    else
      { int32_t NC, NX, NY;
        float_image_get_size(bgrd, &NC, &NX, &NY);
        float_image_check_size(pSel, 1, NX, NY, "bad pixel pick image");
        float_image_t *comp = float_image_new(NC, NX, NY);
        for (int32_t c = 0;  c < NC; c++)
          { float_image_mix_channels(0.75, pSel, 0, 0.25, bgrd, c, comp, c); }
        multifok_image_write_png(comp, outFolder, "pSel", 0.0, 1.0);
      }
  }

void multifok_image_write_selected_pixels(uint32_t NQ, i2_t pix[], float_image_t *sVal, char *outFolder)
  { int32_t NC = (int32_t)sVal->sz[0];
    demand((NC == 3) || (NC == 1), "background must be a color image");
    float_image_t *fimg = float_image_copy(sVal);
    multifok_image_draw_crosses(fimg, 0, NQ, pix, 1.0f);
    multifok_image_write_png(fimg, outFolder, "selected-pixels", 0.0f, 1.0f);
    float_image_free(fimg);
  }

/* FRAME-SPECIFC PNG FILES */

void multifok_image_write_scene_view(float_image_t *sVal, char *frameFolder)
  {
    multifok_image_write_png(sVal, frameFolder, "sVal", 0.0, 1.0);
  }

void multifok_image_write_height_average(float_image_t *hAvg, char *frameFolder, double hMin, double hMax)
  {
    float vMin = +INF;
    float vMax = -INF;
    float_image_update_sample_range(hAvg, 0, &vMin, &vMax);
    fprintf(stderr, "nominal Z range = [ %.6f _ %.6f ]", hMin, hMax);
    fprintf(stderr, " actual in hAvg = [ %.6f _ %.6f ]\n", vMin, vMax);

    multifok_image_write_png(hAvg, frameFolder, "hAvg", (float)hMin, (float)hMax);
  }

void multifok_image_write_height_deviation(float_image_t *hDev, char *frameFolder, double dMax)
  {
    multifok_image_write_png(hDev, frameFolder, "hDev", 0.0, (float)dMax);
  }

void multifok_image_write_normal_average(float_image_t *sNrm, char *frameFolder)
  {
    multifok_image_write_png(sNrm, frameFolder, "sNrm", -1.0, +1.0);
  }

void multifok_image_write_sharpness(float_image_t *shrp, char *frameFolder)
  {
    multifok_image_write_png(shrp, frameFolder, "shrp", 0.0, 1.0);
  }

void multifok_image_write_window_average(float_image_t *sAvg, char *frameFolder)
  { 
    multifok_image_write_png(sAvg, frameFolder, "sAvg", 0.0f, 1.0f);
  }

void multifok_image_write_window_gradient(float_image_t *sGrd, char *frameFolder)
  { int32_t NC = (int32_t)sGrd->sz[0];
    demand(NC == 3, "gradient image must have {NC=3}"); 
    float vMin = +INF; 
    float vMax = -INF;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(sGrd, 0, &vMin, &vMax);
    float_image_update_sample_range(sGrd, 1, &vMin, &vMax);
    vMax = (float)fmax(1.0e-38, fmax(fabs(vMin), fabs(vMax)));
    
    multifok_image_write_png(sGrd, frameFolder, "sGrd", -vMax, +vMax);
  };

void multifok_image_write_window_deviation(float_image_t *sDev, char *frameFolder)
  { float vMin = 0.0f; 
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(sDev, 0, &vMin, &vMax);
    demand(vMin == 0.0f, "negative samples in deviation image");

    multifok_image_write_png(sDev, frameFolder, "sDev", vMin, vMax);
  };

void multifok_image_write_focus_score(float_image_t *fVal, char *frameFolder)
  {
    float vMin = 0.0f; 
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(fVal, 0, &vMin, &vMax);
    demand(vMin == 0.0f, "negative samples in score image");

    multifok_image_write_png(fVal, frameFolder, "fVal", 0.0, 1.0);
  }

void multifok_image_write_focus_score_error(float_image_t *fErr, char *frameFolder)
  {
    multifok_image_write_png(fErr, frameFolder, "fErr", -1.0, +1.0);
  }
    
void multifok_image_write_estimated_height(float_image_t *zEst, char *frameFolder, double hMin, double hMax)
  {
    multifok_image_write_png(zEst, frameFolder, "zEst", (float)hMin, (float)hMax);
  }
 
void multifok_image_write_height_error(float_image_t *zErr, char *frameFolder, double hMin, double hMax)
  {
    double hMag = hMax - hMin;
    assert(hMag >= 1.0e-8);
    multifok_image_write_png(zErr, frameFolder, "zErr", -(float)hMag, +(float)hMag);
  }
      
void multifok_image_write_merged_scene_view(float_image_t *sMrg, char *frameFolder)
  {
    multifok_image_write_png(sMrg, frameFolder, "sMrg", 0.0, 1.0);
  }

void multifok_image_write_merged_scene_view_error(float_image_t *sErr, char *frameFolder)
  { 
    multifok_image_write_png(sErr, frameFolder, "sErr", -1.0f, +1.0f);
  }

/* FRAME-RELATED FNI FILES WITH WEIGHTS */

void multifok_image_write_fni_height_average(float_image_t *hAvgWht, char *frameFolder)
  {
    int32_t NC = (int32_t)(hAvgWht->sz[0]);
    demand(NC == 2, "wrong channel count in height map (no weights?)");
    multifok_image_write_fni(hAvgWht, frameFolder, "hAvg");
  }

void multifok_image_write_fni_normal_average(float_image_t *sNrmWht, char *frameFolder)
  {
    int32_t NC = (int32_t)(sNrmWht->sz[0]);
    demand(NC == 4, "wrong channel count in normal map (no weights?)");
    multifok_image_write_fni(sNrmWht, frameFolder, "sNrm");
  }

/* OPERATOR-RELATED PNG FILES */

void multifok_image_write_normalized(float_image_t *sNrm, char *frameFolder)
  { float vMin = -1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(sNrm, 0, &vMin, &vMax);
    double vR = fmax(fabs(vMin), fabs(vMax));
    vMin = (float)-vR;
    vMax = (float)+vR;

    multifok_image_write_png(sNrm, frameFolder, "sNrm", vMin, vMax);
  };

void multifok_image_write_basis_coeffs(uint32_t NB, float_image_t *bVal[], char *frameFolder, double bMax)
  {
    float vMax = (float)bMax;
    for (uint32_t kb = 0;  kb < NB; kb++)
      { char *imgName = jsprintf("bVal-%03d", kb);
        multifok_image_write_png(bVal[kb], frameFolder, imgName, -vMax, +vMax);
        free(imgName);
      }
  }

void multifok_image_write_basis_coeffs_squared(uint32_t NB, float_image_t *bSqr[], char *frameFolder, double bMax)
  {
    float vMin = 0.0; 
    float vMax = (float)(bMax*bMax);
    
    for (uint32_t kb = 0;  kb < NB; kb++)
      { char *imgName = jsprintf("bSqr-%03d", kb);
        multifok_image_write_png(bSqr[kb], frameFolder, imgName, vMin, vMax);
        free(imgName);
      }
  }

void multifok_image_write_quadratic_terms(uint32_t NT, float_image_t *tVal[], char *frameFolder, double tMax)
  {
    float vMin = 0.0; /* Assumes we only work with positive quadratic terms. */
    float vMax = (float)tMax;
    for (uint32_t kt = 0;  kt < NT; kt++)
      { char *imgName = jsprintf("tVal-%03d", kt);
        multifok_image_write_png(tVal[kt], frameFolder, imgName, vMin, vMax);
        free(imgName);
      }
  }

/* GENERIC IMAGE FILE WRITING */

void multifok_image_write_png(float_image_t *img, char *dir, char *name, float vMin, float vMax)
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

void multifok_image_write_fni(float_image_t *img, char *dir, char *name)
  { char *fname = jsprintf("%s/%s.fni", dir, name);
    float_image_write_named(fname, img);
    free(fname);
  }

#define multifok_image_write_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

