/* See {multifok_image_read.h}. */
/* Last edited on 2025-04-11 09:01:17 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsprintf.h>

#include <float_image.h>
#include <float_image_read_gen.h>

#include <multifok_image.h>

#include <multifok_image_read.h>

/* FRAME-SPECIFC PNG FILES */
      
float_image_t *multifok_image_read_scene_view(char *frameFolder)
  {
    return multifok_image_read_png(frameFolder, "sVal", 0.0, 1.0);
  } 
  
float_image_t *multifok_image_read_height_average(char *frameFolder, double hMin, double hMax)
  {
    return multifok_image_read_png(frameFolder, "hAvg", (float)hMin, (float)hMax);
  }

float_image_t *multifok_image_read_height_deviation(char *frameFolder, double dMax)
  {
    return multifok_image_read_png(frameFolder, "hDev", 0.0, (float)dMax);
  }

float_image_t *multifok_image_read_normal_average(char *frameFolder)
  {
    return multifok_image_read_png(frameFolder, "sNrm", -1.0, +1.0);
  }

float_image_t *multifok_image_read_sharpness(char *frameFolder)
  {
    return multifok_image_read_png(frameFolder, "shrp", 0.0, 1.0);
  }

/* FRAME-SPECIFIC FNI FILES */

float_image_t *multifok_image_read_fni_height_average(char *frameFolder)
  {
    float_image_t *hAvg = multifok_image_read_fni(frameFolder, "hAvg", 1);
    return hAvg;
  }

float_image_t *multifok_image_read_fni_normal_average(char *frameFolder)
  {
    float_image_t *sNrm = multifok_image_read_fni(frameFolder, "sNrm", 3);
    return sNrm;
  }

/* GENERIC IMAGE FILE READING */

float_image_t *multifok_image_read_png(char *dir, char *name, float vMin, float vMax)
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

float_image_t *multifok_image_read_fni(char *dir, char *name, int32_t chns)
  { char *fname = jsprintf("%s/%s.fni", dir, name);
    float_image_t *img = float_image_read_named(fname);
    int32_t NC = (int32_t)(img->sz[0]);
    if ((chns >= 0) && (NC != chns) && (NC != chns+1))
      { fprintf(stderr, "read map %s with %d channels, should be %d or %d\n", fname, NC, chns, chns+1);  
        demand(FALSE, "aborted");
      }
    free(fname);
    return img;
  }

#define multifok_image_read_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

