/* See {neuromat_image_png.h}. */
/* Last edited on 2021-08-29 01:05:06 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <float_image.h>
#include <float_image_to_uint16_image.h>
#include <affirm.h>

#include <jsfile.h>
#include <uint16_image.h>
#include <uint16_image_write_png.h>
#include <jspnm.h>
#include <sample_conv.h>

#include <neuromat_image_png.h>

void neuromat_image_png_write(char *dir, char *name, float_image_t *fim, float vlo, float vhi, double gamma)
  {
    int32_t NC = (int32_t)fim->sz[0];
    int32_t chns = NC;
    demand((chns >= 1) && (chns <= 4), "channel count must be in {1..4}");
    
    /* Copy the image to a work image: */
    float_image_t *wim = float_image_copy(fim);
    
    /* Convert the float image to integer image {pim}: */
    if (isnan(gamma)) { gamma = sample_conv_BT709_ENC_GAMMA; }
    double bias = sample_conv_BT709_BIAS; 
    uint16_t maxval = 65535;
    bool_t yup = TRUE;
    bool_t verbose_cvt = FALSE;
    int32_t cop = ((NC == 1) || (NC == 3) ? -1 : NC-1); /* Index of opacity channel. */
    /* Map the image to {[0_1]} and apply gamma correction: */
    for (uint32_t c = 0;  c < NC; c++)  
       { if (c != cop) 
           { float_image_rescale_samples(wim, c, vlo, vhi, 0.0, 1.0);
             float_image_apply_gamma(wim, c, gamma, bias);
           }
       }
    uint16_image_t *pim = float_image_to_uint16_image(wim, FALSE, chns, NULL, NULL, NULL, maxval, yup, verbose_cvt);
    
    /* Write {pim} to disk: */
    char *fname = jsprintf("%s/%s.png", dir, name);
    bool_t verbose_write = TRUE;
    uint16_image_write_png_named(fname, pim, gamma, verbose_write);
    
    /* Cleanup: */
    float_image_free(wim);
    uint16_image_free(pim);
    free(fname);
  }
