/* See {neuromat_image_png.h}. */
/* Last edited on 2013-12-06 04:51:59 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <float_image.h>
#include <float_pnm_image.h>
#include <affirm.h>

#include <jsfile.h>
#include <jspnm_image.h>
#include <jspng_image.h>
#include <jspnm.h>

#include <neuromat_image_png.h>

void neuromat_image_png_write(char *prefix, char *tag, float_image_t *fim, double vlo, double vhi)
  {
    /* Convert the float image to integer image {pim}: */
    int NC = (int)fim->sz[0];
    int chns = NC;
    assert((chns == 1) || (chns == 3));
    pnm_sample_t maxval = 65535;
    bool_t yup = TRUE;
    bool_t verbose_cvt = FALSE;
    pnm_image_t *pim = float_image_to_pnm_image(fim, FALSE, chns, NULL, NULL, NULL, maxval, yup, verbose_cvt);
    
    /* Write {pim} to disk: */
    char *fname = NULL;
    asprintf(&fname, "%s_%s.png", prefix, tag);
    double gamma = 1.0;
    bool_t verbose_write = TRUE;
    jspng_image_write(fname, pim, gamma, verbose_write);
    
    /* Cleanup: */
    pnm_image_free(pim);
    free(fname);
  }
