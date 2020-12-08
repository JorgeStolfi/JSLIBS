/* See {neuromat_image_png.h}. */
/* Last edited on 2018-12-04 22:01:28 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <float_image.h>
#include <float_image_to_uint16_image.h>
#include <affirm.h>

#include <jsfile.h>
#include <uint16_image.h>
#include <uint16_image_write_png.h>
#include <jspnm.h>

#include <neuromat_image_png.h>

void neuromat_image_png_write(char *prefix, char *tag, float_image_t *fim, double vlo, double vhi)
  {
    /* Convert the float image to integer image {pim}: */
    int NC = (int)fim->sz[0];
    int chns = NC;
    assert((chns == 1) || (chns == 3));
    uint16_t maxval = 65535;
    bool_t yup = TRUE;
    bool_t verbose_cvt = FALSE;
    uint16_image_t *pim = float_image_to_uint16_image(fim, FALSE, chns, NULL, NULL, NULL, maxval, yup, verbose_cvt);
    
    /* Write {pim} to disk: */
    char *fname = NULL;
    asprintf(&fname, "%s_%s.png", prefix, tag);
    double gamma = 1.0;
    bool_t verbose_write = TRUE;
    uint16_image_write_png_named(fname, pim, gamma, verbose_write);
    
    /* Cleanup: */
    uint16_image_free(pim);
    free(fname);
  }
