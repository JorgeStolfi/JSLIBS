/* See {float_image_align_local.h}. */
/* Last edited on 2022-10-30 19:45:20 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
 
#include <bool.h>
#include <r2.h>
#include <rmxn.h>
#include <jsmath.h>
#include <affirm.h>
#include <float_image.h>

#include <float_image_align_local.h>

void float_image_align_local
  ( int32_t ni,               /* Number of images to align. */
    float_image_t *img[],     /* The images. */
    r2_t p[],                 /* (IN/OUT) Corresponding points in each image. */
    double crad,              /* Radius of comparison neighborhood (pixels). */
    r2x2_t L[],               /* Local linear mapping for each image, or {NULL}. */
    r2_t arad[],              /* Max alignment adjustment for each image. */
    bool_t quadopt,           /* True uses quadratic optimization, false uses exhausive enum. */ 
    double tol,               /* Desired precision. */
    float_image_align_local_report_t *rep  /* Reporting function, or {NULL}.
  )
  {
  
    /* Weight table for comparison: */
    int32_t hw = ceil(2*crad);
    int32_t nw = 2*hw + 1;
    double wt[nw];
    bool_t norm = TRUE;
    wt_table_fill_hann(nw, wt, norm)

    auto double F2(int32_t ni1, r2_t p1[]);
      /* A quadratic mismatch function
        of the images in the neighborhood of {p}. */
    
    if (L != NULL) 
      { /* Create deformed copies of the images. */
    else
      { /* Use the given images. */ } 
  
  }
