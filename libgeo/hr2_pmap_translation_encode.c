/* See {hr2_pmap_translation_encode.h}. */
/* Last edited on 2024-11-09 14:01:06 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <math.h>
 
#include <bool.h>
#include <r2.h>
#include <affirm.h>
#include <hr2_pmap_translation.h>

#include <hr2_pmap_translation_encode.h>

void hr2_pmap_translation_encode(hr2_pmap_t *M, double y[])
  {
    double A00 = M->dir.c[0][0];
    y[0] = M->dir.c[0][1] / A00;
    y[1] = M->dir.c[0][2] / A00;
    return;
  }

void hr2_pmap_translation_decode(double y[], hr2_pmap_t *M)
  { 
    r2_t disp = (r2_t){{ y[0], y[1] }};
    (*M) = hr2_pmap_translation_from_disp(&disp);

    return;
  }
