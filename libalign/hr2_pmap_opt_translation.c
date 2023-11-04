/* See {hr2_pmap_opt_translation.h}. */
/* Last edited on 2023-10-20 15:15:36 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>
#include <math.h>
 
#include <bool.h>
#include <r2.h>
#include <affirm.h>

#include <hr2_pmap_opt_translation.h>

void hr2_pmap_opt_translation_encode(hr2_pmap_t *M, double rad[], double y[])
  {
    y[0] = M->dir.c[0][1] / rad[0];
    y[1] = M->dir.c[0][2] / rad[1];
    return;
  }

void hr2_pmap_opt_translation_decode(double y[], double rad[], hr2_pmap_t *M)
  { 
    M->dir.c[0][1] = y[0] * rad[0];
    M->dir.c[0][2] = y[1] * rad[1];

    M->inv.c[0][1] = - M->dir.c[0][1];
    M->inv.c[0][2] = - M->dir.c[0][2];
    return;
  }
