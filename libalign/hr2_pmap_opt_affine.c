/* See {hr2_pmap_opt_affine.h}. */
/* Last edited on 2023-10-21 14:45:24 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>
#include <math.h>
 
#include <bool.h>
#include <r2.h>
#include <affirm.h>

#include <hr2_pmap_opt_affine.h>

void hr2_pmap_opt_affine_encode(hr2_pmap_t *M, double rad[], double y[])
  { y[0] = M->dir.c[0][1] / rad[0];
    y[1] = M->dir.c[0][2] / rad[1];
    y[2] = (M->dir.c[1][1] - 1.0) / rad[2];
    y[3] = M->dir.c[1][2] / rad[3];
    y[4] = M->dir.c[2][1] / rad[4];
    y[5] = (M->dir.c[2][2] - 1.0) / rad[5];
  }

void hr2_pmap_opt_affine_decode(double y[], double rad[], hr2_pmap_t *M)
  { r2_t v = (r2_t){{y[0] * rad[0], y[1] * rad[1] }};
    
    r2x2_t A;
    A.c[0][0] = y[2] * rad[2] + 1.0;
    A.c[0][1] = y[3] * rad[3];
    A.c[1][0] = y[4] * rad[4];
    A.c[1][1] = y[5] * rad[5] + 1.0;
    
    (*M) = hr2_pmap_aff_from_mat_and_disp(&A, &v);
  }
