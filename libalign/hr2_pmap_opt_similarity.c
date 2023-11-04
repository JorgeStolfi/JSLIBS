/* See {hr2_pmap_opt_similarity.h}. */
/* Last edited on 2023-10-21 12:10:08 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>
#include <math.h>
 
#include <bool.h>
#include <r2.h>
#include <affirm.h>

#include <hr2_pmap_opt_similarity.h>

void hr2_pmap_opt_similarity_encode(hr2_pmap_t *M, double rad[], double y[])
  { 
    y[0] = M->dir.c[0][1] / rad[0];
    y[1] = M->dir.c[0][2] / rad[1];
    y[2] = (M->dir.c[1][1] - 1.0) / rad[2];
    y[3] = M->dir.c[1][2] / rad[3];
  }

void hr2_pmap_opt_similarity_decode(double y[], double rad[], hr2_pmap_t *M)
  {
    M->dir.c[0][1] = y[0] * rad[0];
    M->dir.c[0][2] = y[1] * rad[1];

    r2_t u = (r2_t){{ 1.0 + y[2] * rad[2], y[3] * rad[3] }};
    r2_t v = (r2_t){{ - u.c[1], + u.c[0] }};
    M->dir.c[1][1] = u.c[0];
    M->dir.c[1][2] = u.c[1];
    M->dir.c[2][1] = v.c[0];
    M->dir.c[2][2] = v.c[1];

    double m = r2_norm_sqr(&u);

    M->inv.c[0][1] = - (u.c[0]*M->dir.c[0][1] + u.c[1]*M->dir.c[0][2]) / m;
    M->inv.c[0][2] = - (v.c[0]*M->dir.c[0][1] + v.c[1]*M->dir.c[0][2]) / m;
    M->inv.c[1][1] = u.c[0] / m;
    M->inv.c[1][2] = v.c[0] / m;
    M->inv.c[2][1] = u.c[1] / m;
    M->inv.c[2][2] = v.c[1] / m;
  }
