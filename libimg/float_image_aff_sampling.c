/* See {float_image_aff_sampling.h}. */
/* Last edited on 2020-11-05 23:19:34 by jstolfi */

#define _GNU_SOURCE
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <math.h>
 
#include <bool.h>
#include <sign.h>
#include <r2.h>
#include <r2x2.h>
#include <r2_aff_map.h>
#include <i2.h>
#include <jsmath.h>
#include <affirm.h>
#include <ix.h>
#include <wt_table.h>
#include <gauss_table.h>
#include <float_image.h>
#include <float_image_interpolate.h>

#include <float_image_aff_sampling.h>

r2_t float_image_aff_sampling_choose_step(r2x2_t *mat)
  {
    /* Assume that the matrix is roughly orthogonal. */
    /* !!! Find a better formula for highly non-orthogonal matrices !!! */
    
    /* Compute the sample distance for increment 1 in {x} and {y}: */
    double dx = hypot(mat->c[0][0], mat->c[0][1]);
    double dy = hypot(mat->c[1][0], mat->c[1][1]);

    /* Choose the sampling steps to be half a pixel after mapping: */
    return (r2_t){{ 0.5/dx, 0.5/dy }};
  }

i2_t float_image_aff_sampling_grid_size(r2_t dp, double R)
  {
    int32_t nx = (int32_t)ceil(R/dp.c[0]);
    int32_t ny = (int32_t)ceil(R/dp.c[1]);
    i2_t size = (i2_t){{ 2*nx + 1, 2*ny+1 }};
    return size;
  }
