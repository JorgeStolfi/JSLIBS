/* See {float_image_aff_sampling.h}. */
/* Last edited on 2023-11-25 18:19:23 by stolfi */

#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <math.h>
 
#include <bool.h>
#include <sign.h>
#include <r2.h>
#include <r3x3.h>
#include <i2.h>
#include <jsmath.h>
#include <affirm.h>
#include <ix.h>
#include <gauss_table.h>
#include <float_image.h>
#include <float_image_interpolate.h>

#include <float_image_aff_sampling.h>

r2_t float_image_aff_sampling_choose_step(r3x3_t *M)
  {
    /* Assume that the matrix is roughly orthogonal. */
    /* !!! Find a better formula for highly non-orthogonal matrices !!! */
    
    /* Compute the sample distance for increment 1 in {x} and {y}: */
    double w = fabs(M->c[0][0]);
    demand(w > 1.0e-200, "matrix weight is almost zero");
    demand(M->c[1][0] == 0, "map is not affine (1)");
    demand(M->c[2][0] == 0, "map is not affine (2)");
    double dx = hypot(M->c[1][1], M->c[1][2])/w;
    double dy = hypot(M->c[2][1], M->c[2][2])/w;

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
