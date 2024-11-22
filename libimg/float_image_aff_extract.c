/* See {float_image_aff_extract.h}. */
/* Last edited on 2024-10-12 18:24:15 by stolfi */

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
#include <r2x2.h>
#include <hr2.h>
#include <i2.h>
#include <jsmath.h>
#include <affirm.h>
#include <ix.h>
#include <gauss_table.h>
#include <float_image.h>
#include <float_image_interpolate.h>

#include <float_image_aff_extract.h>

/* INTERNAL PROTOTYPES */

void float_image_aff_extract_show_sample(int32_t ix, int32_t iy, r2_t p, r2_t q, double w, int32_t NC, float vf[]);
  /* Prints to stderr the pixel indices {ix,iy}, the corresponding
    sampling point {p} of {\RR^2}, the mask weight {w}, the image {q} or
    {p} by the affine map, and the weighted interpolated image pixel
    value {vf[0..NC-1]}. */

/* IMPLEMENTATIONS */

float_image_t *float_image_aff_extract(float_image_t *img, hr2_pmap_t *A, r2_t dp, i2_t size)
  {
    bool_t debug_sampling = FALSE;
    
    /* The map must be affine: */
    demand(hr2_pmap_is_affine(A, 1.0e-12), "map {A} is not affine");
    
    /* Get the result image size: */
    int32_t NC = (int32_t)img->sz[0]; /* Number of color channels. */
    int32_t NX = size.c[0];
    int32_t NY = size.c[1];
    if (debug_sampling) { fprintf(stderr, "extracted feature size = (%d %d)\n", NX, NY); } 
    demand((NY % 2) == 1, "result row count must be odd");
    demand((NX % 2) == 1, "result column count must be odd");
    int32_t hx = (NX - 1)/2;
    int32_t hy = (NY - 1)/2;
    
    /* Get the steps in {x} and {y}: */
    double dx = dp.c[0];
    double dy = dp.c[1];
    if (debug_sampling) { fprintf(stderr, "step = (%.6f %.6f)\n", dx, dy); }
    demand(dx > 0.0, "step in X must be positive");
    demand(dy > 0.0, "step in Y must be positive");
    
    /* Create the weight tables: */
    bool_t normSum = FALSE;
    bool_t folded = FALSE;
    double *wx = gauss_table_make(hx+1, 0.0, 1.0/dx, normSum, folded);
    double *wy = gauss_table_make(hy+1, 0.0, 1.0/dy, normSum, folded);
    
    /* Extract the feature: */
    float_image_t *res = float_image_new(NC, NX, NY);
    double vd[NC]; /* Interpolated pixel value, as double. */
    float vf[NC]; /* Interpolated pixel value, as float. */
    int32_t order = 1; /* C1 interpolation. */
    ix_reduction_t red = ix_reduction_EXTEND; /* Replicate border pixels. */
    for (int32_t ix = -hx;  ix <= hx; ix++)
      { double wxi = (ix >= 0 ? wx[ix] : wx[-ix]);
        for (int32_t iy = -hy;  iy <= hy; iy++)
          { double wyi = (iy >= 0 ? wy[iy] : wy[-iy]);
            /* Get the image sample {vd[0..NC-1]}: */
            r2_t p = (r2_t){{ ix*dx, iy*dy }};  /* Grid point in {\RR^2}. */
            r2_t q = hr2_pmap_r2_point(&p, A); /* Affine map image of {p}. */
            float_image_interpolate_pixel(img, q.c[0], q.c[1], order, red, vd);
            /* Apply the mask weight: */
            double wxy = wxi*wyi;
            for (int32_t ic = 0; ic < NC; ic++) { vf[ic] = (float)(wxy*vd[ic]); }
            if (debug_sampling) { float_image_aff_extract_show_sample(ix, iy, p, q, wxy, NC, vf); }
            float_image_set_pixel(res, ix + hx, iy + hy, vf);
        }
      }
    return res;
  }

void float_image_aff_extract_show_sample(int32_t ix, int32_t iy, r2_t p, r2_t q, double w, int32_t NC, float vf[])
  { fprintf(stderr, "    ");
    fprintf(stderr, "   %4d %4d", ix, iy);
    fprintf(stderr, " ( %10.8f %10.8f )", p.c[0], p.c[1]);
    fprintf(stderr, " --> ( %10.8f %10.8f )", q.c[0], q.c[1]);
    fprintf(stderr, " * %10.8f\n", w);
    fprintf(stderr, " =");
    for (int32_t kc = 0; kc < NC; kc++) { fprintf(stderr, " %+8.5f", vf[kc]); }
    fprintf(stderr, "\n");
  }

