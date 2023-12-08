/* See {float_image_aff_compare.h}. */
/* Last edited on 2023-11-25 18:19:44 by stolfi */

#define _GNU_SOURCE
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
#include <float_image_aff_sampling.h>

#include <float_image_aff_compare.h>

/* INTERNAL PROTOTYPES */

void float_image_aff_compare_show_sample(char *label, r2_t *p, int32_t NC, double v[]);
  /* Prints the point {p} and the interpolated sample values {v[0..NC-1]}, prefixed by {label}. */

/* IMPLEMENTATIONS */

double float_image_aff_compare
  ( float_image_t *img1,
    hr2_pmap_t *A1,
    float_image_t *img2,
    hr2_pmap_t *A2,
    r2_t *stepP,
    i2_t *sizeP
  )
  {
    int32_t NC = (int32_t)img1->sz[0]; /* Number of color channels. */
    demand(NC == img2->sz[0], "images must have the same channel count");
    
    /* The maps must be affine: */
    demand(hr2_pmap_is_affine(A1), "map {A1} is not affine");
    demand(hr2_pmap_is_affine(A2), "map {A2} is not affine");
    
    bool_t debug_table = FALSE;
    bool_t debug_sampling = FALSE;
    
    /* Choose the steps in {x} and {y}: */
    r2_t step1 = float_image_aff_sampling_choose_step(&(A1->dir));
    r2_t step2 = float_image_aff_sampling_choose_step(&(A2->dir));
    double dx = fmin(step1.c[0], step2.c[0]);
    double dy = fmin(step1.c[1], step2.c[1]);
    if (debug_sampling) { fprintf(stderr, "step = (%.6f %.6f)\n", dx, dy); }
    r2_t step = (r2_t){{ dx, dy }};
    
    /* Choose the number of samples in each axis: */
    double R = 3.5; /* Enough to cover all significant weights. */
    i2_t size = float_image_aff_sampling_grid_size(step, R);
    int32_t hx = (size.c[0] - 1)/2; assert(size.c[0] == 2*hx +1);
    int32_t hy = (size.c[1] - 1)/2; assert(size.c[1] == 2*hy +1);
    if (debug_sampling) { fprintf(stderr, "nsmp = (%d %d)\n", hx, hy); }
    
    /* Create the weight tables: */
    bool_t normSum = FALSE;
    bool_t folded = FALSE;
    double *wx = gauss_table_make(hx+1, 0.0, 1.0/dx, normSum, folded);
    double *wy = gauss_table_make(hy+1, 0.0, 1.0/dy, normSum, folded);
    if (debug_table)
      { int32_t k = 0;
        while((k <= hx) || (k <= hy))
          { fprintf(stderr, "  %8.6f", (k <= hx ? wx[k] : 0));
            fprintf(stderr, "  %8.6f", (k <= hy ? wy[k] : 0));
            fprintf(stderr, "\n");
            k++;
          }
      }
    
    /* Compute the discrete integral: */
    double sum_wd2 = 0.0;
    double sum_w = 0.0;
    double v1[NC], v2[NC]; /* Interpolated pixel values. */
    int32_t order = 1; /* C1 interpolation. */
    ix_reduction_t red = ix_reduction_EXTEND; /* Replicate border pixels. */
    for (int32_t kx = -hx;  kx <= hx; kx++)
      { double wxk = (kx >= 0 ? wx[kx] : wx[-kx]);
        for (int32_t ky = -hy;  ky <= hy; ky++)
          { double wyk = (ky >= 0 ? wy[ky] : wy[-ky]);
            r2_t p = (r2_t){{ kx*dx, ky*dy }};
            
            r2_t p1 = hr2_pmap_r2_point(&p, A1); 
            float_image_interpolate_pixel(img1, p1.c[0], p1.c[1], order, red, v1);
            
            r2_t p2 = hr2_pmap_r2_point(&p, A2);
            float_image_interpolate_pixel(img2, p2.c[0], p2.c[1], order, red, v2);
            
            double w = wxk*wyk;
            double d2 = 0.0;
            for (int32_t kc = 0; kc < NC; kc++) { double d = v1[kc]-v2[kc]; d2 += d*d; }
            
            if (debug_sampling)
              { float_image_aff_compare_show_sample("C1(img1)", &p1, NC, v1);
                float_image_aff_compare_show_sample("C2(img2)", &p2, NC, v2);
                fprintf(stderr, "   %4d %4d", kx, ky);
                fprintf(stderr, " d2 = %14.12f w = %10.8f\n", d2, w);
              }
            
            sum_wd2 += d2*w;
            sum_w += w;
        }
      }
    assert(sum_w > 0.0);
    if (debug_sampling) { fprintf(stderr, "sum_w = %.4f\n", sum_w); }

    if (stepP != NULL) { (*stepP) = step; }
    if (sizeP != NULL) { (*sizeP) = size; }

    return sum_wd2/sum_w;
  }

void float_image_aff_compare_show_sample(char *label, r2_t *p, int32_t NC, double v[])
  { fprintf(stderr, "    %s", label);
    fprintf(stderr, " ( %10.8f %10.8f )", p->c[0], p->c[1]);
    fprintf(stderr, " =");
    for (int32_t kc = 0; kc < NC; kc++) { fprintf(stderr, " %+8.5f", v[kc]); }
    fprintf(stderr, "\n");
  }
