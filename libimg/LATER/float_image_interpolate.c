/* See float_image_interpolate.h */
/* Last edited on 2024-12-04 23:21:14 by stolfi */ 

#include <limits.h>
#include <float.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <ix_types.h>
#include <ix_reduce.h>
#include <float_image.h>
#include <spline_interp.h>

#include <float_image_interpolate.h>

/* INTERNAL PROTOTYPES */

void float_image_interpolate_get_samples_and_weights
  ( float_image_t *A, 
    int32_t c, 
    double x, 
    double y,
    int32_t m,
    ix_reduce_mode_t red, 
    float *p[],
    double wx[],
    double wy[]
  );
  /* Stores in {p[0..m*m-1]} the addresses of the {m×m} block samples
    of channel {c} that are involved in the interpolation of image {A}
    at the point {(x,y)}. Indices {i} outside the range {0..N-1} are
    reduced according to {ix_reduce(i, N, red)}; If a reduced index is
    {-1}, the corresponding pointers are set to NULL. */

/* EXPORTED FUNCTIONS */

double float_image_interpolate_sample
  ( float_image_t *A, 
    int32_t c, 
    double x, 
    double y, 
    int32_t order, 
    ix_reduce_mode_t red
  )
  { 
    int32_t m = float_image_interpolate_compute_num_samples(order);
    
    /* Get pixel values and weights: */
    float *p[m*m];
    double wx[m];
    double wy[m];
    float_image_interpolate_get_samples_and_weights(A, c, x, y, m, red, p, wx, wy);
    
    /* Apply interpolation formula: */
    double sum = 0;
    int32_t kx, ky;
    for (ky = 0; ky < m; ky++)
      { float **py = &(p[m*ky]);
        double sumy = 0;
        for (kx = 0; kx < m; kx++)
          { float *pxy = py[kx];
            if (pxy != NULL) { sumy += wx[kx]*(*pxy); }
          }
        sum += wy[ky]*sumy;
      }
    return sum;
  }
  
void float_image_interpolate_pixel
  ( float_image_t *A, 
    double x, 
    double y, 
    int32_t order, 
    ix_reduce_mode_t red, 
    double z[]
  )
  { int32_t m = float_image_interpolate_compute_num_samples(order);
    
    /* Get pixel values and weights: */
    float *p[m*m];
    double wx[m];
    double wy[m];
    float_image_interpolate_get_samples_and_weights(A, 0, x, y, m, red, p, wx, wy);
    
    /* Apply interpolation formula to all channels: */
    int32_t NC = A->sz[0];         /* Number of channels. */
    ix_step_t cst = A->st[0];  /* Position increment between channels. */
    int32_t c;
    for (c = 0; c < NC; c++)
      { double sum = 0;
        int32_t kx, ky;
        for (ky = 0; ky < m; ky++)
          { float **py = &(p[m*ky]);
            double sumy = 0;
            for (kx = 0; kx < m; kx++)
              { float *pxy = py[kx];
                if (pxy != NULL) { sumy += wx[kx]*(*pxy); py[kx] += cst; }
              }
            sum += wy[ky]*sumy;
          }
        z[c] = sum;
      }
  }

/* !!! TO DO: Improve efficiency of grid evaluation by reusing samples and weights !!! */

void float_image_interpolate_grid_samples
  ( float_image_t *A, 
    int32_t c, 
    double rx, int32_t hx, double dx,
    double ry, int32_t hy, double dy,
    int32_t order, 
    ix_reduce_mode_t red,
    double z[]
  )
  {
    /* Slow but safe interpolation: */
    demand(hx >= 0, "invalid X grid size");
    demand(hy >= 0, "invalid Y grid size");
    int32_t jx, jy;
    int32_t jxy = 0;
    for (jy = -hy; jy <= +hy; jy++)
      { double y = ry + jy * dy;
        for (jx = -hx; jx <= +hx; jx++)
          { double x = rx + jx * dx;
            z[jxy] = float_image_interpolate_sample(A, c, x, y, order, red);
            jxy++;
          }
      }
  }
  
void float_image_interpolate_grid_pixels
  ( float_image_t *A, 
    double rx, int32_t hx, double dx,
    double ry, int32_t hy, double dy,
    int32_t order, 
    ix_reduce_mode_t red, 
    double z[]
  )  
  {
    /* Slow but safe interpolation: */
    demand(hx >= 0, "invalid X grid size");
    demand(hy >= 0, "invalid Y grid size");
    int32_t NC = A->sz[0];         /* Number of channels. */
    int32_t jx, jy;
    int32_t jxy = 0;
    for (jy = -hy; jy <= +hy; jy++)
      { double y = ry + jy * dy;
        for (jx = -hx; jx <= +hx; jx++)
          { double x = rx + jx * dx;
            float_image_interpolate_pixel(A, x, y, order, red, &(z[jxy]));
            jxy += NC;
          }
      }
  }

/* INTERNAL FUNCTIONS */
  
void float_image_interpolate_get_samples_and_weights
  ( float_image_t *A, 
    int32_t c, 
    double x, 
    double y,
    int32_t m, 
    ix_reduce_mode_t red, 
    float *p[],
    double wx[],
    double wy[]
  )
  {
    int32_t ix[m];
    float_image_interpolate_get_indices_and_weights(x, A->sz[1], m, red, ix, wx);
    
    int32_t iy[m];
    float_image_interpolate_get_indices_and_weights(y, A->sz[2], m, red, iy, wy);
      
    int32_t k = 0;
    int32_t dx, dy;
    for (dy = 0; dy < m; dy++)
      { int32_t iyk = iy[dy];
        for (dx = 0; dx < m; dx++)
          { int32_t ixk = ix[dx];
            p[k] = ((ixk < 0) || (iyk < 0) ? NULL : float_image_get_sample_address(A, c, ixk, iyk));
            k++;
          }
      }
  }
