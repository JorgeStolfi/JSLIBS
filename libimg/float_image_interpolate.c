/* See float_image_interpolate.h */
/* Last edited on 2025-01-14 16:25:01 by stolfi */ 

#include <stdint.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <ix_types.h>
#include <ix_reduce.h>
#include <float_image.h>

#include <float_image_interpolate.h>

/* INTERNAL PROTOTYPES */

int32_t float_image_interpolate_compute_num_samples(int32_t order);
  /* Computes the number of samples {m} needed along each axis for
    interpolation of the specified {order}. */

void float_image_interpolate_get_indices_and_weights
  ( double t, 
    int32_t N, 
    int32_t m, 
    ix_reduce_mode_t red,
    int32_t i[], 
    double w[]
  );
  /* Computes the indices {i[0..m-1]} of the {m} pixels needed to
    interpolate a row (or column) of samples at the fractional
    coordinate {t}, and their respective weights {w[0..m-1]}. 
    
    Indices {i} outside the range {0..N-1} are reduced according to
    {ix_reduce(i, N, red)}. The reduced index may be {-1} to denote `no
    such pixel'.  */

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
  /* Stores in {p[0..m*m-1]} the addresses of the {m�m} block samples
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
    
    /* Get pixel weights and addresses of samples in channel 0: */
    float *p[m*m];
    double wx[m];
    double wy[m];
    float_image_interpolate_get_samples_and_weights(A, 0, x, y, m, red, p, wx, wy);
    
    /* Apply interpolation formula to all channels. */
    int32_t NC = (int32_t)A->sz[0];         /* Number of channels. */
    int32_t cst = (int32_t)A->st[0];  /* Position increment between channels. */
    int32_t c;
    for (c = 0; c < NC; c++)
      { /* Apply interpolation formula to channel {c}, advance pointers to next channel: */
        double sumwp = 0;
        double sumw = 1.0e-200;
        int32_t kx, ky;
        for (ky = 0; ky < m; ky++)
          { float **py = &(p[m*ky]);
            double sumwp_y = 0;       /* Sum of pixels on window row, weighted by {wx}. */
            double sumw_y = 0; /* Sum of weights in {wx}. */
            for (kx = 0; kx < m; kx++)
              { float *pxy = py[kx];
                if (pxy != NULL) 
                  { double wxk = wx[kx];
                    sumwp_y += wxk*(*pxy);
                    sumw_y += wxk;
                    py[kx] += cst; }
              }
            double wyk = wy[ky];
            sumwp += wyk*sumwp_y;
            sumw += wyk*sumw_y;
          }
        z[c] = sumwp/sumw;
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
    int32_t NC = (int32_t)A->sz[0];         /* Number of channels. */
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
 
int32_t float_image_interpolate_compute_num_samples(int32_t order)
  {
    /* !!! Rethink !!! */
    if (order == -1)
      { return 1; }
    else if (order == 0)
      { return 2; }
    else if (order == 1)
      { return 4; }
    else
      { demand(FALSE, "invalid interpolation order"); }
  }
 
void float_image_interpolate_get_indices_and_weights
  ( double t, 
    int32_t N, 
    int32_t m, 
    ix_reduce_mode_t red,
    int32_t i[], 
    double w[]
  )
  {
    bool_t debug = FALSE;

    /* Adjust {t} to be relative to start of first source pixel: */
    t = t - 0.5*(m - 1);
    
    /* Get the raw index of the first source pixel: */
    int32_t iz = (int32_t)floor(t);
    
    /* Get the fraction {fz} fromdata pixel 0: */
    double fz = t - iz;

    /* Make sure that the fraction is in {[0 _ 1]}: */
    if (fz < 0.0) { iz--; fz += 1.0; }
    if (fz > 1.0) { iz++; fz -= 1.0; }
    
    /* Compute the fraction's complement {gz}: */
    double gz = 1.0 - fz;
    
    /* Compute the reduced indices {i[0..m-1]}: */
    for (int32_t k = 0; k < m; k++)
      { i[k] = (int32_t)ix_reduce(iz, (ix_size_t)N, red);
        iz++;
      }

    /* Compute the interpolation weights {w[0..m-1]}: */
    if (m == 1)
      { /* Piecewise-constant interpolation: */
        w[0] = 1.0;
      }
    else if (m == 2)
      { /* Affine interpolation: */
        w[0] = gz;
        w[1] = fz;
      }
    else if (m == 4)
      { /* Cubic C1 interpolant: */
        double a = 2.0/3.0;
        double b = 7.0/3.0;
        double c = b-1;
        w[0] = -a*gz*gz*fz; 
        w[1] = 1.0 - fz*fz*(b - c*fz);
        w[2] = 1.0 - gz*gz*(b - c*gz);
        w[3] = -a*fz*fz*gz;
      }

    if (debug)
      { fprintf(stderr, "t = %7.4f", t);
        for (int32_t k = 0; k < m; k++)
          { fprintf(stderr, "  i[%d] = %d  w[%d] = %10.7f\n", k, i[k], k,w[k]); }
        fprintf(stderr, "\n");
      }
  }
  
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
    float_image_interpolate_get_indices_and_weights(x, (int32_t)A->sz[1], m, red, ix, wx);
    
    int32_t iy[m];
    float_image_interpolate_get_indices_and_weights(y, (int32_t)A->sz[2], m, red, iy, wy);
      
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
