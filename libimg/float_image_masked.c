/* See {float_image_masked.h} */
/* Last edited on 2013-10-21 00:22:15 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <float_image.h>
#include <float_image_masked.h>

void float_image_masked_interpolate
  ( float_image_masked_t *im, 
    int c,
    double x,
    double y,
    int degInter,
    float *val,
    float *wht
  )
  {
    assert(im->msk->sz[0] == 1);
    int nx = (int)im->img->sz[1];
    int ny = (int)im->img->sz[2];

    int m = degInter+1; /* Number of data points needed along each axis. */

    /* The inpterpolation acts on a window of size {m} by {m} approximately centered at {x,y}: */
    double delta = ((double)degInter)/2;
    x -= delta;
    y -= delta;
    int ix = (int)(floor(x)); /* First column of window. */
    int iy = (int)(floor(y)); /* First row of window. */
    double fx = x - ix; /* Interpolation argument in X. */
    double fy = y - iy; /* Interpolation argument in Y. */
    float vc[m]; /* Interpolated values for each window column. */
    float wc[m]; /* Weights of those values. */
    int dx, dy;
    for (dx = 0; dx < m; dx++) {
      float v[m]; /* Image values along window column {dx}. */
      float w[m]; /* Weights of those values. */
      for (dy = 0; dy < m; dy++) {
        int jx = ix + dx;
        int jy = iy + dy;
        if ((jx >= 0) && (jx < nx) && (jy >= 0) && (jy < ny)) {
          v[dy] = float_image_get_sample(im->img, c, jx, jy);
          w[dy] = float_image_get_sample(im->msk, 0, jx, jy);
        } else {
          v[dy] = w[dy] = 0;
        }
      }
      /* Interpolate vertically along column {dx}: */
      interpolate_weighted_values(v, w, degInter, fy, &(vc[dx]), &(wc[dx]));
    }
    /* Interpolate horizontally: */
    interpolate_weighted_values(vc, wc, degInter, fx, val, wht);
  }

void interpolate_weighted_values(float v[], float w[], int degInter, double t, float *val, float *wht)
  {
    switch(degInter) {
    case 0:
      { (*val) = v[0]; (*wht) = w[0]; }
      break;
    case 1:
      { interpolate_weighted_values_linear(v, w, t, val, wht); }
      break;
    case 2:
      { interpolate_weighted_values_quadratic_B(v, w, t, val, wht); }
      break;
    default:
      assert(FALSE);
    }
  }

void interpolate_weighted_values_linear(float v[], float w[], double t, float *val, float *wht)
  {
    /* C0, interpolating, non-negative, unit sum. */
    double f0 = 1 - t;
    double f1 = t;
    (*val) = (float)(f0*v[0] + f1*v[1]);
    (*wht) = (float)(2.0/(1.0/w[0] + 1.0/w[1]));
  }

void interpolate_weighted_values_quadratic_B(float v[], float w[], double t, float *val, float *wht)
  {
    /* C1, non-interpolating, non-negative, unit sum. */
    double v0 = v[0], v1 = v[1], v2 = v[2];
    double w0 = w[0], w1 = w[1], w2 = w[2];
    if ((w0 != 0) && (w1 != 0) && (w2 != 0)) {
      /* Nothing to fake. */
    } else if ((w0 != 0) && (w1 != 0)) {
      /* Extrapolate {v2} from {v0,v1}: */
      v2 = v1 + (v1 - v0);
      w2 = (w0 + w1)/4;
    } else if ((w1 != 0) && (w2 != 0)) {
      /* Extrapolate {v0} from {v1,v2}: */
      v0 = v1 + (v1 - v2);
      w0 = (w1 + w2)/4;
    } else if ((w0 != 0) && (w2 != 0)) {
      /* Interpolate {v1} from {v0,v2}: */
      v1 = (v0 - v2)/2;
      w1 = (w0 + w2)/4;
    } else if (w0 != 0) {
      /* Copy {v0} to {v1,v2}: */
      v1 = v2 = v0;
      w1 = w0/2; w2 = w0/4;
    } else if (w1 != 0) {
      /* Copy {v1} to {v0,v2}: */
      v0 = v2 = v1;
      w0 = w2 = w1/2;
    } else if (w2 != 0) {
      /* Copy {v2} to {v0,v1}: */
      v0 = v1 = v2;
      w0 = w2/4; w1 = w2/2;
    } else {
      /* Give up: */
    }

    /* Now apply biquadratic interpolation to {v0,v1,v2}: */
    double f0 = (1 - t)*(1 - t)/2.0;
    double f1 = (0.5 + t*(1 - t));
    double f2 = t*t/2.0;
    (*val) = (float)(f0*v0 + f1*v1 + f2*v2);

    /* Take a weighted harmonic mean of the weights: */
    (*wht) = (float)(4.0/(1.0/w0 + 2.0/w1 + 1.0/w2));
  }

float_image_masked_t *float_image_masked_new(int nc, int nx, int ny)
  {
    float_image_masked_t *im = malloc(sizeof(float_image_masked_t));
    (*im) = (float_image_masked_t) { 
      .img = float_image_new(nc, nx, ny), 
      .msk = float_image_new(1, nx, ny) 
    };
    return im;
  }

void float_image_masked_free(float_image_masked_t *im)
  {
    float_image_free(im->img); 
    float_image_free(im->msk);
    free(im);
  }
