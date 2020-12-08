/* Tools for images with attached masks. */
/* Last edited on 2013-10-21 00:20:38 by stolfilocal */

#ifndef float_image_masked_H
#define float_image_masked_H

#include <float_image.h>

typedef struct float_image_masked_t 
  { float_image_t *img; /* An image. */
    float_image_t *msk; /* A mask image (with values in [0_1]) for {img}. */
  } float_image_masked_t;
  /* The mask {msk} must be grayscale (sz[0] == 1)
    inependently of the number of channels of {img}. */

float_image_masked_t *float_image_masked_new(int nc, int nx, int ny);
  /* Creates an image pair {img,msk}, with {nx} columns,
    and {ny} rows, where {img} has {nc} channels and {msk} has 1 channel. */

void float_image_masked_free(float_image_masked_t *im);
  /* Frees the storage used by the image/mask pair (including the records {*im}). */

void float_image_masked_interpolate
  ( float_image_masked_t *im, 
    int c,
    double x,
    double y,
    int degInter,
    float *val,
    float *wht
  );
  /* Interpolates channel {c} of the masked image {im} at point {(x,y)},
    using an interpolator of degree {deginter}. Returns the interpolated
    value in {val} and the corresponding weight in {*wht}. */

/* AUXILIARY FUNCTIONS */

void interpolate_weighted_values(float v[], float w[], int degInter, double t, float *val, float *wht);

void interpolate_weighted_values_linear(float v[], float w[], double t, float *val, float *wht);

void interpolate_weighted_values_quadratic_B(float v[], float w[], double t, float *val, float *wht);

#endif
