/* Focus detector for multi-focus stereo. */
/* Last edited on 2017-12-27 16:15:43 by stolfilocal */

#ifndef multifok_focus_op_H
#define multifok_focus_op_H

#define _GNU_SOURCE
#include <stdint.h>

/* 
  FOCUS DETECTOR OPERATORS
  
  A focus detector is a local image operator that returns a measure of
  sharpness of the image at the applied pixel.  Ideally it is zero 
  for a totally unfocused (constant) image, and maximum for a 
  maximally focused one.
  
  The focus detector operators in this module work on a window of {NW}
  by {NW} samples, stored in an array {v[0..NS-1]} of {float}s, where {NS
  = NW*NW}, linearized by rows.
  
  The window size {NW} must be odd. In the current implementation, {NW}
  must be 3 */ 

double multifok_focus_op_score_simple(int32_t NW, double v[], double *phi[], double w[], double noise);
  /* Assumes that {v} contains the samples from a some channel of some image {img} contained
    in an {NW} by {NW} window, stored by rows.  Returns a numeric score that tells how 
    well-focused the image seems to be in that window.
    The samples are destroyed in the process. */ 

int32_t multifok_focus_op_num_samples(int32_t NW);
  /* Returns the number of pixels in a window of {NW} by {NW} pixels, 
    that is, {NW*NW}.  Also checks whether the window size is odd. */

double *multifok_focus_op_prod_weights(int32_t NW);
  /* Returns a newly allocated array {w[0..NS-1}} with the weights of each sample in the window,
    for inner product purposes. */

double multifok_focus_op_prod(int32_t NW, double x[], double y[], double w[]);
  /* Computes the inner product of sample vectors {x[0..NS-1]} and {y[0..NS-1]}
    with weights {w[0..NS-1]}.  That is, {SUM{i \in 0..NS-1 : w[i]*x[i]*y[i]}}. */

double multifok_focus_op_dist_sqr(int32_t NW, double x[], double y[], double w[]);
  /* Computes the squared distance of sample vectors {x[0..NS-1]} and {y[0..NS-1]}
    with weights {w[0..NS-1]}.  That is, {SUM{i \in 0..NS-1 : w[i]*(x[i] - y[i])^2}}. */

double **multifok_focus_op_basis(int32_t NW);
  /* Returns a newly allocated array {phi[0..NS-1][0..NS-1]} with the
    basis used to analyze the window samples. Namely {phi[i][j]} is
    coordnate {j} of basis element {i}. */

void multifok_focus_op_basis_check(int32_t NW, double *phi[], double w[]);
  /* Checks whether the basis {phi[0..NS-1][0..NS-1]} is othonormal
    under {multifok_focus_op_prod} with weights {w[0..NS-1]}. */

void multifok_focus_op_remap_samples(int32_t NW, double v[], double *phi[], double w[], double c[]);
  /* Applies a transformation of the window samples {v[0..NS-1]} to coefficients {c[0..NS-1]}
    using the basis {phi[0..NS-1][0..NS-1]} and inner product weights {w[0..NS-1]}. */

void multifok_focus_op_check(int32_t NW);
  /* Runs various diagnostics on this module. */

void multifok_focus_op_set_samples_3x3
  ( double x[], 
    double scale,
    float x00, float x01, float x02, 
    float x10, float x11, float x12, 
    float x20, float x21, float x22 
  );
  /* Specific for {NW = 3} ({NS = 9}): sets the elements {x[0..NS-1]} to the 
    values {x00,x01,x02,x10,...,x22} multiplied by {scale}. */

void multifok_focus_op_orthize(int32_t NW, int k, double *phi[], double w[]);
  /* Makes row {k} of the basis matrix {phi} orthogonal to all the previous ones and normalized. */

#endif
