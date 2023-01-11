/* Focus detector for multi-focus stereo. */
/* Last edited on 2023-01-09 10:21:02 by stolfi */

#ifndef multifok_focus_op_H
#define multifok_focus_op_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

/* 
  FOCUS DETECTOR OPERATORS
  
  A focus detector is a local image operator that returns a measure of
  sharpness of the image at the applied pixel.  Ideally it is zero 
  for a totally unfocused (constant) image, and maximum for a 
  maximally focused one.
  
  The focus detector operators in this module work on a window of {NW}
  by {NW} samples, stored in an array {x[0..NS-1]} of {float}s, where {NS
  = NW*NW}, linearized by rows.
  
  The window size {NW} must be odd and at least 3. */ 
  
typedef enum {
    multifok_focus_op_basis_type_NONE,
    multifok_focus_op_basis_type_DIFF,
    multifok_focus_op_basis_type_HART
  } multifok_focus_op_basis_type_t;
  /* Basis types. */

int32_t multifok_focus_op_num_samples(int32_t NW);
  /* Returns the number of pixels in a window of {NW} by {NW} pixels, 
    that is, {NW*NW}.  Also checks whether {NW} is odd and at least 3. */

/* WINDOWING WEIGHTS

  The effective window usually has smooth boundaries and apodizing.
  Its effective shape is defined by a list of {NS} non-negative weights, 
  stored in a vector, linearized by rows. */

double *multifok_focus_op_sample_weights(int32_t NW);
  /* Returns a newly allocated array {ws[0..NS-1}} with the weights of each sample in the window,
    for inner product purposes. */

double multifok_focus_op_prod(int32_t NW, double x[], double y[], double ws[]);
  /* Computes the inner product of sample vectors {x[0..NS-1]} and {y[0..NS-1]}
    with weights {ws[0..NS-1]}.  That is, {SUM{i \in 0..NS-1 : ws[i]*x[i]*y[i]}}. */

double multifok_focus_op_dist_sqr(int32_t NW, double x[], double y[], double ws[]);
  /* Computes the squared distance of sample vectors {x[0..NS-1]} and {y[0..NS-1]}
    with weights {ws[0..NS-1]}.  That is, {SUM{i \in 0..NS-1 : ws[i]*(x[i] - y[i])^2}}. */

void multifok_focus_op_normalize_samples(int32_t NW, double x[], double ws[], double noise);
  /* Normalizes the samples {x[0..NS-1]} to have weighted average 0 and weighted deviation 1.
    This correction elimines the effect of brightness and contrast variations.
    
    Assumes that the samples are contaminated with random noise with mean 0 and deviation {noise},
    so variations that are small compared to {noise} are not amplified. */

/* SIMPLE FOCUS SCORE */

/* FOCUS SCORE DERIVED FROM LOCAL BASIS  */

double multifok_focus_op_score_from_basis(int32_t NW, double x[], double noise, double ws[], int32_t NB, double *phi[], double wc[]);
  /* Assumes that {x[0..NS-1]} contains the samples from a some channel
    of some image {img} contained in an {NW} by {NW} window, stored by
    rows, where {NS = NW*NW}. Returns a numeric score that tells how
    well-focused the image seems to be in that window.
    
    Before computing the score, the sampels {x[0..NS-1]} are 
    normalized by shifting and scaling with {multifok_focus_op_normalize_samples},
    so that they have zero mean and unit deviation.
    
    The score is then obtained by computing the coeffs {c[0..NB-1]} of the
    window samples in the given basis {phi[0..NB-1][0..NS-1]}, then
    computing the weighted sum of squares of the coeffs, namely {SUM{
    wc[i]*c[i]^2 : i \in 0..NB-1 }}.
    
    The basis is assumed to be orthonormal with respect to the dot
    product with sample weights {ws}.
    
    The samples are destroyed in the process. */ 

void multifok_focus_op_basis_make
  ( int32_t NW, 
    double ws[],
    multifok_focus_op_basis_type_t bType,
    bool_t ortho,
    int32_t *NB_P, 
    double ***phi_P, 
    double **wc_P,
    char ***name_P
  );
  /* Returns in {*phi_P} a newly allocated array {phi[0..NB-1][0..NS-1]} with a
    basis to be used to analyze the window samples. Namely {phi[i][j]}
    is coordnate {j} of basis element {i}. The number {NB} is chosen by
    the procedure and returned in {*NB_P}.
    
    If {bType} is {multifok_focus_op_basis_type_NONE}, basis element
    {kb} will be 0 except at sample {kb}, where it will be 1.
    
    If {bType} is {multifok_focus_op_basis_type_DIFF}, then {NW} must be
    3, and the basis will consist of all the differential operators of
    order 0, 1, and 2, as well as two 3-saddles and the 3x3 checker.
    
    If {bType} is {multifok_focus_op_basis_type_HART} is false, the
    basis will consists of the Hartley waves whose frequency vectors are
    the coordinates of the samples in the window, relative to the center
    pixel.
    
    If {ortho} is true, the basis will be modified to be orthonormal
    with respect to the the dot product {multifok_focus_op_prod}, with
    sample weights {ws[0..NS-1]}; that is, the dot product of
    {phi[ib][0..NS-1]} and {phi[jb][0..NS-1]} is 1 if {ib=jb} and 0
    otherwise.  During this procedure, some elements may be discarded
    for being too close to linear combinations of the previous ones.
    
    If {ortho} is false, each basis element will be scaled so that 
    the dot product with any window whose samples are in {[0_1]} will
    be in the range {[-1 _ +1]} and hit at least one of the two 
    bounds. 
    
    The procedure also returns in {*wc_P} a vector of {NB} elements
    where {wc[k]} is the weight to be used when combining the squared 
    coefficient of {phi[k]} into the focus indicator.
    
    The procedure also returns in {*name_P} a vector of {NB} newly alocated strings 
    such that {name[k]} is a scrutable name for basis element {phi[k]}. */

void multifok_focus_op_basis_free(int32_t NB, double **phi, double *wc, char **name);
  /* Frees the storage used by {phi[0..NB-1][0..NS-1]}, {wc[0..NB-1]}, and {name[0..NB-1]}. */

void multifok_focus_op_basis_print
  ( FILE *wr,
    int32_t NW, 
    int32_t NB, 
    double *phi[],
    double wc[],
    char *name[]
  );
  /* Prints to {wr} the basis {phi[0..NB-1][0..NS-1]} and the coefficient weights {wc[0..NB-1]}. */

void multifok_focus_op_basis_ortho_check
  ( FILE *wr,
    int32_t NW, 
    double ws[], 
    int32_t NB, 
    double *phi[]
  );
  /* Checks whether basis {phi[0..NB-1][0..NS-1]} is osrthonormal. */

void multifok_focus_op_remap_samples(int32_t NW, double x[], double ws[], int32_t NB, double *phi[], double c[]);
  /* Applies a transformation of the window samples {x[0..NS-1]} to coefficients {c[0..NB-1]}
    using the basis {phi[0..NB-1][0..NS-1]} and inner product weights {ws[0..NS-1]}. */

void multifok_focus_op_check(int32_t NW, multifok_focus_op_basis_type_t bType, bool_t ortho);
  /* Runs various diagnostics on this module. */

#endif
