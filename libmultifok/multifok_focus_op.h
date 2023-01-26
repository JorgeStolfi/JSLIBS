/* Focus detector for multi-focus stereo. */
/* Last edited on 2023-01-25 13:54:57 by stolfi */

#ifndef multifok_focus_op_H
#define multifok_focus_op_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <i2.h>

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

void multifok_focus_op_normalize_window_samples
  ( int32_t NW, 
    double x[], 
    double ws[], 
    double noise, 
    double *avg_P,
    double *dev_P
  );
  /* Computes the weighted average {avg} and deviation {dev} of the window samples {x[0..NS-1]},
    returning them in {*ave_P} and {*dev_P}.
    
    Then normalizes the samples {x[0..NS-1]} by subtracting {avg}
    dividing by {hypot(dev, noise}, so that they have zero mean and unit deviation.
    
    This correction elimines the effect of brightness and contrast variations,
    assuming that the samples are contaminated with random noise with mean 0 and deviation {noise},
    so variations that are small compared to {noise} are not amplified. */

/* FOCUS SCORE DERIVED FROM LOCAL BASIS  */
  
typedef enum {
    multifok_focus_op_basis_type_CANC, /* The canonical basis, 1 at each sample and 0 at others. */
    multifok_focus_op_basis_type_LAPL, /* The order 2 differential operators: "DX", "DY", "DXX", "DXY", "DYY". */
    multifok_focus_op_basis_type_CUBE, /* Same as {LAPL} plus the 3rd order saddles "C3", "S3". */
    multifok_focus_op_basis_type_DIFF, /* The {CUBE} operators plus the checker "Q". */
    multifok_focus_op_basis_type_HART  /* The Hartley basis waves. */
  } multifok_focus_op_basis_type_t;
  /* Basis types. */
  
#define multifok_focus_op_basis_type_FIRST multifok_focus_op_basis_type_CANC
#define multifok_focus_op_basis_type_LAST multifok_focus_op_basis_type_HART 

void multifok_focus_op_basis_make
  ( int32_t NW, 
    double ws[],
    multifok_focus_op_basis_type_t bType,
    bool_t ortho,
    int32_t *NB_P, 
    double ***bas_P, 
    char ***belName_P
  );
  /* Returns in {*bas_P} a newly allocated array {bas[0..NB-1][0..NS-1]} with a
    basis to be used to analyze the window samples. Namely {bas[i][j]}
    is coordnate {j} of basis element {i}. The number {NB} is chosen by
    the procedure and returned in {*NB_P}.
    
    If {bType} is {multifok_focus_op_basis_type_CANC}, basis element
    {kb} will be 0 except at sample {kb}, where it will be 1.
    
    If {bType} is {multifok_focus_op_basis_type_LAPL}, then {NW} must be
    3, and the basis will consist of the five differential operators "DX", "DY",
    "DXX", "DXY", "DYY".
    
    If {bType} is {multifok_focus_op_basis_type_DIFF}, then {NW} must be
    3, and the basis will consist of the same operators as {LAPL},
    as well as the local image value "F", the two 3-saddles, and the 3x3 checker.
    
    If {bType} is {multifok_focus_op_basis_type_HART} is false, the
    basis will consists of the Hartley waves whose frequency vectors are
    the coordinates of the samples in the window, relative to the center
    pixel.
    
    If {ortho} is true, the basis will be modified to be orthonormal
    with respect to the the dot product {multifok_focus_op_prod}, with
    sample weights {ws[0..NS-1]}; that is, the dot product of
    {bas[ib][0..NS-1]} and {bas[jb][0..NS-1]} is 1 if {ib=jb} and 0
    otherwise.  During this procedure, some elements may be discarded
    for being too close to linear combinations of the previous ones.
    
    If {ortho} is false, each basis element will be scaled so that 
    the dot product with any window whose samples are in {[0_1]} will
    be in the range {[-1 _ +1]} and hit at least one of the two 
    bounds. 
    
    The procedure also returns in {*belName_P} a vector of {NB} newly alocated strings 
    such that {belName[k]} is a scrutable name for basis element {bas[k]}. */
 
void multifok_focus_op_compute_basis_coeffs(int32_t NW, double x[], double ws[], int32_t NB, double *bas[], double coeff[]);
  /* Applies a transformation of the window samples {x[0..NS-1]} to coefficients {coeff[0..NB-1]}
    using the basis {bas[0..NB-1][0..NS-1]} and inner product weights {ws[0..NS-1]}.
    
    The basis is assumed to be orthonormal with respect to the dot product with
    sample weights {ws}, so that {coeff[kb]} is the weighted inner product of {bas[kb]}. */

char *multifok_focus_op_basis_type_to_text(multifok_focus_op_basis_type_t bType);
  /* Converts the basis type {bType} to text ("LAPL", "DIFF", etc.). 
    Clients should NOT call {free} on the returned strings. */
    
multifok_focus_op_basis_type_t multifok_focus_op_basis_type_from_text(char *name, bool_t fail);
  /* Converts the {name} ("LAPL", "DIFF", etc.) to a basis type.  
    If the {name} invalid, fails if {fail} is true, and returns {-1} if false. */

void multifok_focus_op_basis_free(int32_t NB, double **bas, char **belName);
  /* Frees the storage used by {bas[0..NB-1][0..NS-1]} and {belName[0..NB-1]}. */

void multifok_focus_op_term_indices_from_names(int32_t NB, char *belName[], int32_t NT, char *tname[], i2_t tix[]);
  /* Sets {tix[0..NT-1]} to the term index pairs implied by the term names {tname[0..NT-1}} and the 
    basis element names {belName[0..NB-1]}. See {multifok_focus_op_score_from_basis_coeffs}. */

double multifok_focus_op_score_from_basis_coeffs
  ( int32_t NB, 
    double coeff[],
    int32_t NT,
    i2_t tix[],
    double wt[],
    double term[], /* OUT */
    bool_t squared
  );
  /* Assumes that {coeff[0..NB-1]} are the coefficients of the
    window samples in some local operator basis.
    
    Then, if {NT} is positive, the procedure forms {NT} term values {term[0..NT-1]} from the
    coeefs {coeff[0..NB-1]} as specified in the vector {tix[0..NT-1]}.
    Namely, let {tix[k]} be {(ib,jb)}. If {ib} and {jb} are
    non-negative, then {term[k]} is the quadratic term {coeff[ib]*coeff[jb]}. If
    only {ib} is non-negative, then {term[k]} is {coeff[ib]}. If both are
    negative, then {term[k]} is 1.
    
    Finally, the procedure computes the sharpness score as the weighted sum 
    {SUM{wt[k]*term[k] : k \in 0..NT-1 }}.  However, if {squared} is true,
    assumes that the combination of the terms is an estimate
    of the squared sharpness, so that the returned {score} is the square root 
    of that combination. In any case the returned value is clipped to the
    range {[0 _ 1]}.
   
    If {NT} is zero, the parameters {tix} and {wt} are ignored, and the
    sharpness score is set to zero.
    
    The samples are destroyed in the process. */ 

void multifok_focus_op_basis_print
  ( FILE *wr,
    int32_t NW, 
    int32_t NB, 
    double *bas[],
    char *belName[]
  );
  /* Prints to {wr} the basis {bas[0..NB-1][0..NS-1]}. */

void multifok_focus_op_basis_ortho_check
  ( FILE *wr,
    int32_t NW, 
    double ws[], 
    int32_t NB, 
    double *bas[]
  );
  /* Checks whether basis {bas[0..NB-1][0..NS-1]} is osrthonormal. */

void multifok_focus_op_check(int32_t NW, multifok_focus_op_basis_type_t bType, bool_t ortho);
  /* Runs various diagnostics on this module. */

#endif
