/* Window operations for multi-focus stereo. */
/* Last edited on 2023-01-30 06:56:05 by stolfi */

#ifndef multifok_window_H
#define multifok_window_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

/* PROCESSING WINDOWS
  
  The procedures in the {multifok} library work on windows of {NW} by
  {NW} samples from a single channel of an image. 
  
  The window size {NW} must be odd and at least 3.
  
  The samples in a window are usually stored in an array {s[0..NS-1]} of
  {double}s, where {NS = NW*NW}, linearized by rows. That is, when the
  window is centered at a pixel in column {ix} and row {iy} the sample
  {s[jy*NW+jx]} is taken from image pixel in column {ix+(jx-HW)} and row
  {iy+(jy-HW)}, where {HW = (NW-1)/2}. */ 

int32_t multifok_window_num_samples(int32_t NW);
  /* Returns the number of pixels in a window of {NW} by {NW} pixels, 
    that is, {NW*NW}.  Also checks whether {NW} is odd and at least 3. */

/* WINDOWING WEIGHTS

  The effective window usually has smooth boundaries and apodizing.
  Its effective shape is defined by a list of {NS} non-negative weights, 
  stored in a vector, linearized by rows. */

double *multifok_window_sample_weights(int32_t NW);
  /* Returns a newly allocated array {ws[0..NS-1}} with the weights of
    each sample in the window, for averaging and apodizing purposes.
    
    Currently uses binomial weights {ws[ks] = choose(NW-1,ix)*choose(NW-1,iy)/A}
    where {ix,iy} vary in {0..NW-1}, {ks} is {iy*NW + ix}, and {A} is such that the
    central weight is 1. */

double multifok_window_prod(int32_t NW, double a[], double b[]);
  /* Computes the inner product of sample vectors {a[0..NS-1]} and {b[0..NS-1]}. 
    That is, {SUM{i \in 0..NS-1 : a[i]*b[i]}}. */

double multifok_window_dist_sqr(int32_t NW, double a[], double b[]);
  /* Computes the squared distance of sample vectors {a[0..NS-1]} and {b[0..NS-1]}. 
    That is, {SUM{i \in 0..NS-1 : (a[i] - b[i])^2}}. */

void multifok_window_normalize_samples
  ( int32_t NW, 
    double s[], 
    double ws[], 
    double noise, 
    double *avg_P,
    double *dev_P
  );
  /* Computes the weighted average {avg} and deviation {dev} of the window samples {s[0..NS-1]},
    with sample weights {ws[0..NS-1]}, returning them in {*ave_P} and {*dev_P}.
    
    Then normalizes the samples {s[0..NS-1]} by subtracting {avg}
    dividing by {hypot(dev, noise}, so that they have zero mean and unit deviation.
    
    This correction elimines the effect of brightness and contrast
    variations, assuming that the samples are contaminated with random
    noise with mean 0 and deviation {noise}. Variations that are
    small compared to {noise} are not amplified. */

void multifok_window_set_samples_3x3
  ( double s[], 
    double scale,
    double s00, double s01, double s02, 
    double s10, double s11, double s12, 
    double s20, double s21, double s22 
  );
  /* Specific for {NW = 3} ({NS = 9}): sets the elements {s[0..NS-1]} to the 
    values {s00,s01,s02,s10,...,s22} multiplied by {scale}. */

/* APLHA SAMPLE NAMES
  
  Samples are sometimes refrerred by indices relative to the center of the window,
  positive, zero, or negative.  The procedures below provide a way to encode such
  indices with letters and digits only without using '+' and '-'.  */

char *multifok_window_mop_code(int32_t i);
  /* Returns the index {i} encoded as a newly allocated string. The 
    string is "o" if {i} is zero, "m" if {i} is {-1}, and "p" if {i} is {+1}.
    Otherwise the result is "m{N}" if {i} is negative, of "p{N}" if
    {i} is positive; where {N} is the absolute value of {i}. */

char *multifok_window_sample_name(char *tag, int32_t ix, int32_t iy);
  /* Returns a newly allocated string with the index pair {ix,iy}
    encoded as "{tag}{X}{Y}", where {X} and {Y} are the indices {ix}
    and {iy} encoded with {multifok_window_mop_code}. */

void multifok_window_sample_names(int32_t NW, char *tag, char *sname[]);
  /* Sets {sname[0..NS-1]} to {NS = NW*NW} newly allocated strings that are scrutable names of
    the samples in the window. The names have the form "{tag}{X}{Y}" where  {X} and {Y} are the indices {ix}
    and {iy} of the sample relative to the center sample, encoded with {multifok_window_mop_code}. */

#endif
