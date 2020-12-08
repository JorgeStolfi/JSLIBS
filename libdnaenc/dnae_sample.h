#ifndef dnae_sample_H
#define dnae_sample_H

/* Encoded of numerical signal samples */
/* Last edited on 2014-06-14 01:42:04 by stolfilocal */

#define dnae_sample_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdint.h>

#include <vec.h>

/* This interface defines the type {dnae_sample_enc_t}, an encoded 
  sample from a real-valued signal. 
  
  The numbers are encoded as signed short ints in the symmstric range {-M..+M}, 
  where {M = 2^15-1}, using a non-linear encoding. The encoding is designed to map
  a normal-distributed variable to a uniformly distributed code 
  in that range. */
  
typedef int16_t dnae_sample_enc_t; 
  /* One coordinate of a sample, packed as a character. */

#define dnae_sample_enc_VALID_MIN (-32767)
#define dnae_sample_enc_VALID_MAX (+32767)
  /* Minimum and maximum valid values of a {dnae_sample_enc_t}. */

double dnae_sample_decode(dnae_sample_enc_t ev, double scale);
dnae_sample_enc_t dnae_sample_encode(double dv, double scale);
  /* Converts between unpacked ({double}) and packed ({dnae_sample_enc_t}) 
    representations of sample values, asuming the scale factor {scale}. */

double dnae_sample_diffsq(dnae_sample_enc_t xs, double xscale, dnae_sample_enc_t ys, double yscale);
  /* Returns the abs difference squared between the samples,
    defined as {(dnae_sample_decode(xs, xscale) - dnae_sample_decode(ys, yscale))^2}. */

vec_typedef(dnae_sample_enc_vec_t,dnae_sample_enc_vec,dnae_sample_enc_t);
/* Vector of {dnae_sample_enc_t}. */

#endif

