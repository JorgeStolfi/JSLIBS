#ifndef conv_filter_H
#define conv_filter_H

/* Convolution and downsampling of a sequence with a filter kernel. */
/* Last edited on 2024-11-15 19:12:13 by stolfi */

#define conv_filter_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#include <stdint.h>

void conv_filter(int32_t nx, double x[], int32_t skip, int32_t step, int32_t nw, double w[], int32_t ny, double y[]);
  /* Computes the convolution of the sample sequence {x[0..nx-1]} with the weight
    sequence {w[0..nw-1]}, and stores the result in {y[0..ny-1]}, downsampled with skip
    {skip} and step {step}.  The weight table must have odd length.
    
    More precisely, sets each {y[i]} to
    
      {SUM{ x[skip+step*i+(j-hw)] * w[j]}} / {SUM{w[j]}}
      
    where {hw = (nw-1)/2}, and the sum ranges over all {j} such that all indices are valid. */

#endif
