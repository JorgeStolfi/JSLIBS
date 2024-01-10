#ifndef conv_filter_H
#define conv_filter_H

/* Convolution and downsampling of a sequence with a filter kernel. */
/* Last edited on 2014-07-26 23:04:31 by stolfilocal */

#define conv_filter_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE

void conv_filter(int nx, double x[], int skip, int step, int nw, double w[], int ny, double y[]);
  /* Computes the convolution of the sample sequence {x[0..nx-1]} with the weight
    sequence {w[0..nw-1]}, and stores the result in {y[0..ny-1]}, downsampled with skip
    {skip} and step {step}.  The weight table must have odd length.
    
    More precisely, sets each {y[i]} to
    
      {SUM{ x[skip+step*i+(j-hw)] * w[j]}} / {SUM{w[j]}}
      
    where {hw = (nw-1)/2}, and the sum ranges over all {j} such that all indices are valid. */

#endif
