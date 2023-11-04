#ifndef neuromat_median_filter_H
#define neuromat_median_filter_H

/* NeuroMat filtering and spectral analysis tools. */
/* Last edited on 2023-11-04 02:37:07 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <bool.h>

void neuromat_median_filter_apply
  ( int32_t nt,
    int32_t ne, 
    double **val, 
    double **med, 
    int32_t nw,
    bool_t verbose
  );
  /* Applies a running weighted median filter to each electrode channel {ie} in {0..ne-1} of the 
    discrete signal set {val[0..nt-1][0..ne-1]}.  Stores the result in {med[0..nt-1][0..ne-1]}.
    The arguments {val} and {med} may be the same array.  Any channels beyond the first {ne},
    in either array is ignored.
    
    Specifically, for each {ie} in {0..ne-1}},
    extracts the sample sequence {x[0..nt-1] = val[0..nt-1][ie]},
    computes a smoothed version {s[0..nt-1]} of it a median filter using a weighted
    window of width {nw}, and stores {s} into {med[0..nt-1][ie]}.
    
    The window width {nw} must be odd, {nw = 2*hw+1}, and at 
    least 3.
    
    The smoothed sample {s[it]} is computed by assigning a bell-like
    weight {wt[k]} to each sample {x[k]=val[it+k][ie]} for {k} in
    {-hw..+hw}, and then setting {s[i]} to the value {v} such that the
    total weight of the signal samples {x[-hw..+hw} which are less than
    {v} is half of the total weight.
    
    The window radius {hw} is reduced as needed when {it} is near 0 or {nt-1},
    so that {it-hw..it+hw} is always a subset of {0..nt-1}.  In particular,
    {s[0]} is just {x[0]} and {s[nt-1]} is just {x[nt-1]}. */

#endif
