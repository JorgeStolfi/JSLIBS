#ifndef neuromat_median_filter_H
#define neuromat_median_filter_H

/* NeuroMat filtering and spectral analysis tools. */
/* Last edited on 2023-11-28 04:19:36 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <affirm.h>
#include <wt_table.h>
#include <wt_table_generic.h>

#include <wt_median.h>
#include <wt_median_window.h>

void neuromat_median_filter_apply
  ( int32_t nt,
    int32_t ne, 
    double **val, 
    double **med, 
    int32_t nw,
    int32_t wt[],
    bool_t verbose
  );
  /* Applies a running weighted median filter to each electrode channel
    {ie} in {0..ne-1} of the discrete signal set
    {val[0..nt-1][0..ne-1]}. Stores the result in
    {med[0..nt-1][0..ne-1]}. The arguments {val} and {med} may be the
    same array. Any channels beyond the first {ne}, in either array is
    ignored.
    
    Specifically, for each {ie} in {0..ne-1}}, extracts the sample
    sequence {x[0..nt-1] = val[0..nt-1][ie]}, computes a smoothed
    version {s[0..nt-1]} of it a median filter using a weighted window
    of width {nw}, and stores {s} into {med[0..nt-1][ie]}.
    
    The window width {nw} must be odd, {nw = 2*hw+1}, and at 
    least 3.
    
    The smoothed sample {s[it]} is computed by assigning a bell-like
    weight {wt[iw]} to each sample {v[iw]=val[it+iw-hw][ie]} for {iw} in
    {0..nw-1}, and then setting {s[it]} to the value {vm} such that the
    total weight of the signal samples {v[0..nw-1} which are less than
    {vm} is half of the total weight. 
    
    Near the ends of the dataset the window extends into indices that are
    negative or beyond {nt-1}.  In those places, the channels samples are 
    implicitly extrapolated by mirroring. */

#endif
