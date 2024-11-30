#ifndef conv_filter_H
#define conv_filter_H

/* Convolution and downsampling of a sequence with a filter kernel. */
/* Last edited on 2024-11-23 05:36:48 by stolfi */

#define conv_filter_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#include <stdint.h>

#include <ix_reduce.h>

void conv_filter
  ( uint64_t nx, double x[],
    ix_reduce_mode_t ixred,
    int64_t skip,
    int64_t step,
    uint64_t nw, double w[],
    uint64_t ny, double y[]
  );
  /* Computes the convolution of the sample sequence {x[0..nx-1]} with
    the weight sequence {w[0..nw-1]}, and stores the result in
    {y[0..ny-1]}, downsampled with skip {skip} and step {step}. The
    weight table length {nw} must be odd (hence positive).
    
    More precisely, sets each {y[i]} to
    
      {SUM_j{x[skip+step*i+(j-hw)] * w[j]}} / {SUM_j{w[j]}}
      
    where {hw = (nw-1)/2}, and the sums range over all {j} such that all
    indices are valid.
    
    Note that the central weight of the window is aligned
    with sample {x[skip]} when computing {y[0]}, with {x[skip+step]}
    when computing {y[1]}, with {x[skip+2*step]} when computing {y[2]},
    and so on.  
    
    In any case, the value of {y[i]} in principle uses {nw} consecutive
    samples of {x}, centered at {x[skip+i*step]}. The parameter {ixred}
    specifies what happens when a sample index {r = skip+step*i+(j-hw)}
    falls outside the range {0..nx-1}. If {ixred} is
    {ix_reduce_mode_SINGLE} that sample and its weight are excluded from
    both sums. Otherwise the sample index {r} is remapped to the range
    {0..nx-1} with {ix_reduce} (quod videt).
    
    In any case, if the sum of weights in the denominator is zero, then
    {y[i]} may be set to {±INF} or {NAN}.  In particular, if {ixred} is
    {ix_reduce_mode_SINGLE} and the enrire window for {y[i]} falls outside
    {0..nx-1} (so that both summations are empty), then {y[i]} is set to
    {NAN}.
    
    The subsampling stride {step} may be negative or zero, although the
    latter generally does not make much sense as it results in {y[i]}
    being set to the same value for all{i}. The {skip} and {step} must
    be such that the index formula {skip+step*i+(j-hw)} does not
    overflow the {int64_t} range for any {i} in {0..ny-1]} and any {j}
    in {0..nw-1}.*/
    

#endif
