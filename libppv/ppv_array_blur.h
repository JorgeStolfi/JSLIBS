/* ppv_array_blur.h --- blurring (local averaging) of a PPV array. */
/* Last edited on 2021-07-03 13:57:19 by jstolfi */

#ifndef ppv_array_blur_H
#define ppv_array_blur_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <ppv_types.h>
#include <ppv_array.h>

ppv_array_t *ppv_array_blur
  ( ppv_array_t *A, 
    ppv_sample_floatize_proc_t *floatize,
    ppv_nbits_t bps, 
    ppv_size_t radius,
    double wt[],
    ppv_size_t hstep,
    ppv_sample_quantize_proc_t *quantize
  );
  /* Computes the local averaged version of the array {A}.
    
    The result is an array {G} with same number of dimensions {d} as
    {A}, but {bps} bits per sample instead of {A.bps}.  Bot of these
    two parameters must be positive.
    
    The {radius} parameter defines the size of the neighborhood
    used in the averaging.  It must be at least 1.  
    The value of each voxel of {G} will be the weighted average of the
    voxels of {A} around some voxel {A[ixA]}, using a Hahn (shifted cosine)
    windowing weight mask that extends {radius} voxels in every
    direction from {A[ixA]}; that is, a total of {t = 2*radius+1} voxels
    along each axis.
    
    The {hstep} parameter defines the spacing of the probe points.   
    Specifically, let {s} be {2*hstep+1}. The size of{G} along
    each axis {k} will be {A.size[k]/s}, and each voxel {G[ixG]} is the
    weighted average of neighborhood of {A[ixA]} where {ixA[k] = s*ixG[k] +
    hstep} for each axis {k}. 
    
    The array {wt} should have {t} elements. The weight of a 
    sample {A[jxA]} in the neighborhood of a voxel {A[ixA]} will be the product
    of {wt[jxA[k]-ixA[k]+radius]} for all axes {k}.
    
    The {hstep} value must be non-negative and must not exceed {radius}.
    In addition, {s = hstep+1} must be a divisor of {t = 2*radius+1}.
    In particular, the values {0} and {radius} are always
    valid. If {hstep = 0}, {G} will have the same size as {A}. 
    
    If {floatize} is not {NULL}, that procedure is used to convert the sample
    values of {A} into floatng-point numbers before computng the local average.
    If it is {NULL}, the sample value is converted with {sample_conv_floatize}
    with {maxval = mA}, {isMask = TRUE}, {lo = 0.0}, {hi = 1.0}.
    
    If {quantize} is not {NULL}, it is used to convert the local average from 
    double-precision to a sample value in {0..mG} where {mG = 2^bps-1}. If {quantize}
    is {NULL}, the samples are converted with {sample_conv_quantize}
    with paramters {maxval = mG}, {isMask = TRUE}, {lo = 0.0}, {hi = 1.0}. */

double ppv_array_blur_single
  ( ppv_array_t *A, 
    const ppv_index_t ix[],
    ppv_size_t radius,
    ppv_sample_floatize_proc_t *floatize,
    double wt[]
  );
  /* Computes the weighted averages of the samples of {A}
    in a window of size {t = 2*radius+1} along each axis,
    centered at the voxel of {A} with indices {ix[0..d-1]}. 
    
    If there are no voxels of {A} in that neighborhood,
    returns {NAN}.
    
    If {floatize} is not {NULL}, it is used toconvert the samples 
    of {A} to double-precision floats. Otherwise they are converted
    with {sample_conv_floatize}.
    
    The array {wt} should have {t} elements. The weight of a 
    sample {A[jx]} in the neighborhood will be the product
    of {wt[jx[k]-ix[k]+radius]} for all axes {k}. */

#endif

