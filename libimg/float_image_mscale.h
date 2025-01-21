#ifndef float_image_mscale_H
#define float_image_mscale_H

/* Tools for multiscale image processing. */
/* Last edited on 2025-01-14 21:03:44 by stolfi */

#include <bool.h>
#include <r2.h>
#include <float_image.h>

/* 
  HALF-SIZE REDUCTION
  
    !!! The old behavior for slope maps is {NXR=(NXA+1)/2} and {NYR=(NYA+1)/2}. !!!
    
    !!! Old behavior is obtained with {dx,dy == (k-1),(k-1)} if {nw=2*k}, and {k,k} !!!
    !!! if {nw=2*k+1}; that is, {dx,dy == (nw-1)/2,(nw-1)/2} rounded down in any case. !!!

    !!! Clients may want to make {R} concentric with {A}, for even or odd size. !!!
    !!! This requires the parity of {NX} and {nwx} be the same, ditto for {NY},{nwy}. !!!
    
    The procedures below are concerned with creating a reduced image
    {R} from an image {A}. The image {R} will have {NXR} columns and {NYR} rows.
    The number of channels is the same.
    
    If {A} has {NXA} columns and {NYA} rows, typically one chooses
    {NXR=NXA/2} and {NYR=NYA/2}, rounded up or down.

    Each pixel {R[x,y]} of {R} is a weighted average of pixels of {A},
    within a window with {nw × nw} pixels, namely
    {A[2*x+r-dx,2*y+s-dy]} for {r,s} in {0..nw-1}. The window weights
    are binomial coefficients {comb(nw,i)*comb(nw,j)}, and have the
    partition of unit property iff {nw >= 2}.

    If {nw = 1}, then {R[x,y]} is merely {A[2*x-dx,2*y-dy]}, ie. the
    {R} array is obtained by subsampling and not averaging, starting
    with sample {A[dx,dy]}. In this case 3/4 of the {A} samples are
    simply ignored and the partition-of-unit property does not hold.
    
    Irrespective of whether one views samples as being located at the
    corners of the integer grid or at the centers of pixels, point
    {(x,y)} in the domain of {R} corresponds to point
    {(2*x+(nw-1)/2-dx,2*y+(nw-1)/2-dy)} of the domain of {A}. */

float_image_t *float_image_mscale_shrink(float_image_t *A, float_image_t *M, int32_t NXR, int32_t NYR, int32_t dx, int32_t dy, int32_t nw);
  /* Returns a version {R} of {A} downscaled by a factor of {1/2}. 
  
    If {M} is not null, it must be a single-channel image with the
    same size as {A}; in that case, each sample of {M} is assumed to
    be the reliability of the corresponding pixel of {A}, and is
    multiplied into the window weight. Samples of {A} which are NAN or
    lie outside {A}'s domain are assumed to have zero reliability
    weight.

    In any case, if all contributing samples for a
    sample of {R} have zero weight, the sample is set to NAN. */

float_image_t *float_image_mscale_mask_shrink(float_image_t *M, int32_t NXR, int32_t NYR, int32_t dx, int32_t dy, int32_t nw, bool_t harm);
  /* Similar to {float_image_mscale_shrink}, but specialized for weight masks.
   
    The average is the weighted arithmetic mean if {harm} is false, or
    the weighted harmonic mean if {harm} is true. */

r2_t float_image_mscale_point_shrink(r2_t *pA, int32_t dx, int32_t dy, int32_t nw);
  /* If {R} is the result of {float_image_mscale_shrink(A,M,dx,dy,nw)}, 
    returns the point of {R}  domain that corresponds 
    to point {p} in {A}'s domain. */

r2_t float_image_mscale_point_expand(r2_t *pR, int32_t dx, int32_t dy, int32_t nw);
  /* If {R} is the result of {float_image_mscale_shrink(A,M,dx,dy,nw)}, 
    returns the point of {A} domain that corresponds 
    to point {p} in {R}'s domain. */

/* !!! Add {float_image_mscale_expand} generalizing {pst_height_map_expand} !!! */

int32_t float_image_mscale_rounding_bias(int32_t n);
  /* A bias bit that may be useful to keep the image centered during reduction.
      For n =  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 
          b =  0  1  0  0  0  1  0  1  0  1  0  0  0  1  0  0  0  1
     !!! Find a use for this procedure !!! 
  */

/* DEBUGGING I/O */

char *float_image_mscale_file_name(char *filePrefix, int32_t level, int32_t iter, char *tag, char *ext);
  /* Returns a file name from the given elements, namely 
    "{filePrefix}-{LL}-{NNNNNNNNN}-{tag}.fni"
    where "{LL}" is the level (2 digits) and "{NNNNNNNNN}" is {iter} 
    (9 digits). If {level} is negative, omits the "-{level}" part. Ditto
    for the "-{iter}" part.  If {filePrefix} is "-", ignores all other 
    arguments and returns a new string containing just "-". */

void float_image_mscale_write_file(float_image_t *M, char *filePrefix, int32_t level, int32_t iter, char *tag);
  /* Writes the float image {M} in a human-readable FNI format.
    The file name is created by {float_image_mscale_file_name} with extension "fni".
    Diagnostic messages are indented by {2*level+2} if {level} is non-negative. */

#endif
