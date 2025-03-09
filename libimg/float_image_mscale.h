#ifndef float_image_mscale_H
#define float_image_mscale_H

/* Tools for multiscale image processing. */
/* Last edited on 2025-03-02 01:31:44 by stolfi */

#include <bool.h>
#include <r2.h>
#include <float_image.h>

float_image_t *float_image_mscale_shrink
  ( float_image_t *A,
    int32_t wch,
    bool_t harm,
    int32_t NXR,
    int32_t NYR,
    int32_t dx,
    int32_t dy
  );
  /* Returns a version {R} of {A} reduced by a factor of 1/2.
    
    More precisely, the image {R} will have {NXR} columns and {NYR}
    rows. The number of channels will be the same. If {A} has {NXA}
    columns and {NYA} rows, typically one chooses {NXR=NXA/2} and
    {NYR=NYA/2}, rounded up or down.
    
    Each pixel {R[X,Y]} of {R} is an average of a {2×2} block of pixels of {A},
    namely {A[X',Y']} where {X'=2*X-dx+r} and {Y'=2*Y-dy+s} for {r,s} in {0..1}. 
    These blocks will be pairwise disjoint.
    
    Irrespective of whether one views samples of {A} and {R} as being
    located at the corners of the integer grid or at the centers of
    grid cells, point {(x,y)} in the domain of {R} corresponds to point
    {(2*x-dx+0.5,2*y-dy+0.5)} of the domain of {A}.
  
    If {wch} is the index of a channel of {A}, each sample {A[wch,X,Y]}
    is taken to be the reliability of the samples {A[c,X,Y]} for all
    {c!=wch}. These reliability weights must be finite and positive. In
    this case each sample {R[c,X,Y]} is a weighted average of the
    samples {A[c,X',Y'] in the corresponding {2×2} block, with these
    weights. In this case, channel {wch} of the {R} will be obtained
    from channel {wch} of {A} by UNWEIGHTED arithmetic mean if {harm} is
    false, or the weighted harmonic mean if {harm} is true.
    
    If {wch} is not the index of a channel of {A}, the reliability 
    weight of every pixel is assumed to be 1.0, and {R} will have no
    reliability weight channel.

    Samples of {A} which are NAN or lie outside {A}'s domain are assumed
    to have zero reliability weight. In any case, if all contributing
    samples for a sample of {R{c,X,Y]} have zero weight, the sample is
    set to NAN, and the reliability weight {R[wch,X,Y]} (if it exists)
    is set to zero. */
    
float_image_t *float_image_mscale_expand
  ( float_image_t *R,
    int32_t NXA,
    int32_t NYA,
    int32_t dx,
    int32_t dy
  );
  /* The approximate inverse of {pst_float_image_mascale_shrink}:
    returns an image {A} that is the image {R} expanded by a factor of 2.
    
    More precisely, the image {A} will have {NXA} columns and {NYA}
    rows. The number of channels will be the same. If {R} has {NXR}
    columns and {NYR} rows, typically one should have {NXR=NXA/2} and
    {NYR=NYA/2}, rounded up or down.
    
    Each pixel {R[c,X,Y]} is replicated into a {2×2} block of samples of
    {A}, namely {A[c,2*X-dx+r,2*Y-dy+s]} for {r,s} in {0,1}. These blocks
    will be pairwise disjoint. */
  
r2_t float_image_mscale_point_shrink(r2_t *pA, int32_t dx, int32_t dy);
  /* If {R} is the result of {float_image_mscale_shrink(A,wch,harm,,NXR,NYR,dx,dy)}
    for some {wch} and {harm}, returns the point of {R}  domain that corresponds 
    to point {p} in {A}'s domain. */

r2_t float_image_mscale_point_expand(r2_t *pR, int32_t dx, int32_t dy);
  /* If {R} is the result of {float_image_mscale_shrink(A,M,dx,dy)}, 
    returns the point of {A} domain that corresponds 
    to point {p} in {R}'s domain. */

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
    (9 digits). 
    
    If {level} is negative, omits the "-{LL}" part. The {NNNNNNNNN} part
    isreplaced by "beg" if {iter} is zero, and "end" if {iter} is
    negative.
    
    If {filePrefix} is "-", ignores all other arguments and returns a
    new string containing just "-". */

void float_image_mscale_write_file(float_image_t *M, char *filePrefix, int32_t level, int32_t iter, char *tag);
  /* Writes the float image {M} in a human-readable FNI format.
    The file name is created by {float_image_mscale_file_name} with extension "fni".
    Diagnostic messages are indented by {2*level+2} if {level} is non-negative. */

#endif
