/* Tools for weight maps. */
/* Last edited on 2025-04-13 13:31:03 by stolfi */

#ifndef multifok_weight_H
#define multifok_weight_H

#include <stdint.h>

#include <float_image.h>

float_image_t *multifok_weight_from_height_dev(float_image_t *hDev, double dhMax);
  /* Computes a reliability weight map {hWht} for the height map the
    visible scene's surface, given the height deviation map {hDev}.
  
    The image {hDev}, which must not be null, should have at least 1
    channel, and channel 0 is interpreted as
    the standard deviation of the height of the scene's surface visible
    in each pixel.  
    
    The deviation {hDev[0,x,y]} for each pixel {[x,y]} must be either
    {NAN} or a finite non-negative value. The resulting weight
    {hWht[0,x,y]} is zero if the {hDev[0,x,y]} at the pixel exceeds
    {dhMax} or is {NAN}. Otherwise {hWht[0,x,y]} is in principle 1.
    
    However, if {hDev} has two or more channels, channel 1 is
    interpreted as an independebt reliability weight map for the data in
    channel 0. These independent weights, if present and not {NAN}, must
    be finite and non-negative numbers, which are multiplied into the
    final weight {hWht[0,x,y]}. In particular, if {hDev[1,x,y]} is zero,
    {nWht[0,x,y]} is set to zero too. */

float_image_t *multifok_weight_from_normal_grad(float_image_t *sNrm, double dnMax);
  /* Computes a reliability weight map {nWht} for the normal map of the
    visible scene's surface, given the the normal map {sNrm}.
  
    The image {sNrm}, which must not be null should have at least 3 channels.
    Channels {0..2} are interpreted as the {X}, {Y}, and {Z}
    coordinates of the average unit vector that is normal to the scene's
    surface visible in each pixel.
    
    For each pixel {[x,y]}, the vector {sNrm[0..2,x,y]} must be either
    {(0,0,0)}, {(NAN,NAN,NAN)}, meaning that the normal is undefined; or
    a nonzero 3-vector with finite components. The resulting weight
    {nWht[0,x,y]} will be zero in the first two cases. or if the modulus
    of the difference between {sNrm[0..2,x,y]} it and {sNrm[0..2,x1,y1]}
    any pixel {[x1,y1]} in its 8 nearest neighbors exceeds {dnMax}
    (excluding neighbors that do not exist or are themselves undefined).
    Both vectors are normalized before the comparison. Otherwise
    {nWht[0,x,y]} is in principle 1.
    
    However, if {sNrm} has four or more channels, channel 3 is
    interpreted as an independent reliability weight for channels
    {0..2}. These prior weights, if present and not {NAN}, should be
    finite non-negative numbers, which are multiplied into
    {nWht[0,x,y]}. In particular, if {sNrm[3,x,y]} is zero,
    {nWht[0,x,y]} is set to zero too. */
    
void multifok_weight_multiply(float_image_t *dst, int32_t cd, float_image_t *src, int32_t cs);
  /* Multiples the weight channel of {dst}, assumed to be {cd},
    by the weight channel of {src}, assumed to be {cs}.  
    
    If channels {cd} and/or {cs} do not exist in the respective maps, does nothing.
    
    Each weight {dst[cd,x,y]} and {src[cs,x,y]} must be either {NAN} or a
    finite non-negative value.
    
    In particular, if {src[cs,x,y]} is {NAN}, {dst[cd,x,y]} is left
    unchanged).  If {src[cs,x,y]} is zero, then {dst[cd,x,y]} is set to
    zero (even if it is {NAN}). Otherwise, if {dst[cd,x,y]} is {NAN}, it
    is left unchanged. */

#endif    
