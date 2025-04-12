/* Tools for weight maps. */
/* Last edited on 2025-04-10 22:27:38 by stolfi */

#ifndef multifok_weight_H
#define multifok_weight_H

#include <stdint.h>

#include <float_image.h>

float_image_t *multifok_weight_from_height_and_normal(float_image_t *htd, double dhMax, float_image_t *nrm, double dnMax);
  /* Computes a reliability weight map for the height and normal of the
    visible scene's surface, given the height deviation map {htd} and
    the normal map {nrm}.
  
    The images {htd} and {nrm}, if not null, shoud have the same col
    and row counts. The image {htd} should have at least 1 channel,
    with non-negative values, and channel 0 is interpreted as the
    standard deviation of the height of the scene's surface visible in
    each pixel. The image {nrm} should have at least three channels,
    and channels {0..2} are interpreted as the {X}, {Y}, and {Z}
    coordinates of the average unit vector that is normal to the scene's
    surface visible in each pixel.
    
    The weight {wht[0,x,y]} of a pixel is zero if (1) the value of {htd}
    at the pixel exceeds {dhMax}, or (2) the normal vector at the pixel
    is a null vector, or (3) the modulus of the difference between the
    {nrm} vector at a pixel and at any of its 8 neighbors exceeds
    {dnMax}, or (4) any of the samples in {htd[0,x,y]} and {nrm[0..2,x,y]}
    is not finite. Otherwise {wht[0,x,y]} is 1.
    
    However, if {htd} has two or more channels, channel 1 is
    interpreted as an 'a priori' reliability weight for the data in
    channel 0. Likewise, if {nrm} has four or more channels, channel 3
    is interpreted as a reliability weight for channels {0..2}. These
    prior weights, if they exist, should be finite non-negative numbers.
    If any of them is zero, the weight {wht[0,x,y]} is set to zero too.
    
    If {htd} is null, the tests involving it are skipped, and {dhMax} is ignored.
    Likewise if {nrm} is null, its tests are skipped and {dnMax} is ignored. */
    
#endif    
