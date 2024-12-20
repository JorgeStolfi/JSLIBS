/* Choosing sampling points and rays. */
/* Last edited on 2024-12-14 08:43:07 by stolfi */

#ifndef multifok_sampling_H
#define multifok_sampling_H

#include <stdint.h>

#include <interval.h>
#include <r3.h>
#include <hr3.h>
#include <bool.h>
#include <frgb.h>
#include <float_image.h>

#include <multifok_image.h>
#include <multifok_frame.h>
  
typedef struct multifok_sampling_t
  { uint32_t NS;   /* Number os sampoints (sub-sampling points) per pixel. */
    i2_t *iSmp;    /* Indices of sampoints to use; {(0,0)} is pixel center. */
    double step;   /* Distance between sampoints. */
    double *wSmp;  /* Weights of sampoints in pixel averages. */
    uint32_t KR;   /* Numbe of rays to throw per sampoint. */
  } multifok_sampling_t;
  /* Features of the simulated camera lens.
    
    The arrays {iSmp[0..NS-1]} and {wSmp[0..NS-1]} and the parameter
    {step} define {NS} sampoints (sub-sampling points) for each pixel.
    Sampoint {k} is {ctr + iSmp[k]*step} where {ctr} is the the pixel's
    center, all in image coordinates; and its weight is {wSmp[k]}.
    
    The parameter {KR} specifies the number of rays to cast for each sampoint. */

multifok_sampling_t *multifok_sampling_choose(uint32_t HS, uint32_t KR, bool_t verbose);
  /* Returns a {multifok_sampling_t} record {samp} with sampling parameters 
    for ray tracing.   Allocates the arrays and returns the choices in
    {*NS_P}, {*iSmp_P}, {*wSmp_P}. The weights {wSmp[0..NS-1]} will be all positive.
    
    The non-negative parameter {HS} determines the maximum number of
    sampoints. The number {NS} of sampoints will be at most
    {(2*HS+1)^2}, and their indices {iSmp[0..NS-1]} will be a subset of
    the square grid {{-HS..+HS}×{-HS..+HS}}. For any {HS}, the first
    sample {iSmp[0]} will be {(0,0)} and {wSmp[0]} will be the largest
    weight. Any other rays will have {iSmp} distinct from {(0,0)}. Both
    arrays will be sorted so that {|iSmp[k]|} increases with {k}.
    
    IN particular, if {HS} is zero, {NS} will be 1, so there will be just the {iSmp[0]=(0,0)}
    point with weight {wSmp[0]=1.0}. */

void multifok_sampling_free(multifok_sampling_t *samp);
  /* Reclaims the storage used by {samp}, including
    internal tables. */

#endif
