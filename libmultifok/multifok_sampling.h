/* Choosing sampling points and rays. */
/* Last edited on 2024-10-22 13:05:13 by stolfi */

#ifndef multifok_sampling_H
#define multifok_sampling_H

#define _GNU_SOURCE
#include <stdint.h>

#include <interval.h>
#include <r3.h>
#include <hr3.h>
#include <bool.h>
#include <frgb.h>
#include <float_image.h>

#include <multifok_image.h>
#include <multifok_frame.h>

void multifok_sampling_choose_pixel_sampling_points_and_weights
  ( int32_t HS, 
    int32_t *NS_P,
    r2_t **uSmp_P, 
    double **wSmp_P,
    bool_t verbose
  );
  /* Chooses a number {NS} pixel sub-sampling points, their deviations
    {uSmp[0..NS-1} from the pixel's center, and their weights
    {wSmp[0..NS-1]}. Allocates the arrays and returns the choices in
    {*NS_P}, {*uSmp_P}, {*wSmp_P}.

    The weights {wSmp[0..NS-1]} will be all positive.
    
    For any {HS}, the first sample {uSmp[0]} will be {(0,0)} and
    {wSmp[0]} will be the largest weight. Any other rays will have
    {uSmp} distinct from {(0,0)}. Both arrays will be sorted so that
    {|uSmp[ks]|} increases with {ks}.
    
    If {HS} is zero, {NS} will be 1, so there will be just the {(0,0)}
    point. Otherwise the number of sampling points will be greater than 1,
    about {(2*HS+1)^2}. */

void multifok_sampling_choose_ray_tilts_and_weights
  ( int32_t KR, 
    int32_t NS,
    int32_t *NR_P,
    r2_t **tRay_P, 
    double **wRay_P,
    bool_t verbose
  );
  /* Chooses a number {NR} of rays, their relative deviations from the
    vertical {tRay[0..NR-1}, and their weights {wRay[0..NR-1]}.
    Allocates the arrays and returns the choices in {*NR_P}, {*tRay_P},
    {*wRay_P}.   
    
    If {KR=0}, the procedure will generate a single ray direction ({NR =
    1}) irrespective of {NS}. If {KR>=1}, then {NS} should be the number
    of subsampling points per pixel, and must be odd. The number of rays
    {NR} will be {KR*NS} rounded up to the smallest odd square that is a
    multiple of {NS}.

    The weights {wRay[0..NR-1]} will be positive.
    
    For any {NR}, the first ray will be strictly vertical, with {tRay[0]
    = (0,0)} and {wRay[0]} maximum. Any other rays will have {tRay}
    distinct from {(0,0)}. Both arrays will be sorted so that
    {|tRay[ks]|} increases with {ks}.
    
    If {NR} is 1, there will be just that vertical ray. Otherwise the
    RMS radius of the points {tRay}, (namely the square root of the
    weighted average of {|tRay[ir]|^2} with the {wRay} weights,
    including {tRay[0]) will be 1.0. */

#endif
