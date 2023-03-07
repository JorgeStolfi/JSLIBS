/* frgb_interp_vis.h - interpolation of colors. */
/* Last edited on 2023-03-07 13:57:48 by stolfi */

#ifndef frgb_interp_vis_H
#define frgb_interp_vis_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>

#include <frgb.h>

frgb_t frgb_interp_vis(double r, frgb_t *p0, frgb_t *p1, bool_t logY);
  /* Interpolates a point {r} of the way between RGB triples {p0} and {p1}
    in a visual-sensitive way.  The parameter {r} must be in {[0_1]}.
    
    Basically converts {p0,p1} to 
    HTY coordinate triples {q0,q1}, calls {frgb_interp_vis_HTY(r,q0,q1,logY)},
    and converts the result back to RGB coordinates. */
    
frgb_t frgb_interp_vis_HTY(double r, frgb_t *q0, frgb_t *q1, bool_t logY);
  /* Interpolates a point {r} of the way between colors {q0} and {q1} in
    a visual-sensitive way.  The parameter {r} must be in {[0_1]}.
    
    The arguments {q0,q1} and the result {q} are assumed to be in the
    UV-based HTY coordinate system. See {frgb_to_HTY_UV} for its
    description.
    
    The luminosity {Y} of the returned color will be interpolated
    between the luminosities {Y0,Y1} of {q0} and {q1}, which should be
    non-negative. If {logY} is false, the interpolation will be affine
    ("linear"). If {logY} is true the interpolation will be in log
    scale, with a small bias.
    
    The hue {H} will be interpolated affinely ("linearly") between the
    hues {H0,H1} of {q0,q1}. These hues need not be in the range
    {[0_1)}. If their difference is {dH}, the hue {H} will make {dH}
    turns around the hue circle as {r} goes from 0 to 1. In any case,
    the resulting {H} will be reduced modulo 1 to the range {[0_1)}.
    
    The relative saturation {T} will be interpolated exponentially
    between the saturations {T0,T1} of {q0,q1}. These coordinates should
    be in {[0_1]}. */

#endif

