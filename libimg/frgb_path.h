#ifndef frgb_path_H
#define frgb_path_H

/* frgb_path.h - mapping real numbers to paths in the RGB cube */
/* Last edited on 2023-03-07 16:22:57 by stolfi */
/* See the rights and conditions notice at the end of this file. */

#define _GNU_SOURCE
#include <stdint.h>

#include <frgb.h>
#include <bool.h>

/* COLOR SCALES */

/* The procedures in this section generate pseudocolor
  scales.  They map a real parameter {z} to an {R,G,B} 
  triple by a more or less convoluted path, usually continuous.
  
  The {unsigned} versions assume that the parameter {z} varies over
  the unsigned unit interval {U = [0_1]}. They usually map {0} to
  black and 1 to white. The brightness {Y} of the result usually
  increases monotonically with {z}.
  
  The {signed} versions assume that {z} varies over the signed unit
  interval {V = [-1_+1]}. They usually map {0} to a neutral color, and
  map {+z} and {-z} to colors with complementary hues. For some of
  them, {-1} is mapped to black, {+1} to white, and {Y} increases
  monotonically with {z} over the whole range {V}. For others, {Y} is
  a monotonic function of {|z|}.
  
  Some procedures take a {cycles} parameters that determines how varied
  are the hues (and, sometimes, saturation). The simplest color
  gradation is obtained with {cycles=1}. */ 

frgb_t frgb_path_map_unsigned_0(double z, int32_t cycles);
  /* Chooses a pseudocolor appropriate for a pixel with value {z}.
    Maps 0 (or less) to black, 1 (or greater) to white, and
    intermediate values to intermediate colors.
    
    The returned color will be maximally saturated, that is, on the
    surface of the unit RGB cube. The {R}, {G}, and {B} coordinates are
    continuous functions of {z}. The hue and brightness are also smooth
    (C1) functions of {z}, but the saturation has kinks (derivative
    jumps).
    
    The hue will make {|cycles|-1} full sweeps over the hue spectrum,
    plus one partial sweep over the arc from purple to yellow. If
    {cycles} is negative, the hues will run in the opposite sense.  If
    {cycles} is zero, the result is a gray scale. */ 

int32_t frgb_path_map_unsigned_max_style(void);
  /* The largest value of the {style} parameter for an unsigned path. */

frgb_t frgb_path_map_unsigned(double z, int32_t cycles, int32_t style);
  /* Same as {frgb_path_map_unsigned_{style}}. */ 
    
frgb_t frgb_path_map_signed_0(double z, int32_t cycles);
frgb_t frgb_path_map_signed_1(double z, int32_t cycles);
frgb_t frgb_path_map_signed_2(double z, int32_t cycles);
  /* Chooses a pseudocolor appropriate for a pixel with value {z},
    assumed to range in [-1_+1]. Value {0} is mapped to a gray tone,
    values in [0_+1] are mapped to light warm hues, and values in
    [-1_0] to cold hues.
    
    In each half-range, the hue will make {|cycles|-1} times sweeps over
    the appropriate half of the hue spectrum, plus one partial sweep. If
    {cycles} is negative, the hues will run in the opposite sense. If
    {cycles} is zero, the result is a gray scale. */  

int32_t frgb_path_map_signed_max_style(void);
  /* The largest value of the {style} parameter for a signed path. */

frgb_t frgb_path_map_signed(double z, int32_t cycles, int32_t style);
  /* Same as {frgb_path_map_signed_{style}}. */ 

/* PATH PARAMETER ADJUSTMENT */

double frgb_path_reduce_unsigned(double z);
  /* Reduces the unsigned path parameter {z} modulo 1 to a number in
    {[0_1]}. That is, returns the fractional part {f} of {z}, as a
    non-negative number in {[0_1]}.
    
    More precisely, if {z} is positive, returns {f} in {(0_1]},
    otherwise returns {f} in {[0,1)}; in any case, {z-f} is
    integer. In particular, if {z} is {0}, returns {f=0}. */

double frgb_path_reduce_signed(double z);
  /* Reduces the signed path parameter {z} modulo 1 to a number in
    {[-1_+1]}. That is, returns the fractional part {f} of {z},
    preserving its sign.
    
    More precisely, if {z} is {0}, returns {f=0}. If {z} is positive,
    returns {f} in {(0_+1]); if {z} is negative, returns {f} in
    {[-1_0)}. In any case, {z-f} will be an integer. */

#endif

