#ifndef frgb_path_H
#define frgb_path_H

/* frgb_path.h - mapping real numbers to paths in the RGB cube
** Last edited on 2013-05-24 22:43:39 by stolfilocal
**
** Copyright (C) 2007 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
*/

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
  
  Some procedures take a {cycles} parameters that determines
  how varied are the hues (and, sometimes, saturation). */ 

frgb_t frgb_path_map_unsigned_0(double z);
  /* Chooses a pseudocolor appropriate for a pixel with value {z}.
    Maps 0 (or less) to black, 1 (or greater) to white, and
    intermediate values to intermediate colors. */ 

frgb_t frgb_path_map_unsigned_1(double z, double H0, double Y0, double H1, double Y1);
  /* A color path consisting of a polygonal spiral of saturated colors
    that sweeps from hue {H0} and brightness {Y0} to hue {H1} and brightness {Y1}.

    The point {R,G,B} corresponding to a given {z} is the maximally
    saturated color with brightness {Y = (1-z)*Y0 + z*Y1} and hue {H =
    (1-z)*H0 + z*H1}, where {H0} and {H1} are two client-given hues.
    Note that by giving {H1 = H0 + k} one obtains a path that cycles
    {k} times over the whole hue spectrum. */

frgb_t frgb_path_map_signed_0(double z, int cycles);
frgb_t frgb_path_map_signed_1(double z, int cycles);
frgb_t frgb_path_map_signed_2(double z, int cycles);
  /* Chooses a pseudocolor appropriate for a pixel with value {z},
    assumed to range in [-1_+1]. Value {0} is mapped to 50% gray,
    values in [0_+1] are mapped to light warm colors, and values in
    [-1_0] to the complementary (bluish) colors. In each half-range, 
    the hues will sweep through the spectrum {cycles} times. */ 
    
frgb_t frgb_path_map_signed(double z, int cycles, int style);
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
