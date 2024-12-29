#ifndef pst_basic_H
#define pst_basic_H

/* pst_basic.h -- basic data types for gauge-based photostereo. */
/* Last edited on 2024-12-29 01:53:32 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <values.h>

#include <float_image.h>
#include <vec.h>
#include <r2.h>
#include <r3.h>
#include <frgb.h>
#include <argparser.h>

#define INF INFINITY
  /* Plus infinity. */

vec_typedef(name_vec_t,name_vec,char *);
  /* Defines the type {name_vec_t} as a vector of {char*} elems. */

vec_typedef(image_vec_t,image_vec,float_image_t *);
  /* Defines the type {image_vec_t} as a vector of {float_image_t*} elems. */

typedef r3_t pst_normal_func_t (r2_t *p);
  /* A procedure that computes the normal direction {nrm} at a visible point
    {P} of a some surface, given the projection {p} of that point in some
    plane. 
    
    Both {p} and {nrm} are given in some orthogonal {U,V,W} coordinate
    system such that the projection of point {(u,v,w)} has coordinates
    {(u,v)} (i.e., such that the {W} axis is parallel to the direction
    of projection). The returned normal should have a non-negative {W}
    component.
    
    The procedure should return {(NAN,NAN,NAN)} if 
    the normal direction is not defined at the point
    {P} (e.g. if {P} is at infinity).  Otherwise it should 
    return a valid unit-length vector. */ 
 
typedef frgb_t pst_albedo_func_t (r2_t *p);
  /* A procedure that computes the albedo (intrinsic color) {alb} at a
    visible point {P} of a some surface, given the projection {p} of
    that point in some plane. The components of {alb} must be finite and
    usually (but not necessarily) between 0 and 1.  
    
    The procedure should return the invalid color
    {(NAN,NAN,NAN)} if the albedo is not defined at {p}.
    Otherwise it should return an {frgb_t} color 
    with finite components, usually (but not necessarily)
    between 0 and 1. */ 

/* TUPLES OF VALUES FOR CHANNELS */

void pst_double_vec_regularize(double_vec_t *v, int32_t NC, double defval);
  /* Make sure that the vector {v} has {NC} elements. If
    {v.ne == 0}, expands it to {NC} elements, and sets all elements
    to {defval}. If {v.ne == 1}, expands it to {NC} elements,
    replicating the first element. If {v.ne == NC},
    does nothing. Otherwise fails. */

void pst_double_vec_uniformize(double_vec_t *v, double defval);
  /* Make sure that the vector {v} has one element. If
    {v.ne == 0}, expands it to one element and sets it to {defval}.
    If {v.ne > 1}, requires that all elements be equal, and
    truncates it to one element. Otherwise does nothing. */

double_vec_t pst_double_vec_parse(argparser_t *pp, int32_t *NC);
  /* Parses a tuple of values from the command line, with optional 
    denominator, in the format described by {pst_double_vec_spec_HELP} and 
    {pst_double_vec_spec_INFO}.  See {argparser.h} for an explanation 
    of the {pp} parameter.
    
    If {NC} is is NULL, or only one numeric argument is present (with
    optional denominator), ignores {NC}. Otherwise, if {*NC} is
    negative, sets {*NC} to the number of elements read. Otherwise
    demands and parses exactly {*NC} numeric arguments (with an
    optional denominator). */
  
#define pst_double_vec_spec_HELP \
  "{NUM} .. " pst_double_vec_spec_den_HELP
  
#define pst_double_vec_spec_den_HELP \
  "[ / {DEN} ]"
  
#define pst_double_vec_spec_INFO \
  "The argument consists of one {NUM} value" \
  " for each channel, or by a single {NUM} that" \
  " applies to all channels.  " \
  pst_double_vec_spec_den_INFO
  
#define pst_double_vec_spec_den_INFO \
  "If the \"/ {DEN}\" part is present," \
  " the given values are divided by {DEN}."

/* TUPLES OF INTEGERS FOR CHANNELS */

void pst_int_vec_regularize(int32_vec_t *v, int32_t NC, int32_t defval);
  /* Make sure that the color vector {v} has {NC} color channels. If
    {v.ne == 0}, expands it to {NC} channels, and sets all elements
    to {defval}. If {v.ne == 1}, expands it to {NC} channels,
    replicating the first element. Otherwise does nothing. */

int32_vec_t pst_int_vec_parse(argparser_t *pp, int32_t *NC);
  /* Parses a tuple of values from the command line, in the format
    described by {pst_int_vec_spec_HELP} and {pst_int_vec_spec_INFO}.
    See {argparser.h} for an explanation of the {pp} parameter.
    
    If {NC} is is NULL, or only one numeric argument is present (with
    optional denominator), ignores {NC}. Otherwise, if {*NC} is
    negative, sets {*NC} to the number of elements read. Otherwise
    demands and parses exactly {*NC} numeric arguments (with an
    optional denominator). */
  
#define pst_int_vec_spec_HELP \
  "{NUM} .. " pst_int_vec_spec_den_HELP
  
#define pst_int_vec_spec_den_HELP \
  "[ / {DEN} ]"
  
#define pst_int_vec_spec_INFO \
  "The argument consists of one {NUM} value" \
  " for each channel, or by a single {NUM} that" \
  " applies to all channels.  " \
  pst_int_vec_spec_den_INFO
  
#define pst_int_vec_spec_den_INFO \
  "If the \"/ {DEN}\" part is present," \
  " the given values are divided by {DEN}."

#endif
