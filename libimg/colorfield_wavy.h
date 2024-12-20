/* colorfield_wavy.h - definitions for wavy fields 
** Last edited on 2004-11-06 03:13:33 by stolfi
**
** Copyright (C) 2003 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
*/

#ifndef colorfield_wavy_H
#define colorfield_wavy_H

#include <frgb.h>
#include <frgb_ops.h>
#include <colorfield.h>

/* USER ARGUMENTS FROM THE COMMAND LINE */

typedef struct cfld_wave_args_t  /* User specs for a single wave. */
  { cfld_int_pair_t bot;           /* Wave-bottom point closest to origin. */
    frgb_t botColor;       /* Raw pixel value at wave bottom. */
    struct cfld_wave_args_t *next;  /* Next wave in list, or {NULL} */
  } cfld_wave_args_t;
  /* Each {cfld_wave_args_t} record describes a wave-like adjustment 
    to the basic color value.

    An ordinary wave has, at pixel {p} and channel {c}, the value
      { ampl[c] * wfunc(p) }
    whereas a logarithmic wave has value
      { exp(ampl[c] * wfunc(p)) }
    In either case, the wave is defined in terms of
      { wfunc(p) = (1 - cos(2*Pi*phase(p)))/2 }
    where
      { phase(p) = (p.col - org[0])*fr[0] + (p.row - org[1])*fr[1] }
    Here {org} is the wave's origin and {fr = 2(bot - org)/L2} is the
    wave's frequency vector, where {L2} is the norm squared of the
    period vector {2(bot - org)}. In either case, the amplitudes
    {ampl[c]} are computed so as to produce the final colors specified
    by the user.

    Note that the function {wfunc} is null at the "origin" and 
    varies symmetrically away from that point in the direction
    of the {period} vector.  Note also that the {row} index is
    0 at the *top* row of the image, and increases *downwards*. */

typedef cfld_wave_args_t *cfld_wave_args_list_t;  /* Pointer to a list of {cfld_wave_args_t} records. */

typedef struct cfld_wavy_args_t  /* User specs for a wavy field */
  { cfld_int_pair_t org;         /* Origin pixel for wave system. */
    frgb_t orgColor;    /* Base color (where all waves are at maximum). */
    cfld_wave_args_list_t waves;  /* Wave-like corrections for {orgColor}. */
  } cfld_wavy_args_t;

cfld_wavy_args_t *cfld_wavy_parse_uniform(argparser_t *pp);
  /* Parses the next few command line arguments, starting at
    {argv[*argn]}, as a color specification {COLOR}. Assumes that the
    "-field uniform" arguments have already been parsed. 
    Increments {*argn}.  Returns the resultpackaged as a
    {cfld_wavy_args_t} with no waves. */ 

cfld_wavy_args_t *cfld_wavy_parse_simple(argparser_t *pp, int32_t uX, int32_t uY);
  /* Parses the next few command line arguments, starting at
    {argv[*argn]}, as an horizontal or vertical color-wave
    specification. Assumes that the "-field hWave" or 
    "-field vWave" arguments have been scanned. Increments {*argn}. 
    
    The arguments should have the form {POS0 COLOR0 POS1 COLOR1},
    where {POS0},{POS1} are coordinates of the wave's maximum and of
    the nearest minimum, and {COLOR0},{COLOR1} are the respective
    color values (in the {frgb_parse_color} format).
    
    The direction of the wave is defined by the unit vector {(uX,uY)},
    namely {(1,0)} for horizontal, {(0,1)} for vertical (downwards).
    
    The returned record {wa} has the wave's origin pixel in {wa->org} 
    the maximum value {COLOR0} in the {wa->color},
    and a single {cfld_wave_args_t} in {wa->waves}, containing
    the minimum value {COLOR1} and the wave's period. */ 
    
cfld_wavy_args_t *cfld_wavy_parse_general(argparser_t *pp, int32_t n);
  /* Parses the next few command line arguments, starting at
    {argv[*argn]}, as the specifications of a variable color field.
    Assumes that the "-field wave" or "-field wavePair" arguments have
    already been parsed. Increments {*argn}.
    
    The field consists of {n} waves (either 1 or 2) in
    general directions, multiplying a reference color value. The
    arguments should have the form 
      { X0 Y0 COLOR0 X1 Y1 COLOR1 ... Xn Yn COLORn }
    where {X0 Y0} is a pixel where all waves are maximum, and {Xk Yk},
    for {k>0}, is a pixel where wave {k} is minimum and the other wave
    (if any) is still maximum. The triplets {COLOR0}, {COLOR1}, ...
    {COLORn} are the color values of those pixels.
    
    The returned record {wa} has the wave's origin pixel in {wa->org} 
    the maximum value {COLOR0} in the {wa->color},
    and either one or two {cfld_wave_args_t} in {wa->waves}, containing
    the minimum values {COLOR1,COLOR2} and the respective periods. */

/* PREPROCESSED DATA */

typedef struct cfld_wave_params_t
  { cfld_int_pair_t freqNum;         /* Numerator of frequency vector {fr}. */
    int32_t freqDen;             /* Denominator of frequency vector {fr}. */
    frgb_t *tb;             /* Wave adjutment for each {iphase}. */
    struct cfld_wave_params_t *next;  /* Next wave table. */
  } cfld_wave_params_t;
  /* A {cfld_wave_params_t} is a precomputed table {wp} of wave-like adjustments
    implied by some {cfld_wave_args_t} {wa} (e.g. black-level or white-level),
    indexed by the {iphase} parameter, which is defined as 
      { iphase(p) = ((p.col - org[0])*frN[0] + (p.row - org[1])*frN[1]) MOD frD }
    where {frN = freqNum} is the wave's {period} vector reduced to its lowest
    terms, and {frD = fewDen} is the dot product of {period} and {frN}.  Note that 
    {iphase} ranges in {[0..frD-1]}. */

cfld_wave_params_t *cfld_wavy_compute_wave_params
  ( cfld_wave_args_t *wa,
    cfld_int_pair_t *org,
    frgb_t *orgColor, 
    frgb_t *botColor,
    bool_t logarithmic
  );
  /* Precomputes a wave adjustment table {wp} for the wave specs {*wa},
    given the properly corrected values {orgColor} and {botColor}
    at the wave's origin {org} and bottom points {wa->bot}, respectively. */
    
typedef cfld_wave_params_t *cfld_wave_params_list_t; /* A list of {cfld_wave_params_t} records. */

typedef struct cfld_wavy_params_t  /* Precomputed data for a wavy field */
  { cfld_int_pair_t org;         /* Origin pixel for wave system. */
    frgb_t orgColor;    /* Base color (where all waves are at maximum). */
    cfld_wave_params_list_t waves; /* Wave-like corrections for {orgColor}. */
  } cfld_wavy_params_t;

cfld_wavy_params_t *cfld_wavy_compute_params
  ( cfld_wavy_args_t *wfa, 
    frgb_adjuster_t *adjust,
    bool_t logarithmic
  );
  /* Computes the reference color {wfp->orgColor} 
    and its wave adjustment tables {wfp->tables}. */

/* EVALUATION */

void cfld_wavy_eval
  ( cfld_wavy_params_t *wfp,
    bool_t logarithmic, 
    int32_t col, 
    int32_t row,
    frgb_t *fv,
    int32_t chns
  );
  /* Evaluates the wavy field described by {wfp} at the pixel 
    in column {col} and row {row}. */

void cfld_wavy_apply_waves
  ( frgb_t *locColor,
    int32_t chns,
    int32_t dCol, 
    int32_t dRow,
    cfld_wave_params_list_t wp,
    bool_t logarithmic
  );
  /* Applies the wave corrections {wp} to the color
    {locColor[0..chns-1]}, evaluated at the pixel in column {dCol} and
    row {dRow} relative to the wave's origin. Assumes that the wave
    corrections for all values of {iphase} have been precomputed in
    {wp}. If {logarithmic} is true, the corrections are multiplied,
    otherwise they are added. */

#endif
