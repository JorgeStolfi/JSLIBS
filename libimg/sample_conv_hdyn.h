#ifndef sample_conv_hdyn_H
#define sample_conv_hdyn_H

/* {sample_conv_hdyn.h} - conversion between floating-point and integer samples. */
/* Last edited on 2017-06-16 01:55:17 by stolfilocal */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include <interval.h>
#include <sample_conv.h>

/* CAMERA MODEL */

#define sample_conv_hdyn_INFO \
  "  The integer pixel value {Y} is assumed to be derived from an internal " \
  " brightness value {Z} by the formula\n" \
  "\n" \
  "    {Y = round((w-k)*(max{0, min{1, c*(Z + G(s)) + b}}^g) + k)},\n" \
  "\n" \
  " where {c} is the contrast, {b} is the brightness, {s} is the noise" \
  " level, {G(s)} is a normal random variable with mean 0 and" \
  " deviation {s}, {g} is the gamma exponent, {k} is the black" \
  " offset, and {w} is the white limit."

interval_t sample_conv_hdyn_floatize
  ( sample_uint32_t iv,      /* Input integer sample ({Y} in formula). */
    double brght,          /* Brightness setting ({b}). */
    double ctrst,          /* Contrast setting ({c}). */ 
    double sigma,          /* Noise level ({s}). */
    double gamma,          /* Power law exponent ({g}). */ 
    sample_uint32_t black,   /* Black offset ({k}). */ 
    sample_uint32_t white,   /* White limit ({w}). */ 
    /* Statistical accumulators: */
    sample_uint32_t *imin,   /* Minimum {iv} seen. */
    sample_uint32_t *imax,   /* maximum {iv} seen. */
    int *clo,              /* Count of underexposed pixels. */
    int *chi,              /* Count of overexposed pixels. */
    float *vmin,           /* Minimum finite {LO(fv)} seen. */
    float *vmax            /* Maximum finite {HI(fv)} seen. */
  );
  /* Converts an integer sample {iv} to an abstract brightness
    value interval {fv} according to the camera model above.
    
    If {iv <= black}, returns the interval {[-INF,Zmax]}, where {Zmax}
    is the {Z} largest value that gives {Y == black} in the formula,
    and increment {*clo}; if {iv >= white}, returns {[Zmin,+INF]},
    where {Zmin} is the least {Z} that gives {Y == white}and increment
    {*chi}. Otherwise computes {fv = Z} by inverting the formula
    {sample_conv_hdyn_INFO} above and assuming that the noise term
    {G(s)} ranges in {[-3*sigma _ +3*sigma]}, and that the quantization
    error ranges in {[-0.5 _ +0.5]}.
    
    In any case, updates the range {imin,imax} to enclose the input
    value {iv}, and updates the range {vmin,vmax} to enclose the range
    {fv}, if finite. */
    
void sample_conv_hdyn_print_floatize_stats
  ( int iChan,           /* Channel index in input image. */
    int oChan,           /* Channel index in output image. */
    double brght,        /* Brightness setting ({b}). */
    double ctrst,        /* Contrast setting ({c}). */ 
    sample_uint32_t black, /* Black offset ({k}). */ 
    sample_uint32_t white, /* White limit ({w}). */ 
    sample_uint32_t imin,  /* Minimum integer sample seen. */
    sample_uint32_t imax,  /* Maximum integer sample seen. */
    int clo,             /* Count of underexposed pixels. */
    int chi,             /* Count of overexposed pixels. */
    float vmin,          /* Minimum finite float sample seen. */
    float vmax           /* Maximum finite float sample seen. */
  );
  /* Prints statistics for floatizing channel {iChan} of a PGM/PPM image
    into channel {oChan} of a float image. */

interval_t sample_conv_hdyn_merge_intervals(int n, interval_t fv[]);
  /* Combines a list of interval estimates {fv[0..n-1]} for the brightness
    of a pixel into a single interval {rv}. 
    
    Assumes that a finite interval {fv[k]} is actually a gaussian
    distribution with mean {mid(fv[k])} and standard deviation
    {rad(fv[k])/3}. For certain input data, the returned interval may
    be semi-infinite, infinite, or empty. */

sample_uint32_t sample_conv_hdyn_quantize
  ( float fv, 
    double brght,          /* Brightness setting ({b}). */
    double ctrst,          /* Contrast setting ({c}). */ 
    double sigma,          /* Noise level ({s}). */
    double gamma,          /* Power law exponent ({g}). */ 
    sample_uint32_t black,   /* Black offset ({k}). */ 
    sample_uint32_t white,   /* White limit ({w}). */ 
    /* Statistical accumulators: */
    float *vmin,
    float *vmax, 
    int *clo,
    int *chi,
    sample_uint32_t *imin, 
    sample_uint32_t *imax
  );
  /* Converts a float sample {fv = Z} to an integer sample {iv = Y}
    according to the formula {sample_conv_hdyn_INFO}. Also updates the
    range {vmin,vmax} to enclose the input value {fv}, and the range
    {imin,imax} to enclose the output value {iv}. Increments {clo} or
    {chi} if the computed pixel is underexposed or overexposed,
    respectively. */

void sample_conv_hdyn_print_quantize_stats
  ( int iChan,           /* Channel index in input image. */
    int oChan,           /* Channel index in output image. */
    double brght,        /* Brightness setting ({b}). */
    double ctrst,        /* Contrast setting ({c}). */ 
    sample_uint32_t black, /* Black offset ({k}). */ 
    sample_uint32_t white, /* White limit ({w}). */ 
    float vmin,          /* Minimum float sample seen. */
    float vmax,          /* Maximum float sample seen. */
    int clo,             /* Number of samples seen below {lo}. */
    int chi,             /* Number of samples seen above {hi}. */
    sample_uint32_t imin,  /* Minimum integer sample seen. */
    sample_uint32_t imax   /* Maximum integer sample seen. */
  );
  /* Prints statistics for quantizing channel {iChan} of a {sample_conv_hdyn_t}
    into channel {oChan} of a PGM/PPM image. */

#endif
