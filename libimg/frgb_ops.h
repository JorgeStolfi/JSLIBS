/* frgb_ops.h - basic operations on colors. */
/* Last edited on 2024-12-20 17:34:57 by stolfi */

#ifndef frgb_ops_H
#define frgb_ops_H

#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include <argparser.h>
#include <bool.h>

#include <frgb.h>

frgb_t frgb_mix(double ca, frgb_t *a, double cb, frgb_t *b);
  /* Computes the linear combination {ca*a + cb*b}. */

frgb_t frgb_scale(double s, frgb_t *a);
  /* Multiplies all components of {a} by {s}. */

frgb_t frgb_shift(double d, frgb_t *a);
  /* Adds {d} to all components of {a}. */

frgb_t frgb_add(frgb_t *a, frgb_t *b);
frgb_t frgb_sub(frgb_t *a, frgb_t *b);
frgb_t frgb_mul(frgb_t *a, frgb_t *b);
  /* Evaluates the specified operation to {a} and {b}, componentwise. */

bool_t frgb_is_all_zeros(frgb_t *a);
  /* TRUE iff {a} is the triple {(0,0,0)}. */

bool_t frgb_is_all_ones(frgb_t *a);
  /* TRUE iff {a} is the triple {(1,1,1)}. */

bool_t frgb_eq(frgb_t *a, frgb_t *b);
  /* TRUE iff colors {a} and {b} are identical. */

double frgb_gamma_encoding_gray(double y, double expo_dec, double bias);
  /* Applies gamma encoding to the intensity {y}, roughly
    {sgn(y)*|y|^(1/expo_dec)}. Equivalent to
    {sample_conv_gamma(y,1/expo_dec,bias)}.
    
    Assumes that {y} is in {[-1 _ +1]} and that {expo_dec}
    is positive.  Typically used when {y} is a computed 
    physical intensity that is about to be converted to 
    an integer for output to an image file. */
		
double frgb_gamma_decoding_gray(double y, double expo_dec, double bias);
  /* Undoes the gamma encoding of the intensity {y}, roughly
    {sgn(y)*|y|^(expo_dec)}. Equivalent to
    {sample_conv_gamma(y,expo_dec,bias)}.
    
    Assumes that {expo_dec} is positive and {y} is in
    {[-1 _ +1]}. Typically used when {y} is an integer sample read
    from an image file that was linearly scaled to {[0-1]}. */

typedef frgb_t frgb_adjuster_t(frgb_t *p, int32_t col, int32_t row);
  /* A function called by other procs to adjust a user-supplied 
     color argument {p}, which presumably applies to pixel 
     {col,row}, before internal use. (A typical use of this
     proc is to apply gamma-correction on argument colors.) */

frgb_t frgb_correct_arg(frgb_t *p, frgb_t *inGamma, int32_t gray);
  /* First computes a new triplet {x} by applying {frgb_undo_gamma} to
    each component of the RGB triplet {*p} with {gamma} set to the
    corresponding component of {inGamma}. If {inGamma == NULL},
    assumes {inGamma = (1,1,1)} which means no gamma correction. Then,
    if {gray = FALSE}, returns {x}; if {gray = TRUE}, returns a gray
    triplet with same luminosity as {x}. */

/* Minimum intensity for logarithmic computations: */
#define VAL_EPS 1.0e-6

double frgb_log_scale_gray(double x);
  /* Computes the logarithm of {x}, after ensuring that it 
    is not less than {VAL_EPS}. */

void frgb_log_scale(frgb_t *p, int32_t chns);
  /* Computes the logarithm of {p[0..chns-1]}, 
    after ensuring that it is not less than {VAL_EPS}. */

double frgb_clip_gray(double p);
  /* Clips {p} to the unit interval. */

void frgb_clip_rgb(frgb_t *p);
  /* Clips {*p} to the unit cube, preserving its luminosity (if
    possible) and its hue. If the luminosity is greater than 1, the
    result is white. If the luminosity is less than 0, the result is
    black. */

void frgb_clip_rgb_towards(frgb_t *p, frgb_t *q);
  /* Clips {*p} to the unit cube, by moving it towards
    the color {*q} (which must be strictly inside the unit
    cube) if necessary. */

void frgb_clip_rgb_towards_grey(frgb_t *p);
  /* Clips {*p} to the unit cube, by moving it towards
    middle gray {(0.5,0.5,0.5)} if necessary. */

double frgb_apply_kappa_gray(double y, double kappa);
  /* Applies a nonlinear correction {f(y,kappa)} to the intensity {y}.
    For {y} in {[0 _ 1]}, {f(y,k) = k*y/((k-1)*y + 1)}.
    Note that  {f(y,1)=y}, {f'(0,k) = k}, and {f(f(y,k),1/k) = y}.
    Outside that domain, {f(y,k)} is 0 for {y<0} and 1 for {y>1}. */

void frgb_apply_glob_kappa_sat_clip(frgb_t *p, double kap, double satf);
  /* Applies the global brightness (`kappa') correction to 
    pixel {p.c[0..3]}, then adjusts the saturation of {p} 
    by the factor {satf}, and clips it to the [0_1]^3 cube:  */

int32_t frgb_dequal(double *a, double *b, int32_t chns);
int32_t frgb_fequal(float *a, float *b, int32_t chns);
  /* TRUE iff {a[i]} equals {b[i]} for {i} in {0..chns-1}. */

double frgb_floatize(int32_t ival, int32_t maxval, double zero, double scale);
  /* Converts an integer pixel value {ival} (read from the input file) to
    an interval of floating-point intensities. Maps pixel value {zero}
    to 0, {zero+scale} to 1. */

int32_t frgb_quantize(double fval, double zero, double scale, int32_t maxval);
  /* Converts an interval of floating-point intensities to an integer 
    intensity value, suitable for output.  Maps 0 to {zero}, 1 to 
    {zero+scale}, clipping the result to the interval {[0..maxval]}. */
    
/* PARSING TRIPLETS */
  
#define frgb_parse_HELP \
  "{R_VALUE} {G_VALUE} {B_VALUE}"
  
#define frgb_parse_INFO \
  "three real numbers.  E.g.\n" \
  "\n" \
  "      \"1 1 1\"\n" \
  "      \"1.000 0.500 0.600\""

frgb_t frgb_parse(argparser_t *pp, double lo, double hi);
 /* Parses three consecutive numbers from the command line,
   as described by {frgb_parse_INFO}. Checks whether each
   component lies in {[lo _ hi]}.  Increments {*argn}. */
 
frgb_t frgb_read(FILE *rd, double lo, double hi);
 /* Parses three consecutive numbers from file {rd},
   checks whether each component lies in {[lo _ hi]}. */
 
/* PARSING COLOR VALUES */
  
#define frgb_parse_color_HELP \
  "{R_VALUE} {G_VALUE} {B_VALUE} [ / {DENOM} ]"
  
#define frgb_parse_color_INFO \
  "three numbers, optionally" \
  " followed by a slash (\"/\") and a common denominator.  All these" \
  " tokens (including the slash) should be separated by spaces.  E.g.\n" \
  "\n" \
  "      \"1 1 1\"\n" \
  "      \"1.000 0.500 0.600\"\n" \
  "      \"75.0 80.2 57.1 / 100\" (same as \"0.750 0.802 0.571\")\n" \
  "      \"128 196 255 / 255\""

frgb_t frgb_parse_color(argparser_t *pp);
  /* Parses a color specification from the command line,
    as described by {frgb_parse_color_INFO}. If a denominator
    is present, divides it into all three components. Returns the color
    as three floats. Increments {*argn}. */
    
frgb_t frgb_read_color(FILE *rd);
  /* Parses a color specification from file {rd},
    as described by {frgb_parse_color_INFO}. 
    Returns the color as three doubles in {[0 _ 1]}. */
    
/* COLORSPACE CONVERSIONS */

/* The following XYZ coordinates of the RGB primaries are used by the procedure 
  {frgb_to_CIE_XYZrec601_1} below. Not clear where these numbers
  came from. They are claimed to use the "European TV RGB standard.
  according to CIE XYZ Rec. 601-1" */

#define frgb_YR (+0.298911)
#define frgb_YG (+0.586611)
#define frgb_YB (+0.114478)
  
#define frgb_XR (+0.606881)
#define frgb_XG (+0.173505)
#define frgb_XB (+0.200336)

#define frgb_ZR (00.000000)
#define frgb_ZG (+0.066097)
#define frgb_ZB (+1.116157)

void frgb_to_CIE_XYZrec601_1(frgb_t *p);
void frgb_from_CIE_XYZrec601_1(frgb_t *p);
double frgb_luminance_CIE_XYZrec601_1(frgb_t *p);
  /* Converts the RGB triple {p} to/from CIE XYZ Rec. 601-1 coordinates. */

void frgb_to_CIE_XYZccir709(frgb_t *p);
void frgb_from_CIE_XYZccir709(frgb_t *p);
double frgb_luminance_CIE_XYZccir709(frgb_t *p);
  /* Converts the RGB triple {p} to/from CIE XYZ CCIR 709 coordinates. */

void frgb_to_CIE_XYZitu_D65(frgb_t *p);
void frgb_from_CIE_XYZitu_D65(frgb_t *p);
double frgb_luminance_CIE_XYZitu_D65(frgb_t *p);
  /* Converts the RGB triple {p} to/from CIE XYZ ITU (D65) coordinates. */
     
void frgb_to_YUV(frgb_t *p);
void frgb_from_YUV(frgb_t *p);
  /* Converts the RGB triple {p} to/from European TV YUV coordinates. */

double frgb_get_Y(frgb_t *p);
  /* The luminance of the RGB triple {p} (the Y coord
    of European TV YUV coords). */

double frgb_get_Y_pbm(frgb_t *p);
  /* The luminance of the RGB triple {p} by the formula
    used in the PBMplus package (almost, but not exactly,
    the same as {frgb_get_Y}). */

void frgb_to_HSV_CG(frgb_t *p);
void frgb_from_HSV_CG(frgb_t *p);
  /* Converts the RGB triple {p} to/from the HSV (hue, saturation,
    value) system used in computer graphics (A. R. Smith,
    SIGGRAPH'78). Namely the hue H is a number in the range [0_1],
    defined for the primary and secondary colors as follows
    
      {(1,0,0) -> 0/6 = 0.0000}, 
      {(1,1,0) -> 1/6 = 0.1667},
      {(0,1,0) -> 2/6 = 0.3333},
      {(0,1,1) -> 3/6 = 0.5000},
      {(0,0,1) -> 4/6 = 0.6667},
      {(1,0,1) -> 5/6 = 0.8333},
      {(1,0,0) -> 6/6 = 1.0000} (in the limit from the magenta side),
    
    
    Any `pure' color (one whose min and max coordinates are 0 and 1,
    respectively) lies on some edge of the color cube connecting those
    six colors; its hue is then defined by linear interpolation. By
    definition, the hue is not changed by mixture with white or black.
    The hue of black, white, and grays is undefined; the function
    returns an arbitrary value in that case.
    
    The saturation S is the percentage of `pure' color of that hue
    that, mixed with the proper shade of gray, gives {p}. The value V
    is the maximum of the R,G, and B coordinates. */

double frgb_get_H(frgb_t *p);
  /* The hue of the RGB triple {p}, as implied by 
    its UV chroma coordinates in the European TV standard. The
    hue {H} is direction of the vector {(U,V)}, converted from the
    range {[0 _ 2*PI)} to {[0 _ 1)} with a shift so that 
    red {(1,0,0)} has {H = 0}. If {U} and {V} are 0 the
    result is 0. */

double frgb_H_from_UV(double U, double V);
  /* Same as {frgb_get_H} but from the given European TV {U,V}
    coordinates.  Only the direction of {(U,V)} matters. */

void frgb_H_to_uv(double H, double *u, double *v);
  /* The partial inverse of {frgb_H_from_UV}. The result {u,v} is a vector {U,V}
    with hue {H}, but  normalized so that {u^2+v^2=1}. */

double frgb_T_from_YUV(double Y, double U, double V);
  /* Computes the relative saturation {T} of a color {p} given its YUV
    coordinates {Y,U,V}.
    
    The relative saturation {T} of a color {p} is the linear position of 
    {p} along the segment of colors with the same luminance
    that goes trough {p} and extends from the gray diagonal
    to the boundary of the unit RGB cube. Thus {T} is 0 for a
    gray color, and 1 for any color on the boundary of the RGB
    cube. 
    
    By convention, {T} is zero if {Y} is outside the interval {[0_1]}.
    
    The relative saturation {T} is a continuous function of the RGB
    coordinates, but it is not smooth (C1): it has a kink whenever the
    distal end of the segment crosses an edge of the cube. */
    
void frgb_to_HTY(frgb_t *p);
void frgb_from_HTY(frgb_t *p);
  /* Converts the RGB triple {p} to/from the HTY color system, whose
    coordinates are hue {H}, relative saturation {T}, and luminance {Y}.
    The hue is defined by {frgb_get_H}. */

void frgb_YUV_to_YHS(frgb_t *p);
void frgb_YHS_to_YUV(frgb_t *p);
  /* Converts the triple {p} betweem the {YUV}
    and the {YHS} coordinate systems.  In the latter, the {Y} coordinate 
    is the same as in {YUV}; the {H} coordinate is given by {frgb_H_from_UV(U,V)},
    and the "Euclidean" saturation {S} is {hypot(U,V)/Y}. 
    
    Note that {S} is not silply related to the relative saturation {T} of the {HTY} 
    system. In particular, a color with {S==1} may liee inside or outside the {RGB}
    unit cube. */

void frgb_to_YIQ(frgb_t *p);
void frgb_from_YIQ(frgb_t *p);
  /* Converts the RGB triple {p} to/from American TV YIQ coordinates. */

void frgb_to_YCbCr_601_1(frgb_t *p);
void frgb_from_YCbCr_601_1(frgb_t *p);
  /* Converts the RGB triple {p} to/from Y/Cb/Cr coordinates (JFIF?). */

void frgb_to_YUV_a(frgb_t *p);
void frgb_from_YUV_a(frgb_t *p);
  /* Converts the RGB triple {p} to/from alternate YUV coords,
    used by old {frgb_to_yuv}. */

void frgb_to_YUV_b(frgb_t *p);
void frgb_from_YUV_b(frgb_t *p);
  /* Converts the RGB triple {p} to/from alternate YUV coords,
    used by old {ppmtoyuvx}. */

void frgb_YUV_to_yuv(frgb_t *p, double ybias);
void frgb_YUV_from_yuv(frgb_t *p, double ybias);
  /* Maps any luminance/chroma linear coordinates through
    a projective transformation that provides more
    visually accurate distances between colors. */

void frgb_YUV_to_Yuv(frgb_t *p, double ybias);
void frgb_YUV_from_Yuv(frgb_t *p, double ybias);
  /* Like {frgb_YUV_to_yuv} and {frgb_YUV_from_yuv}, but 
    applies the transformation to the U and V coordinates
    only, leaving {Y} unchanged. */

/* PRINTING */

void frgb_print(FILE *f, char *pref, frgb_t *p, int32_t chns, char *fmt, char *suff);
  /* Prints the first {chns} channels of {*p} to file {f}, surrounded
    by the given {pref} and {suff} strings. Each component is printed
    with the format {fmt}, which should be suitable for a float
    value. */

void frgb_print_int_pixel(FILE *f, char *pref, int32_t *p, int32_t chns, char *suff);
  /* Prints {p[0..chns-1]} to file {f}, surrounded by the
    given {pref} and {suff} strings.  (This function should be moved
    to some other interface.) */

/* DEBUGGING */

extern int32_t frgb_DEBUG;    /* Set this variable to activate the following procs. */

void frgb_debug(char *label, int32_t col, int32_t row, frgb_t *p, int32_t chns, char *tail);
  /* Prints {col}, {row}, and the first {chns} channels of {*p} to
    {stderr}, if the global flag {frgb_DEBUG} is true. The printout is
    preceded by the string {label} and followed by {tail}. */

void frgb_debug_int_pixel(char *label, int32_t col, int32_t row, int32_t *p, int32_t chns, char *tail);
  /* Prints {col}, {row}, and {p[0..chns-1]} to {stderr}, if the global flag 
    {frgb_DEBUG} is true.  The printout is preceded by the string {label}
    and followed by {tail}.  (This function should be moved
    to some other interface.)  */

#endif
