#ifndef sample_conv_H
#define sample_conv_H

/* {sample_conv.h} - conversion between floating-point and integer samples. */
/* Last edited on 2024-12-18 22:49:54 by stolfi */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <bool.h>

typedef uint32_t sample_uint32_t;
  /* A quantized sample value. */

/* THE "BTU709" ENCODING/DECODING */

float sample_conv_BT709_encode(float X);
  /* Applies the luminance encoding (a.k.a. `gamma correction')
    according to ITU-R Recommendation BT.709. Namely, returns the
    sample value {V} that must be stored in an image in order to
    produce intensity {X} on the idealized monitor assumed by BT.709.
    If {X} is {NAN} or {±INF}, returns {X} itself.
    See {sample_conv_BT709_INFO} above for details. */

float sample_conv_BT709_decode(float V);
  /* The inverse of {sample_conv_encode_BT709}. Namely, returns the 
    intensity {X} produced on the idealized monitor assumed by BT.709
    when displaying a sample value {V} taken from an image.
    If {V} is {NAN} or {±INF}, returns {V} itself.
    See {sample_conv_BT709_INFO} above for details. */

#define sample_conv_BT709_INFO \
  "The ITU-R BT709 encoding and decoding functions are" \
  " strictly monotonic for all arguments, positive or" \
  " negative, and map the values {-1}, 0, and {+1} to" \
  " themselves.\n" \
  "\n" \
  "  The encoding map {V = E(X)} is defined as the linear" \
  " function {V = 4.5*X} for {X} between 0 and {X0}, and" \
  " {V = (1+CD)*X^0.45 - CD} for {X > X0}; where {X0 = 0.018} and" \
  " {CD = 0.099}.\n" \
  "\n" \
  "  The decoding map {X=D(V)} is  defined as the linear" \
  " function {X = V/4.5} for {V} between 0 and {V0}, and" \
  " {X = ((V + CE)/(1 + CE))V^(1/0.45)} for {V > V0}; where {V0 = 0.081} and" \
  " {CE = 0.0992965876298}.\n" \
  "\n" \
  "  These functions are extended to" \
  " negative arguments by the equations" \
  " {E(-X) = -E(X)} and {D(-V) = -D(V)}.\n" \
  "\n" \
  "  There are slight" \
  " discontinuities in {D} and {E} at the transition points, and they" \
  " are not exactly the inverses of each other."
    
/* THE "sRGB" ENCODING/DECODING */

float sample_conv_sRGB_encode(float X);
  /* Applies the luminance encoding (a.k.a. `gamma correction')
    according to the sRGB standard (IEC 61966-2-1:1999) Namely, returns the
    sample value {V} that must be stored in an image in order to
    produce intensity {X} on the idealized monitor assumed by the sRGB standard.
    If {X} is {NAN} or {±INF}, returns {X} itself.
    See {sample_conv_sRGB_INFO} above for details. */

float sample_conv_sRGB_decode(float V);
  /* The inverse of {sample_conv_encode_sRGB}. Namely, returns the 
    intensity {X} produced on the idealized monitor assumed by the sRGB standard
    when displaying a sample value {V} taken from an image.
    If {V} is {NAN} or {±INF}, returns {V} itself.
    See {sample_conv_sRGB_INFO} above for details. */

#define sample_conv_sRGB_INFO \
  "The sRGB encoding and decoding functions are" \
  " strictly monotonic for all arguments, positive or" \
  " negative, and map the values {-1}, 0, and {+1} to" \
  " themselves.\n" \
  "\n" \
  "  The encoding map {V = E(X)} is\n" \
  " the linear function {V = 12.92*X} for\n" \
  " {X} beteen 0 and {X0}, and {(1 + CE)*X^(1/2.4)} for\n" \
  " {X > X0}; where {X0 = 0.0031308} and {CE = 0.055}.\n" \
  "\n" \
  "  The decoding map  {X = D(V)} is the linear\n" \
  " function {X = V/12.92} for {V} between 0 and {V0}, and\n" \
  " ((V + CD)/(1 + CD))^2.4} for {V > V0}; where\n" \
  " {V0 = 0.04045} and {CD = 0.055}.\n" \
  "\n" \
  "  These functions are extended to" \
  " negative arguments by the equations" \
  " {E(-X) = -E(X)} and {D(-V) = -D(V)}.\n" \
  "\n" \
  "  There are slight" \
  " discontinuities in {D} and {E} at the transition points, and they" \
  " are not exactly the inverses of each other."

float sample_conv_log(float u, double bias, double uref, double logBase);
  /* Converts {u} from linear to logarithmic scale, relative to the
    reference value {uref} and the base {exp(logBase)}. In particular,
    {logBase == 1} gives natural logarithms, {logBase == M_LOG2} gives
    result in octaves, {logBase == M_LOG10} gives result in decades, etc.
    
    More precisely, returns {+INF} if {u == +INF}, {NAN} if {u} is
    negative or {NAN}, and log(hypot(u,bias)/uref)/logBase} otherwise.
    
    In particular, if {bias} and {u} are both zero, returns {-INF}. If
    {bias = uref > 0}, the result will be non-negative, and will be zero
    iff {u} is zero.
    
    Requires {bias} to be finite and non-negative, {uref} to be finite
    and positive, and {logBase} to be finite and nonzero; otherwise
    returns {NAN} for any {u}. */

float sample_conv_undo_log(float u, double bias, double uref, double logBase);
  /* Converts {u} from log scale to linear scale. The inverse of {sample_conv_log}.
    May return {NAN} if {u} is not a valid result of {sample_conv_log} for those
    parameters. */
    
float sample_conv_interp(float u, int32_t np, double U[], double V[]);
  /* Computes a piecewise affine function of {u} defined by the 
    nodal points {(0,0)}, {(U[i],V[i])} for {i = 0..np-1}, and {(1,1)}. 
    
    The values {U[0..np-1]} must be strictly increasing and 
    in {(0 _ 1)}, and ditto for {V[0..np-1]}.  The function 
    is anti-symmetric across zero, i.e. {f(u) = -f(-u)}. */
    
#define sample_conv_0_1_isMask_false_INFO \
  "each integer sample {IV} is assumed to be the result of rounding some float value" \
  " in the interval {[IV_IV+1]}, within the full interval {[0_MAXVAL+1]}.  Therefore, it" \
  " is converted to {(IV+0.5)/(MAXVAL+1)}.  In this case the integer samples 0" \
  " and {MAXVAL} correspond to the float values {0.0+EPS} and" \
  " {1.0-EPS}, respectively, where {EPS=0.5/(MAXVAL+1)}.  This choice is preferrable when the" \
  " float values 0.0 and 1.0 have no special significance, e.g. for properly" \
  " exposed digital photos."

#define sample_conv_0_1_isMask_true_INFO \
  "the integer sample value 0 corresponds to the float value 0.0, and" \
  " {MAXVAL} to 1.0, exactly.   For other float values, the" \
  " correspondence is defined by linear interpolation and" \
  " rounding to the nearest integer value.  This choice is" \
  " appropriate when the the float values 0.0 and 1.0 are especially" \
  " significant and must be preserved when reading or writing; for instance," \
  " when converting opacity masks where the sample values 0 and {MAXVAL} mean" \
  " fully transparent and fully opaque, respectively."
  /* Documentation for the {isMask} parameter of {sample_conv_floatize} and 
    {sample_conv_quantize}, for image processing programs that 
     use float samples in the range {[0_1]} internally. */

float sample_conv_floatize
  ( sample_uint32_t iv,      /* Integer sample value to convert. */
    sample_uint32_t maxval,  /* Max output sample value. */
    bool_t isMask,         /* The precise meaning of integer sample values mean. */
    double lo,             /* Nominal output for input {0}. */
    double hi,             /* Nominal output for input {maxval}. */
    sample_uint32_t *imin,   /* (IN/OUT) Min integer input value seen, or NULL. */
    sample_uint32_t *imax,   /* (IN/OUT) Max integer input value seen, or NULL. */
    float *vmin,           /* (IN/OUT) Min float output value seen, or NULL. */
    float *vmax            /* (IN/OUT) Max float output value seen, or NULL. */
  );
  /* Converts an integer sample {iv} in {0..maxval},
    to a float value {fv} in {[0_1]}, by an affine function.  
    
    If {isMask} is TRUE, the function simply takes sample value 0 to
    {lo}, and {maxval} to {hi}. This mode is appropriate when
    disproportionally many pixels are expected to be equal to {lo} or
    {hi}; for instance, when converting opacity masks where {0} means
    fully transparent and {maxval} means fully opaque.
    
    If {isMask} is FALSE, any integer sample value {iv} is assumed to
    represent the real interval {[iv_iv+1]} within the total range
    {[0_maxval+1]}. Therefore {iv} is first replaced by {iv+0.5}, and
    then the function maps 0 to {lo} and {maxval+1} to {hi}. Note that
    in this mode the input 0 returns {lo + 0.5*d} and {maxval} returns
    {hi - 0.5*d}, where {d = 1/(maxval+1)}. This mode is more
    appropriate when the original real values before quantization were
    uniformly distributed in {[lo_hi]}, with no particular preference
    for the endpoints; for instance, for unprocessed digital
    photographs and document scans.
    
    In either case, the mapping from {iv} to {fv} interpolates and
    extrapolates linearly between the two points above. This is true
    even if {lo==hi} or {lo > hi}, but {maxval} must be positive. Note
    that if {iv} exceeds {maxval}, the output {fv} will lie beyond the {hi}
    limit.  The values {lo} and {hi} should be finite, otherwise the
    result may be {+INF}, {-INF}, or {NAN}.
    
    Normally, if {iv} was the output of {sample_conv_quantize}, the
    {isMask} parameter should be the same that was used in that conversion.
    This will minimize the rounding error of the two conversions.
    
    The procedure also updates the range {imin,imax} (if not NULL) to
    enclose the input value {iv}, and the range {vmin,vmax} (if not
    NULL) to enclose the output value {fv}. */
    
void sample_conv_print_floatize_stats
  ( int32_t iChan,           /* Channel index in input image. */
    int32_t oChan,           /* Channel index in output image. */
    sample_uint32_t imin,       /* Minimum integer sample seen. */
    sample_uint32_t imax,       /* Maximum integer sample seen. */
    sample_uint32_t maxval,     /* Maximum possible integer sample. */
    double lo,           /* Low end of float scaling range. */
    double hi,           /* High end of float scaling range. */
    float vmin,          /* Minimum float sample seen. */
    float vmax           /* Maximum float sample seen. */
  );
  /* Prints statistics for floatizing channel {iChan} of a PGM/PPM image
    into channel {oChan} of a float image. */

sample_uint32_t sample_conv_quantize
  ( float fv,             /* Float sample value to convert. */
    sample_uint32_t maxval, /* Max output sample value. */ 
    bool_t isMask,        /* The precise meaning of integer sample values mean. */
    double lo,            /* Input value to map to 0. */
    double hi,            /* Input value to map to {maxval}. */
    float *vmin,          /* (IN/OUT) Min float input value seen, or NULL. */
    float *vmax,          /* (IN/OUT) Max float input value seen, or NULL. */
    int32_t *clo,             /* (IN/OUT) Count of input values below {lo}, or NULL. */
    int32_t *chi,             /* (IN/OUT) Count of input values above {hi}, or NULL. */
    sample_uint32_t *imin,  /* (IN/OUT) Min output integer value seen, or NULL. */
    sample_uint32_t *imax   /* (IN/OUT) Max output integer value seen, or NULL. */
  );
  /* Converts a float value {fv} in {[0_1]} to an integer 
    sample {iv} in {0..maxval}, by an affine function.  
    
    If {isMask} is TRUE, the function takes the input value {lo} to
    sample value 0, {hi} to {maxval}, and interpolates linearly in
    between, rounding the results to the nearest integer. Halfway
    values are rounded down if {fv} is in the range {[lo_md)}, and up
    if {fv} is in {[md_hi]}, where {md} is the center of {[lo_hi]}.
    This mode is appropriate when disproportionally many pixels are
    expected to be equal to {lo} or {hi}; for instance, when
    quantizing opacity masks where {lo} means fully transparent and
    {hi} means fully opaque.
    
    If {isMask} is FALSE, any integer sample value {iv} is assumed to
    represent the real interval {[iv_iv+1)} within the total range
    {[0_maxval+1)}. Therefore {lo} is mapped to 0, {hi} is mapped to
    {maxval+1}, and the result is truncated down to an integer.
    As a special case, the value {hi} is mapped to {maxval} rather than
    {maxval+1}. This mode is more appropriate
    if the sample {fv} is expected to be smoothly distributed in {[lo_hi]}, 
    with no particular preference for the endpoints; for instance,
    when {fv} is light intensity or other physically-based quantity.
    
    Normally, if {fv} was the output of {sample_conv_floatize}, the
    {isMask} parameter should be the same that was used in that conversion.
    This will reproduce the original integer sample value without 
    additional rounding error.
    
    If {lo==hi}, the output will be {maxval/2}, irrespective of {fv}.
    Otherwise the resulting sample value is always clipped to the
    range {0..maxval}; if {clo} or {chi} are not NULL, they are
    incremented whenever the input value {fv} lies before {lo} or
    beyond {hi} respectively (except that values outside but very
    close to those limits may fail to be counted because of roundoff
    errors). The procedure works even if {lo>hi}, in which case
    {clo} counts input values *greater* than {lo} and {chi}
    counts values *less* than {hi}.
    
    The value {fv} must not be {NAN}. The values {lo} and {hi} should
    be finite, otherwise the procedure may fail.
    
    The procedure also updates the range {vmin,vmax} (if not NULL) to
    enclose the input value {fv}, and the range {imin,imax} (if not
    NULL) to enclose the output value {iv}. */

void sample_conv_print_quantize_stats
  ( int32_t iChan,           /* Channel index in input image. */
    int32_t oChan,           /* Channel index in output image. */
    float vmin,          /* Minimum float sample seen. */
    float vmax,          /* Maximum float sample seen. */
    double lo,           /* Low end of float scaling range. */
    double hi,           /* High end of float scaling range. */
    int32_t clo,             /* Number of samples seen below {lo}. */
    int32_t chi,             /* Number of samples seen above {hi}. */
    sample_uint32_t maxval,     /* Maximum possible integer sample. */
    sample_uint32_t imin,       /* Minimum integer sample seen. */
    sample_uint32_t imax        /* Maximum integer sample seen. */
  );
  /* Prints statistics for quantizing channel {iChan} of a {sample_conv_t}
    into channel {oChan} of a PGM/PPM image. */

void sample_conv_choose_maxval(uint32_t chns, sample_uint32_t imaxval[], sample_uint32_t maxmaxval, sample_uint32_t *omaxvalP);
  /* Assumes that integers {imaxval[0..chns-1]} are 
    the max sample value for channels {0..chns-1} of some quantized image.
    This procedure chooses a single max sample value {omaxval}
    for all channels, so that the scaling of the samples 
    from {0..imaxval[i]} to {0..omaxval} is as precise as possible.
    
    Each {imaxval[i]} must be in {1..maxmaxval}. The chosen {omaxval} is
    returned in {*omaxvalP}. It will be no less than {imaxval[i]}, for
    all {i}, and will not exceed {maxmaxval}. */

#endif
