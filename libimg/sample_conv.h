#ifndef sample_conv_H
#define sample_conv_H

/* {sample_conv.h} - conversion between floating-point and integer samples. */
/* Last edited on 2023-01-14 12:04:13 by stolfi */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <bool.h>

typedef uint32_t sample_uint32_t;
  /* A quantized sample value. */

float sample_conv_encode_BT709(float Y);
  /* Applies the luminance encoding (a.k.a. `gamma correction')
    according to ITU-R Recommendation BT.709. Namely, returns the
    sample value {V} that must be stored in an image in order to
    produce intensity {Y} on the idealized monitor assumed by BT.709.
    If {Y} is {NAN} or {±INF}, returns {Y} itself.
    See {sample_conv_BT709_INFO} above for details. */

float sample_conv_decode_BT709(float V);
  /* The inverse of {sample_conv_encode_BT709}. Namely, returns the 
    intensity {Y} produced on the idealized monitor assumed by BT.709
    when displaying a sample value {V} taken from an image.
    If {V} is {NAN} or {±INF}, returns {V} itself.
    See {sample_conv_BT709_INFO} above for details. */

/* !!! Should change my gamma+bias encoding so that BT709 is a special case. !!! */

#define sample_conv_BT709_INFO \
  "The ITU-R BT709 encoding and decoding functions are" \
  " strictly monotonic for all arguments, positive or" \
  " negative, and map the values {-1}, 0, and {+1} to" \
  " themselves.\n" \
  "\n" \
  "  The encoding map is approximately {V = Y^0.45} for" \
  " positive {Y}, and {V = -((-Y)^0.45)} for negative {Y}; but" \
  " is replaced by a linear function {V = 4.5*Y} near" \
  " zero (for {|Y| < 0.01805}).\n" \
  "\n" \
  "  The decoding map is" \
  " approximately {Y = V^(1/0.45) = V^(2.222...)} for" \
  " positive {V},  and {V = -((-V)^(1/0.45))} for" \
  " negative {V}; but is replaced by a linear" \
  " function {V = 0.2222*z} near" \
  " zero (for {|V| < 0.08124})."

float sample_conv_gamma(float z, double gamma, double bias);
  /* Applies a modified power-law correction, with exponent
    {gamma} and offset {bias}, to sample value {z}.
    
    If {z} is {NAN} or {±INF}, returns {z} itself.
    See {sample_conv_gamma_INFO} for details. See
    {sample_conv_gamma_BT709_equiv_INFO} for the relationship between
    {sample_conv_gamma}, {sample_conv_encode_BT709}, and
    {sample_conv_decode_BT709}.
    
    !!! Should take and return a {double} rather than {float}. !!!
  */

#define sample_conv_gamma_INFO \
  "The sample encoding function {sample_conv_gamma} depends" \
  " on two parameters, {gamma} (which must be positive) and" \
  " {bias} (which must be between 0 and 1). The function" \
  " is strictly monotonic for any {gamma} and {bias} and for" \
  " all arguments, positive or negative; and takes the" \
  " values {-1}, {0}, and {+1} to themselves.\n" \
  "\n" \
  "  If {gamma} is 1, the function is the identity, for" \
  " any {bias}.  If the parameter {bias} is zero, the function is a" \
  " simple power-law encoding that maps {V} to {V^gamma} for" \
  " positive {V}, and to {-((-V)^gamma)} for negative {V}.  If {bias}" \
  " is positive, affine  corrections are applied before" \
  " and after the power-law map so that the slope at" \
  " the origin is neither zero or infinity.  To get the" \
  " inverse mapping, use {1/gamma} instead of {gamma}, with" \
  " the same offset {bias}." 
  
#define sample_conv_gamma_BT709_equiv_INFO \
  "The combination {gamma=0.450} and {bias=0.0327} provides" \
  " a good approximation to the ITU-R BT.709 encoding" \
  " function, with maxmum discrepancy {-0.0165} in the" \
  " low-intensity range and about {+0.002} in the" \
  " mid-to-high range.\n" \
  "\n" \
  "  The combination {gamma=1/0.450} and {bias=0.0327} provides" \
  " a good approximation to the ITU-R BT.709 decoding function, with" \
  " maxmum discrepancy {+0.00365} overall."

#define sample_conv_BT709_ENC_GAMMA (0.450)
#define sample_conv_BT709_DEC_GAMMA (1/0.450)
  /* The values of {gamma} that make {sample_conv_gamma}
    approximate the ITU-R BT.709 sample encoding and decoding
    functions, when used with {sample_conv_BT709_BIAS}. 
    See {sample_conv_gamma_BT709_equiv_INFO} for details. */

#define sample_conv_BT709_BIAS (0.0327)
  /* The value of {bias} that makes {sample_conv_gamma}
    approximate the ITU-R BT.709 sample encoding and decoding
    functions, when used with {sample_conv_BT709_ENC_GAMMA} and
    {sample_conv_BT709_DEC_GAMMA} respectively. 
    See {sample_conv_gamma_BT709_equiv_INFO} for details. */

float sample_conv_log(float u, double uref, double logBase);
  /* Converts {u} from linear to logarithmic scale, relative to the
    reference value {uref} and the base {exp(logBase)}. In particular,
    {logBase == 1} gives natural logarithms, {logBase == M_LOG2} gives
    result in octaves, {logBase == M_LOG10} gives result in decades, etc.
    
    More precisely, returns {-INF} if {u} is zero, {+INF} if {u == +INF}, {NAN} if {u} is negative
    or {NAN}, and log(u/uref)/logBase} otherwise.
    
    Requires {uref} to be finite and positive, and {logBase}
    to be finite and nonzero; otherwise returns {NAN} for any {u}. */

float sample_conv_undo_log(float u, double uref, double logBase);
  /* Connverts {u} from log scale to linear scale. The inverse of {sample_conv_log( */
    
float sample_conv_interp(float u, int np, double U[], double V[]);
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
  ( int iChan,           /* Channel index in input image. */
    int oChan,           /* Channel index in output image. */
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
    int *clo,             /* (IN/OUT) Count of input values below {lo}, or NULL. */
    int *chi,             /* (IN/OUT) Count of input values above {hi}, or NULL. */
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
  ( int iChan,           /* Channel index in input image. */
    int oChan,           /* Channel index in output image. */
    float vmin,          /* Minimum float sample seen. */
    float vmax,          /* Maximum float sample seen. */
    double lo,           /* Low end of float scaling range. */
    double hi,           /* High end of float scaling range. */
    int clo,             /* Number of samples seen below {lo}. */
    int chi,             /* Number of samples seen above {hi}. */
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
