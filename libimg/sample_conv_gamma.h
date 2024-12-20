#ifndef sample_conv_gamma_H
#define sample_conv_gamma_H

/* {sample_conv_gamma.h} - generic gamma-like encoding/decoding. */
/* Last edited on 2024-12-20 16:59:29 by stolfi */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <bool.h>

float sample_conv_gamma(float z, double expo, double bias);
  /* Applies a modified power-law correction, with exponent
    {expo} and offset {bias}, to sample value {z}.
    
    If {z} is {NAN} or {±INF}, returns {z} itself.
    See {sample_conv_gamma_INFO} for details. See
    {sample_conv_gamma_BT709_equiv_INFO} for the relationship between
    {sample_conv_gamma}, {sample_conv_encode_BT709}, and
    {sample_conv_decode_BT709}.
    
    !!! Should take and return a {double} rather than {float}. !!!
  */

#define sample_conv_gamma_INFO \
  "The sample encoding function {sample_conv_gamma} depends" \
  " on two parameters, {expo} (which must be positive) and" \
  " {bias} (which must be between 0 and 1). The function" \
  " is strictly monotonic for any {expo} and {bias} and for" \
  " all arguments, positive or negative; and takes the" \
  " values {-1}, {0}, and {+1} to themselves.\n" \
  "\n" \
  "  If {expo} is 1, the function is the identity, for" \
  " any {bias}.  If the parameter {bias} is zero, the function is a" \
  " simple power-law encoding that maps {V} to {V^expo} for" \
  " positive {V}, and to {-((-V)^expo)} for negative {V}.  If {bias}" \
  " is positive, affine  corrections are applied before" \
  " and after the power-law map so that the slope at" \
  " the origin is neither zero or infinity.  To get the" \
  " inverse mapping, use {1/expo} instead of {expo}, with" \
  " the same offset {bias}." 
  
#define sample_conv_gamma_BT709_equiv_INFO \
  "The combination {expo=0.45} and {bias=0.0415} provides" \
  " a good approximation to the ITU-R BT.709 encoding" \
  " function, with maxmum discrepancy {0.0097}.\n" \
  "\n" \
  "  The combination {expo=1/0.45} and {bias=0.0415} provides" \
  " a good approximation to the ITU-R BT.709 decoding function, with" \
  " maxmum discrepancy {0.0092}."

#define sample_conv_gamma_BT709_ENC_EXPO (0.45)
#define sample_conv_gamma_BT709_DEC_EXPO (1/0.45)
  /* The values of the exponennt {expo} that make {sample_conv_gamma}
    approximate the ITU-R BT.709 sample encoding and decoding
    functions, when used with {sample_conv_gamma_BT709_BIAS}. 
    See {sample_conv_gamma_BT709_equiv_INFO} for details. */

#define sample_conv_gamma_BT709_BIAS (0.0415)
  /* The value of {bias} that makes {sample_conv_gamma}
    approximate the ITU-R BT.709 sample encoding and decoding
    functions, when used with {sample_conv_gamma_BT709_ENC_EXPO} and
    {sample_conv_gamma_BT709_DEC_EXPO} respectively. 
    See {sample_conv_gamma_BT709_equiv_INFO} for details. */
  
#define sample_conv_gamma_sRGB_equiv_INFO \
  "The combination {expo=0.4} and {bias=0.02} provides" \
  " a good approximation to the IEC sRGB encoding" \
  " function, with maxmum discrepancy {0.0064}.\n" \
  "\n" \
  "  The combination {expo=1/0.4} and {bias=0.02} provides" \
  " a good approximation to the sRGB decoding function, with" \
  " maxmum discrepancy {0.0021}."

#define sample_conv_gamma_sRGB_ENC_EXPO (0.4)
#define sample_conv_gamma_sRGB_DEC_EXPO (1/0.4)
  /* The values of {expo} that make {sample_conv_gamma}
    approximate the IEC sRGB sample encoding and decoding
    functions, when used with {sample_conv_gamma_sRGB_BIAS}. 
    See {sample_conv_gamma_sRGB_equiv_INFO} for details. */

#define sample_conv_gamma_sRGB_BIAS (0.02)
  /* The value of {bias} that makes {sample_conv_gamma}
    approximate the IEC sRGB sample encoding and decoding
    functions, when used with {sample_conv_gamma_sRGB_ENC_EXPO} and
    {sample_conv_gamma_sRGB_DEC_EXPO} respectively. 
    See {sample_conv_gamma_sRGB_equiv_INFO} for details. */

#endif
