#ifndef frgb_H
#define frgb_H

/* frgb.h - floating-point RGB color data type */
/* Last edited on 2024-12-28 18:13:18 by stolfi */

#include <stdint.h>
#include <limits.h>

/* Channels in a color image: */
#define frgb_CHANNELS 3

typedef struct frgb_t { float c[frgb_CHANNELS]; } frgb_t;
  /* An color value represented as intensities in three channels.
    
    The three components usually specify the intensities of 
    red, green, and blue primaries that are to be mixed to give the 
    desired color.  Typically 0 means no amount of the primary, and 
    1 means the maximum allowed in the situation.  
    
    However, the colors of the primaries themselves are 
    not specified by this interface and should be defined by the
    context.  Also the intensities may be linear (proportional 
    to physical luminous power) or encoded in some non-linear
    scale, such as gamma-encoding; or may be differences, or 
    points in some other color space, scuh as CIE XYZ or YCbCr.
    As such they may be negative or greater than 1.  These 
    details too should be determined by the context of use.
    
    See {frgb_ops.h} for functions that convert {frgb_t} values
    between different encodings and some standard color spaces. */

#ifndef INF
#define INF INFINITY 
#endif

#define frgb_Zeros   (frgb_t){{0.0, 0.0, 0.0}}
#define frgb_Black   (frgb_t){{0.0, 0.0, 0.0}}
#define frgb_White   (frgb_t){{1.0, 1.0, 1.0}}
#define frgb_Ones    (frgb_t){{1.0, 1.0, 1.0}}
#define frgb_NoColor (frgb_t){{NAN, NAN, NAN}}

#endif
