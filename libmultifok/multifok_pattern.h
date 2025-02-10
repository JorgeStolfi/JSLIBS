/* 3D patterns for texturizing objects. */
/* Last edited on 2025-02-08 17:28:18 by stolfi */

#ifndef multifok_pattern_H
#define multifok_pattern_H

#include <stdint.h>

#include <r3.h>
#include <bool.h>
#include <interval.h>
#include <frgb.h>
    
typedef double multifok_pattern_double_proc_t(r3_t *q);
  /* Type of a function that computes a grayscale value as a function of 
    the point {q}. */
    
typedef frgb_t multifok_pattern_frgb_proc_t(r3_t *q);
  /* Type of a function that computes an {frgb_t} tuple as a function of 
    the point {q}. */

#endif
