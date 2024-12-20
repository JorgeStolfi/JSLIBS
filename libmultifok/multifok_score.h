/* Focus detector for multi-focus stereo. */
/* Last edited on 2024-12-05 10:36:58 by stolfi */

#ifndef multifok_score_H
#define multifok_score_H

#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_term.h>

/* 
  FOCUS DETECTOR OPERATORS
  
  A /focus sharpness score/ is a local image operator that estimates the
  sharpness of the image at the applied pixel.  Ideally it is close to 0
  for a very unfocused (blurred) image, and close to 1 for a well-focused (sharp) one.
  
  In practice, it is sufficient that, when applied at the same point of 
  the same image blurred by various amounts, it has the maximum value for 
  the version with least blurring.
  
  The focus detector operators in this module work on a window of {NW}
  by {NW} samples, stored in an array {x[0..NS-1]} of {float}s, where {NS
  = NW*NW}, linearized by rows.  See {multifok_window.h} for 
  mreo details. */ 

typedef double multifok_score_op_t(int32_t NC, double wt[]);
  /* Given an array of window samples ?? */ 

double multifok_score_from_terms(int32_t NT, double wt[], double term[]);
  /* Computes the focus sharpness score as the weighted sum 
    {SUM{wt[k]*term[k] : k \in 0..NT-1 }}. */ 

#endif
