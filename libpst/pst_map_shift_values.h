#ifndef pst_map_shift_values_H
#define pst_map_shift_values_H

/* pst_map_shift_values.h -- procedures for working with height maps. */
/* Last edited on 2025-02-27 12:08:26 by stolfi */

#include <stdint.h>

#include <sign.h>
#include <float_image.h>

void pst_map_shift_values(float_image_t *A, int32_t wch, sign_t dir, double shift[]);
  /* Adds to every sample of every channel {c}, except {wch},
    the amount {dir*shift[c]}. 
    
    If {wch} is a valid channel index of {A}, every sample {A[wch,X,Y]}
    is interpreted as a reliability weight for the samples {A[c,X,Y]} in
    channels {c != wch}. If any input or output sample {A[c,X,Y]} is not
    finite or its weight is zero, then {A[wch,X,Y]} (if it exists) is
    set to zero and other samples {A[c,X,Y]} of that pixels are set to
    {NAN}. Otherwise the weights are not affected. */
 
#endif
