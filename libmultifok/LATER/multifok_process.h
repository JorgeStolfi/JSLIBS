/* Process a stack of images estimating the sharp image and scene depth. */
/* Last edited on 2024-12-05 10:35:51 by stolfi */

#ifndef multifok_process_H
#define multifok_process_H

#include <stdint.h>

#include <float_image.h>

#include <multifok_stack.h>
#include <multifok_result.h>

multifok_result_t *multifok_process
  ( multifok_stack_t *stack,
    multifok_score_op_t *score,
    float_image_t *
  );
  /* Processes a stack of multi-focus images
    

#endif
