/* Process a stack of images estimating the sharp image and scene depth. */
/* Last edited on 2024-08-02 15:55:18 by stolfi */

#ifndef multifok_process_H
#define multifok_process_H

#define _GNU_SOURCE
#include <stdint.h>

#include <float_image.h>

#include <multifok_image_stack.h>
#include <multifok_result.h>

multifok_result_t *multifok_process
  ( multifok_image_stack_t *stack,
    multifok_score_op_t *score,
    float_image_t *
  );
  /* Processes a stack of multi-focus images
    

#endif
