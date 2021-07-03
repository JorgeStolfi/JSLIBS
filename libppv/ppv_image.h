/* ppv_image.h --- converting {ppv_array_t} to/from {uint16_image_t}. */
/* Last edited on 2021-07-03 14:53:21 by jstolfi */

#ifndef ppv_image_H
#define ppv_image_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <ppv_types.h>
#include <ppv_array.h>
#include <uint16_image.h>

uint16_image_t *ppv_image_from_array(ppv_array_t *A, ppv_sample_t maxval);
  /* Converts {A} into a {uint16_image_t} {J}.  
    
    The array {A} must have {A->d == 2} or {A->d == 3}.
    The first two axes of {A} are interpreted as column 
    and row indices of {J}.  If {A->d} is 3, the third axis
    is interpreted as color channel index; otherwise the 
    image will have a single channel.
    
    The samples of {A} must be in the range {0..maxval}, and
    {maxval} must not exceed {uint16_image_MAX_SAMPLE}. */

ppv_array_t *ppv_image_to_array(uint16_image_t *J);
  /* Converts {J} into a {ppv_array_t} {A}.  
    
    If {J} has a single color channel, {A->d} will be 2, and the two
    axes will be will be the column and row indices of {J}. If {J} has a
    two or more channels. {A->d} will be 3, and the third axis will be
    the channel index.
    
    The samples are copied without change. The bit count {A->bps} will
    be the smallest value such that {2^A->bps > J.maxval}.
    Note that the {maxval} information will be lost 
    if {J.maxval} is not of the form {2^k-1}. */

#endif

