/* See {multifok_image_stack.h}. */
/* Last edited on 2024-08-02 15:44:38 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <multifok_image_stack.h>
  
#define DASHES "------------------------------------------------------------"
     
multifok_image_stack_t* mfmi_read_stack
  ( char *inPrefix,
    int32_t NI,
    double *zFoc
  )
  {
    
    float_image_t *csimg = talloc(NI, float_image_t*);
    float_image_t *shimg = talloc(NI, float_image_t*);
    
    int32_t NC, NX, NY; /* Image dimensions. */
 
    for (int32_t ki = 0; ki < NI; ki++)
      { char *zTag = NULL;
        asprintf(&zTag, "-fd%08.4f-z%08.4f", zDep, zFoc[ki]);
        csimg[ki] = multifok_image_stack_read_color_image(inPrefix, zTag); 
        if (ki == 0)
          { float_image_get_size(csimg[ki], &NC, &NX, &NY); }
        else
          { float_image_check_size(csimg[ki], NC, NX, NY); }

        shimg[ki] = multifok_image_stack_read_sharpness_image(inPrefix, zTag); 
        float_image_check_size(shimg[ki], 1, NX, NY);
      } 

    mfok_image_stack_t *stack = talloc(1, mfok_image_stack_t);
    stack->NI = NI;
    stack->NC = NC;
    stack->NX = NX;
    stack->NY = NY;
    stack->csimg = csimg,
    stack->shimg = shimg;
    stack->zFoc = zFoc
    
    return stack;
  }
    /* Reads a stack of scene images {csimg[0..NI-1]} and the corresponding 
    sharpness score maps {shimg[0..NI-1]}, from files 
    "{inprefix}-???", and resturns them in a {multifok_image_stack_t}
    object {stack}.  Also sets {stac.zFoc} to the vector {zFoc},
    assumed to have {NI} elements. */
 
#define multifok_image_stack_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

