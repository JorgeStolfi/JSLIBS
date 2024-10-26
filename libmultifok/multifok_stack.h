/* Stack of images with limited depth of focus. */
/* Last edited on 2024-10-22 08:25:43 by stolfi */

#ifndef multifok_stack_H
#define multifok_stack_H

#define _GNU_SOURCE
#include <stdint.h>

#include <float_image.h>

#include <multifok_frame.h>
  
typedef struct multifok_stack_t 
  { int32_t NI;                /* Number of images in stack. */
    int32_t NC;                /* Number of channels in the scene images. */
    int32_t NX, NY;            /* Image dimensions of all images. */
    multifok_frame_t **frame;  /* Frames of the stack. */
  } multifok_stack_t;
  /* A stack of input frames {frame[ki]} for {ki} in {0..NI-1}. All frames must have the 
    same number of columns {NX} and of rows {NY}.  All scene images {frame[ki].sVal}
    must have the same number of channels {NC}. */
     
multifok_stack_t* multifok_stack_new(int32_t NI, int32_t NC, int32_t NX, int32_t NY);
  /* Creates a new stack record with the given paramenters.  The vector 
    {stack.frame} is allocated with space for {NI} frames, but is filled with {NULL}. */

multifok_stack_t* multifok_stack_read
  ( char *stackDir,
    bool_t gray,
    int32_t NI,
    double zFoc[],
    double zDep[],
    double hMin,
    double hMax
  );
  /* Reads a stack of simulated multifocus images from files
    "{stackDir}/{frameSub}/{name}.png", where {ki} ranges in {0..NI-1},
    {frameSub} is "frame-zf{FFF.FFFF}-df{DDD.DDDD}", {FFF.FFFF} and {DDD.DDDD} are
    respectively the values of {zFoc[ki]} and {zDep[ki]}, formatted as
    "%08.4f", and {name} is "sVal", "hAvg", "hDev", or "shrp".
    However, if {zDep[ki]} is {+INF}, the {frameSub} will be "frame-sharp" instead.
    
    Each frame is read with {multifok_frame_read}.  The
    parameters {gray}, {hMin} and {hMax} are passed to it.
    
    Returns those frames in a {multifok_stack_t} object {stack}. */

void multifok_stack_write
  ( multifok_stack_t *stack,
    char *stackDir,
    double hMin,
    double hMax
  );
  /* Writes the frames {stack.frame[0..NI-1]}, where {NI} is {stack->NI},
    to files "{stackDir}/{frameSub}/{name}.png", using {multifok_frame_write}.
    See that procedure for the file names and formats. */

#endif
