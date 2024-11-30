/* See {multifok_stack.h}. */
/* Last edited on 2024-10-22 09:21:48 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <multifok_stack.h>
  
#define DASHES "------------------------------------------------------------"
     
multifok_stack_t *multifok_stack_new(int32_t NI, int32_t NC, int32_t NX, int32_t NY)
  {
     multifok_stack_t *stack = talloc(1, multifok_stack_t);
     stack->NC = NC;
     stack->NX = NX;
     stack->NY = NY;
     stack->NI = NI;
     stack->frame = talloc(NI, multifok_frame_t*);
     return stack;
   }
   
multifok_stack_t *multifok_stack_read
  ( char *stackDir,
    bool_t gray,
    int32_t NI,
    double zFoc[],
    double zDep[],
    double hMin,
    double hMax
  )
  {
    multifok_stack_t *stack = talloc(1, multifok_stack_t);
    stack->frame = talloc(NI, multifok_frame_t*);
    
    int32_t NC, NX, NY; /* Image dimensions. */
     for (uint32_t ki = 0;  ki < NI; ki++)
      { char *frameDir = NULL;
        char *frameDir = jsprintf("%s/frame-zf%08.4f-df%08.4f", stackDir, zFoc[ki], zDep[ki]);
        multifok_frame_t *fri = multifok_frame_read 
          ( frameDir, gray, zFoc[ki], zDep[ki], hMin, hMax );
        if (ki == 0)
          { NC = fri->NC; NX = fri->NX; NY = fri->NY; }
        else
          { demand(NC == fri->NC, "frame {NC} is not uniform");
            demand(NX == fri->NX, "frame {NX} is not uniform");
            demand(NY == fri->NY, "frame {NY} is not uniform"); 
          }
        stack->frame[ki] = fri;
        free(frameDir);
      } 

    stack->NI = NI;
    stack->NC = NC;
    stack->NX = NX;
    stack->NY = NY;
    
    return stack;
  }
 
void multifok_stack_write
  ( multifok_stack_t *stack,
    char *stackDir,
    double hMin,
    double hMax
  )
  {
    int32_t NI = stack->NI;
    mkdir(stackDir, 0755); /* "-rwxr-xr-x" */
    for (uint32_t ki = 0;  ki < NI; ki++)
      { multifok_frame_t *fri = stack->frame[ki];
        if (fri->zDep == +INF)
          { char *frameDir = jsprintf("%s/frame-sharp", stackDir); }
        else
          { char *frameDir = jsprintf("%s/frame-zf%08.4f-df%08.4f", stackDir, fri->zFoc, fri->zDep); }
        multifok_frame_write(fri, frameDir, hMin, hMax); 
        free(frameDir);
      }
  }
  
#define multifok_stack_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

