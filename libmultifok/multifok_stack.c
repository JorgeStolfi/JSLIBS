/* See {multifok_stack.h}. */
/* Last edited on 2025-02-08 11:09:37 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsprintf.h>

#include <multifok_stack.h>
  
#define DASHES "----------------------------------------------------------------------"
     
multifok_stack_t *multifok_stack_new(uint32_t NI, int32_t NC, int32_t NX, int32_t NY)
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
  ( char *stackFolder,
    bool_t gray,
    uint32_t NI,
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
      { char *frameFolder = jsprintf("%s/frame-zf%08.4f-df%08.4f", stackFolder, zFoc[ki], zDep[ki]);
        multifok_frame_t *fri = multifok_frame_read 
          ( frameFolder, gray, zFoc[ki], zDep[ki], hMin, hMax );
        if (ki == 0)
          { NC = fri->NC; NX = fri->NX; NY = fri->NY; }
        else
          { demand(NC == fri->NC, "frame {NC} is not uniform");
            demand(NX == fri->NX, "frame {NX} is not uniform");
            demand(NY == fri->NY, "frame {NY} is not uniform"); 
          }
        stack->frame[ki] = fri;
        free(frameFolder);
      } 

    stack->NI = NI;
    stack->NC = NC;
    stack->NX = NX;
    stack->NY = NY;
    
    return stack;
  }
 
void multifok_stack_write
  ( multifok_stack_t *stack,
    char *stackFolder,
    double hMin,
    double hMax
  )
  {
    uint32_t NI = stack->NI;
    mkdir(stackFolder, 0755); /* "-rwxr-xr-x" */
    for (uint32_t ki = 0;  ki < NI; ki++)
      { multifok_frame_t *fri = stack->frame[ki];
        char *frameFolder = NULL;
        if (fri->zDep == +INF)
          { frameFolder = jsprintf("%s/frame-sharp", stackFolder); }
        else
          { frameFolder = jsprintf("%s/frame-zf%08.4f-df%08.4f", stackFolder, fri->zFoc, fri->zDep); }
        multifok_frame_write(fri, frameFolder, hMin, hMax); 
        free(frameFolder);
      }
  }
  
#define multifok_stack_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

