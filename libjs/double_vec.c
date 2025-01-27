/* See double_vec.h. */
/* Last edited on 2025-01-22 19:02:24 by stolfi */

/* Copyright © 2003 Jorge Stolfi, Unicamp. See note at end of file. */
/* Based on Params.m3 by J.Stolfi, DEC-SRC, 1988.  */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <float.h>

#include <affirm.h>
#include <vec.h>

#include <double_vec.h>

void double_vec_uniformize(double_vec_t *v, double defval)
  { uint32_t KC = v->ne;
    if (KC == 0) 
      { double_vec_trim(v, 1);
        v->e[0] = defval;
      }
    else if (KC > 1)
      { defval =  v->e[0];
        for (int32_t c = 1; c < KC; c++) 
          { demand(v->e[c] == defval, "inconsistent values"); }
        double_vec_trim(v, 1);
      }
  }

void double_vec_regularize(double_vec_t *v, int32_t NC, double defval)
  { uint32_t KC = v->ne;
    if (KC != NC) 
      { demand(KC <= 1, "range info specified for the wrong number of channels");
        /* Provide default value: */
        double val = (KC == 1 ? v->e[0] : defval);
        double_vec_trim(v, (vec_size_t)NC);
        for (int32_t c = 0; c < NC; c++) { v->e[c] = val; }
      }
  }
