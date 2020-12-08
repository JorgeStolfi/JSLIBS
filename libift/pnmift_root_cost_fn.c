/* pnmift_root_cost_fn.c - implementation of pnmift_root_cost_fn.h */
/* Last edited on 2010-06-06 16:54:32 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <affirm.h>
#include <bool.h>
#include <frgb.h>
#include <frgb_ops.h>
  
#include <ift.h>
#include <pnmift_root_cost_fn.h>

pnmift_root_cost_fn_t *pnmift_root_cost_fn_from_name(char *name)
  {
    if (strcmp(name, "zero") == 0)
      { return &pnmift_root_cost_fn_zero; }
    else if (strcmp(name, "lum") == 0)
      { return &pnmift_root_cost_fn_lum; }
    else 
      { demand(FALSE, "bad function name"); }
    return NULL;
  }

pnmift_root_cost_t pnmift_root_cost_fn_zero(frgb_t q, int chns)
  {
    return 0;
  }

pnmift_root_cost_t pnmift_root_cost_fn_lum(frgb_t q, int chns)
  {
    return frgb_Y(&q);
  }
