/* pnmift_path_cost_fn.c - implementation of pnmift_path_cost_fn.h */
/* Last edited on 2010-06-06 14:51:37 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <affirm.h>
#include <bool.h>
  
#include <ift.h>
#include <pnmift_arc_cost_fn.h>
#include <pnmift_path_cost_fn.h>

/* FUNCTION SELECTION BY NAME */

pnmift_path_cost_fn_t *pnmift_path_cost_fn_from_name(char *name)
  {
    if (strcmp(name, "max") == 0)
      { return &pnmift_path_cost_fn_max; }
    else if (strcmp(name, "sum") == 0)
      { return &pnmift_path_cost_fn_sum; }
    else if (strcmp(name, "mono") == 0)
      { return &pnmift_path_cost_fn_mono; }
    else 
      { demand(FALSE, "bad function name"); }
    return NULL;
  }

/* PATH COST FUNCTIONS */

pnmift_path_cost_t pnmift_path_cost_fn_max(pnmift_path_cost_t sC, pnmift_arc_cost_t aC)
  {
    return fmax(sC, aC);
  }

pnmift_path_cost_t pnmift_path_cost_fn_sum(pnmift_path_cost_t sC, pnmift_arc_cost_t aC)
  {
    return sC + aC;
  }

pnmift_path_cost_t pnmift_path_cost_fn_mono(pnmift_path_cost_t sC, pnmift_arc_cost_t aC)
  {
    if (aC < sC)
      { return +INFINITY; }
    else
      { return aC; }
  }
