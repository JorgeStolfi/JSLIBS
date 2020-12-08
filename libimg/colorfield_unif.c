/* see colorfield_unif.h 
** Last edited on 2008-05-25 03:23:06 by stolfi
**
** Copyright (C) 2003 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
*/

#include <colorfield_unif.h>

#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include <argparser.h>

#include <frgb.h>
#include <frgb_ops.h>
#include <colorfield.h>
  
cfld_unif_args_t *cfld_unif_parse(argparser_t *pp)
  {
    cfld_unif_args_t *rfa = (cfld_unif_args_t *)malloc(sizeof(cfld_unif_args_t));
    frgb_t color = frgb_parse_color(pp);
    rfa->color = color;
    return rfa;
  }

cfld_unif_params_t *cfld_unif_compute_params
  ( cfld_unif_args_t *rfa, 
    frgb_adjuster_t *adjust,
    int logarithmic
  )
  { cfld_unif_params_t *rfp = (cfld_unif_params_t *)malloc(sizeof(cfld_unif_params_t));
    /* Note: user-input colors are all RGB, even when {o->gray} is true. */
    frgb_t fv = adjust(&(rfa->color), 0, 0);
    /* {logarithmic} is irrelevant for uniform fields. */
    rfp->color = fv;
    return rfp;
  }

void cfld_unif_eval
  ( cfld_unif_params_t *rfp,
    int logarithmic, 
    int col, 
    int row,
    frgb_t *fv,
    int chns
  )
  { (*fv) = rfp->color; }
  
cfld_unif_args_t *cfld_unif_make_args(frgb_t *color)
  { cfld_unif_args_t *args = malloc(sizeof(cfld_unif_args_t *));
    args->color = *color;
    return args;
  }
