/* colorfield_unif.h - uniform color field
** Last edited on 2004-11-06 02:59:51 by stolfi
**
** Copyright (C) 2003 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
*/

#ifndef colorfield_unif_H
#define colorfield_unif_H

#include <frgb.h>
#include <frgb_ops.h>
#include <colorfield.h>

/* USER ARGUMENTS FROM THE COMMAND LINE */

typedef struct cfld_unif_args_t  /* User specs for a unif field */
  { frgb_t color;  /* Raw color. */
  } cfld_unif_args_t;

cfld_unif_args_t *cfld_unif_parse(argparser_t *pp);
  /* Parses the next few command line arguments, starting at
    {argv[*argn]}, as a uniform color field specification --- i.e. 
    a single color value, in the {frgb_parse_color} format. */ 

cfld_unif_args_t *cfld_unif_make_args(frgb_t *color);
  /* Makes a default argument record, as if the user had specified 
    the given {color}. */
    
/* PREPROCESSED DATA */

typedef struct cfld_unif_params_t  /* Preprocessed params for unif field */
  { frgb_t color;       /* Adjusted color. */
  } cfld_unif_params_t;

cfld_unif_params_t *cfld_unif_compute_params
  ( cfld_unif_args_t *rfa, 
    frgb_adjuster_t *adjust,
    int logarithmic
  );
  /* Computes the gamma-corrected reference color {wp->color} from the
    user-specified one {rfa->color}, by applying the color adjustments embodied
    in {adjust} and {logarithmic}. */

void cfld_unif_eval
  ( cfld_unif_params_t *rfp,
    int logarithmic, 
    int col, 
    int row,
    frgb_t *fv,
    int chns
  );
  /* Evaluates the uniform field described by {rfp} at the pixel 
    in column {col} and row {row}. I.e., just returns {rfp->color}.
    The arguments {logarithmic}, {col}, and {row} are ignored. */
    
#endif
