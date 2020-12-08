/* colorfield_ramp.h - definitions for ramp fields
** Last edited on 2004-11-06 02:59:42 by stolfi
**
** Copyright (C) 2003 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
*/

#ifndef colorfield_ramp_H
#define colorfield_ramp_H

#include <frgb.h>
#include <frgb_ops.h>
#include <colorfield.h>

/* USER ARGUMENTS FROM THE COMMAND LINE */

typedef struct cfld_ramp_args_t  /* User specs for a ramp field */
  { cfld_int_pair_t p[3];      /* Indices of the reference pixels. */
    frgb_t color[3];  /* Raw reference colors. */
  } cfld_ramp_args_t;

cfld_ramp_args_t *cfld_ramp_parse_general(argparser_t *pp);
  /* Parses the next few command line arguments, starting at
    {argv[*argn]}, as a color ramp specification. Assumes that the
    "-field ramp" arguments have already been parsed. Increments {*argn}.
    
    The field arguments have the form {H0 V0 COLOR0 H1 V1 COLOR1 H2 V2
    COLOR2}, specifying the coordinates of three non-collinear pixels,
    and the respective are color values (in the {frgb_parse_color}
    format). */ 

/* PREPROCESSED DATA */

typedef struct cfld_ramp_params_t  /* Preprocessed params for ramp field */
  { frgb_t forg;        /* Color at point {(0,0)}. */
    frgb_t fcol;        /* Coefficients of pixel column. */
    frgb_t frow;        /* Coefficients of pixel row. */
  } cfld_ramp_params_t;

cfld_ramp_params_t *cfld_ramp_compute_params
  ( cfld_ramp_args_t *rfa, 
    frgb_adjuster_t *adjust,
    int logarithmic
  );
  /* Computes the gamma-corrected color ramp parameters from the 
    user arguments {rfa}, by applying the color adjustments 
    embodied in {adjust} and {logarithmic}. */

/* EVALUATION */

void cfld_ramp_eval
  ( cfld_ramp_params_t *rfp,
    int logarithmic, 
    int col, 
    int row,
    frgb_t *fv,
    int chns
  );
  /* Evaluates the ramp field described by {rfp} at the pixel 
    in column {col} and row {row}. */

#endif
