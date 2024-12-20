/* colorfield.h - procedurally-defined images
** Last edited on 2008-05-25 03:22:55 by stolfi
**
** Copyright (C) 2003 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
*/

#ifndef colorfield_H
#define colorfield_H

#include <stdint.h>
#include <limits.h>
#include <math.h>

#include <bool.h>
#include <i2.h>

#include <frgb.h>
#include <frgb_ops.h>
#include <argparser.h>

/* GENERIC FIELDS */

typedef enum { cfld_UNIF, cfld_RAMP, cfld_WAVY } cfld_kind_t;

typedef struct cfld_params_t /* Preprocessed data for the color field */
  { cfld_kind_t kind;   /* Which kind of field? */
    void *params;       /* Specific field parameters; type depends on {kind}. */
    bool_t logarithmic;    /* TRUE to compute the field in log scale. */
  } cfld_params_t;
     
void cfld_eval(cfld_params_t *fld, int32_t col, int32_t row, frgb_t *fv, int32_t chns);
  /* Evaluates the first {chns} components of 
    the color field {fld} at the given pixel. */

/* COMMAND LINE ARGUMENTS */

/* Pixel indices {(column,row)}, from top left: */
typedef i2_t cfld_int_pair_t;

typedef struct cfld_args_t /* User specs for a the color field */
  { cfld_kind_t kind;   /* Which kind of field? */
    void *args;         /* Specific raw ccommand line args; type depends on {kind}. */
  } cfld_args_t;

cfld_args_t *cfld_parse(argparser_t *pp);
  /* Parses the next few command line arguments, starting at
    {argv[*argn]}, as a color field specification (what
    follows the "-field" keyword). */ 
    
cfld_args_t *cfld_make_args_uniform(frgb_t *color);
  /* Creates a description of an uniform field of the given color. */ 

cfld_params_t *cfld_compute_params
  ( cfld_args_t *fargs,
    frgb_adjuster_t *adjust,
    bool_t logarithmic
  );
  /* Computes the field parameters from the user specs {fargs}. 
    Any color parameter {v} specified by {fargs} as applying to a
    point {p} gets adjusted by calling {adjust(&v,p)}. */

#endif
