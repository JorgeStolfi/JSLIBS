#ifndef boap_model_H
#define boap_model_H

/* Geometric models for one item of the Boaretto 113 side porch. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <sign.h>
#include <frgb.h>
#include <epswr.h>

#include <boap_plate.h>
#include <boap_dim.h>
#include <boap_bar.h>

typedef struct boap_model_t 
  { /* Bounding box: */
    double cmin[3]; 
    double cmax[3];
    /* Bars: */
    int32_t nb; /* Number of bars in model */
    boap_bar_ref_vec_t barv;
    /* Plates: */
    int32_t np; /* Number of plates in model */
    boap_plate_ref_vec_t pltv;
    /* Dimension lines: */
    int32_t nd; /* Number of dimension lines in model */
    boap_dim_ref_vec_t dimv;
    
  } boap_model_t;
  /* A set of 3D geometric elements.  Only plates for now. */

boap_model_t *boap_model_new(void);
  /* Creates a new model, initially empty. */

void boap_model_add_plate(boap_model_t *mod, boap_plate_t *plt);
  /* Adds a new plate to the model {mod}. */
  
void boap_model_add_bar (boap_model_t *mod, boap_bar_t *bar);
  /* Adds the {bar} to the  model {mod}, by smashing it into plates. */

void boap_model_add_dim(boap_model_t *mod, boap_dim_t *dim);
  /* Adds a new dimension spec to the model {mod}. */

void boap_model_draw
  ( epswr_figure_t *epsf, 
    boap_model_t *mod, 
    int8_t uax,
    double pos
  );
  /* Draws the model {mod} onto figure {epsf}, as seen from {+oo}
    in the direction of axis {uax}.   If {pos} is not {NAN},
    draws a cross-section at coordinate {pos} on that axis. */

#endif
