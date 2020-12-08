#ifndef boap_plate_H
#define boap_plate_H

/* A box-like plate with edges parallel to the axes. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <sign.h>
#include <frgb.h>
#include <r3.h>
#include <epswr.h>

/* All dimensions are in mm. */

typedef struct boap_plate_t
  { r3_t cmin;      /* World coordinates of nominal center of plate. */
    r3_t cmax;      /* World axes that correspond to local axes. */
    frgb_t *color;  /* Fill color for plotting. */
  } boap_plate_t;
  /* A record that defines a steel plate with edges parallel to the 
    world coordinate axes, of the given {color}. */

boap_plate_t *boap_plate_new(r3_t cmin, r3_t cmax, frgb_t *color);
  /* Allocates a new {boap_plate_t} record,  NOTE: the color record is 
    shared, not copied. */

/* BASIC DRAWING */

/* The following procedures assume that the observer is at {+oo}
    in the direction of axis {uax}. 
    
    If {pos} is {NAN}, draws a projection on a plane perpendicular to
    {uax}. If {pos} is not {NAN}, draws a cross-section by a plane
    perpendicular to {uax} at coordinate{pos}. In this latter case, only
    plates that intersect the plane are shown.
    
    On the projection, If {uax} is 0, the projection is {x=Y,y=Z}. If
    {uax} is 1, it is {x=X,y=Z}. If {uax} is 2, it is {x=X,y=Y}.
    
    The {fill} and {draw} parameters have the same meanning as in
    {epswr_rectangle}. */

void boap_plate_draw
  ( epswr_figure_t *epsf, 
    boap_plate_t *plt, 
    double penwd,
    bool_t fill, bool_t draw, 
    int8_t uax,
    double pos
  );
  /* Draws the plate {*plt} on {epsf} */

/* DRAWING WITH VISIBILITY */

typedef boap_plate_t *boap_plate_ref_t;

vec_typedef(boap_plate_ref_vec_t,boap_plate_ref_vec,boap_plate_ref_t);

void boap_plate_ref_vec_draw
  ( epswr_figure_t *epsf, 
    int32_t np,
    boap_plate_ref_vec_t *pltv, 
    bool_t fill, bool_t draw, 
    int8_t uax,
    double pos
  );
  /* Draws plates {pltv.e[0..np-1]} on {epsf}.  If {pos} is not NAN,
    sorts the plates according to the "painting order" visibility method. */

void boap_plate_print(FILE *wr, boap_plate_t *plt);
  /* Prints to {wr} a legible description of the plate {plt}. */

#endif
