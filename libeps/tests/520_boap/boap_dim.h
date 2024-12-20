/* Dimension display elements for the Boaretto 113 side gate drawings. */
/* Last edited on 2024-12-05 10:16:37 by stolfi */

#ifndef boap_dim_H
#define boap_dim_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <sign.h>
#include <vec.h>
#include <r3.h>
#include <frgb.h>
#include <epswr.h>

#include <boap_plate.h>

typedef struct boap_dim_t
  { r3_t a;                   /* The {a} reference point (world coords). */
    r3_t b;                   /* The {b} reference point (world coords). */
    double agap;              /* Distance for start of extension segment for {a} (world mm). */
    double bgap;              /* Distance for start of extension segment for {b} (world mm). */
    double elen;              /* Min length of extension segments. */
    bool_t inner;             /* {TRUE} to show dimension inside. */
    int8_t decimals;          /* Number of decimal fractions to show. */
  } boap_dim_t;
  /* A dimension display element. */
  
boap_dim_t *boap_dim_new
  ( double Xa, double Ya, double Za,
    double Xb, double Yb, double Zb,
    double agap,
    double bgap, 
    double elen,
    bool_t inner,
    int8_t decimals
  );
  /* Creates a new {boap_dim_t} record. */

void boap_dim_draw(epswr_figure_t *epsf,  boap_dim_t *dim, int8_t uax);
  /* Draws the dimension {dim} on {epsf} assuming that the observer
    is at infinity in the direction of axis {uax}. 
    Omits the dimension if it is directed parallel to that axis. */
 
typedef boap_dim_t *boap_dim_ref_t;

vec_typedef(boap_dim_ref_vec_t, boap_dim_ref_vec, boap_dim_ref_t);

void boap_dim_print(FILE *wr, boap_dim_t *dim);
  /* Prints to {wr} a legible description of the given {dim}. */

#endif

