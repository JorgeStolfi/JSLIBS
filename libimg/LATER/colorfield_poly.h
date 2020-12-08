/* polyfield - a polynomial color field
** Last edited on 2007-12-27 14:49:50 by stolfi
**
** Copyright (C) 2003 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
*/

#ifndef colorfield_poly_H
#define colorfield_poly_H

#include <frgb.h>
#include <frgb_ops.h>
#include <colorfield.h>

/* TO BE COMPLETED

  !!! Should use a tensor basis instead of the 
  product of two unidimensional splines. */
    
typedef struct cfld_poly_args_t  /* Color field specifications from command line. */
  { cfld_int_pair_t deg;      /* Degree of spline along each axis. */
    cfld_int_pair_t nknots;   /* Number of knots along each axis. */
    int npt;                  /* Number of sample points. */
    cfld_int_pair_t *pt;      /* Sample points. */
    frgb_t *color;            /* Corresponding colors. */
  } cfld_poly_args_t;
  /* A {cfld_poly_args_t} record describes a smoothly varying color field.
    The field is defined by a set of sample points {pt[i]}
    and the corresponding color values {color[i]}.  
    
    The color field is the product of two polynomial splines 
    of the specified degrees. */

/* PREPROCESSED DATA */

typedef struct cfld_poly_params_t
  { int deg[2];    /* Degree along each axis. */
    frgb_t *w;     /* Bézier coefficients. */
    int additive;  /* TRUE to add, FALSE to multiply. */
  } cfld_poly_params_t;
  /* A {cfld_poly_params_t} record describes a variable color field {f}.
    Each component of {f(p)} at a pixel {p} at column {col}
    and row {row} is a polynomial on the coordinates {h,v},
    where {h = col/cols} and {v = row/rows}.
    
    If {additive = TRUE}, a field of degree {d} is determined by an
    array of {{d+1}^2} coefficients {w[i,j]} where {i} and {j} range
    in {0..d}. The field value is {f(p) = wfun(p)} where
      { wfun(p) = \sum w[i,j] B^{deg[0]}_{i}(h) B^{deg[1]}_{j}(v) }
    where {B} is the Bernstein-Bézier polynomial
      { B^d_{k}(x) = \choose{d}{k} (1 - x)^k x^{d-k} }
    Note that the field is always positive if all coefficients are positive.
    
    If {additive = FALSE}, the field value is {f(p) = exp(wfun(p))},
    where {wfun(p)} is defined as above.
    
    Note also that the {row} index is 0 at the *top* row of the image,
    and increases *downwards*.*/ 

cfld_poly_params_t *cfld_poly_compute_bezier_coeffs
  ( cfld_poly_args_t *w, 
    frgb_t *orgColor, 
    frgb_t *botColor
  );
  /* Computes the Bézier coefficients from the user-given field
    specs {w}, already corrected for gamma and other applicable
    factors. */

frgb_t cfld_poly_eval_bezier(cfld_poly_params_t wtb, double h, double v);
  /* Evaluates the Bézier color field {wtb} at the pixel 
    with relativ coordinates {h,v}. */

#endif
