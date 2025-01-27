/* double_vec.h -- tools for vectors of doubles. */
/* Last edited on 2025-01-22 19:01:09 by stolfi */

#ifndef double_vec_H
#define double_vec_H

/* Copyright © 2005 by the State University of Campinas (UNICAMP).*/
/* Last edited on 2006-04-01 10:31:38 by stolfi*/

#include <stdint.h>

#include <vec.h>
#include <affirm.h>

void double_vec_uniformize(double_vec_t *v, double defval);
  /* Make sure that the vector {v} has one element. If
    {v.ne == 0}, expands it to one element and sets it to {defval}.
    If {v.ne > 1}, requires that all elements be equal, and
    truncates it to one element. Otherwise does nothing. */

void double_vec_regularize(double_vec_t *v, int32_t NC, double defval);
  /* Make sure that the vector {v} has {NC} elements. If
    {v.ne == 0}, expands it to {NC} elements, and sets all elements
    to {defval}. If {v.ne == 1}, expands it to {NC} elements,
    replicating the first element. If {v.ne == NC},
    does nothing. Otherwise fails. */

#endif
