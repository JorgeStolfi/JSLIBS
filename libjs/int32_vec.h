/* int32_vec.h -- toold for vectors of {int32_t}. */
/* Last edited on 2025-01-22 19:05:35 by stolfi */

#ifndef int32_vec_H
#define int32_vec_H

/* Copyright © 2005 by the State University of Campinas (UNICAMP).*/
/* Last edited on 2006-04-01 10:31:38 by stolfi*/

#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <affirm.h>

void int32_vec_regularize(int32_vec_t *v, int32_t NC, int32_t defval);
  /* Make sure that the color vector {v} has {NC} color channels. If
    {v.ne == 0}, expands it to {NC} channels, and sets all elements
    to {defval}. If {v.ne == 1}, expands it to {NC} channels,
    replicating the first element. Otherwise does nothing. */

#endif
