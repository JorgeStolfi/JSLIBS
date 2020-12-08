#ifndef jsdithers_H
#define jsdithers_H

/* dithers.h - ordered dither matrices. */
/* Last edited on 2017-06-21 00:58:27 by stolfilocal */

#include <stdint.h>
#include <bool.h>

uint8_t *get_dither_matrix(int n, bool_t cluster);
  /* Returns a pointer to a dither matrix with {n} rows and {n} columns,
     stored linearized by rows.
     
     If {cluster} is false, the matrix contains all the integers in
     {0..n^2-1}, each appearing only once, in some optimized random-like
     order. Allowed values of {n} are 8 and 16.
     
     If {cluster} is true, the values are all integers in {0..n^2/2-1},
     each appearing exactly twice, roughly sorted by distance
     from element {[0,0]}.  Allowed values of {n} are 6, 8, and 16. */

#endif
