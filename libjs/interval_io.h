#ifndef interval_io_H
#define interval_io_H

/* Input/output for {interval_t} values. */
/* Last edited on 2012-01-04 00:31:46 by stolfi */ 

#define _GNU_SOURCE
#include <stdio.h>

#include <jsmath.h>
#include <interval.h>

void interval_print(FILE *wr, interval_t *X);
  /* Writes {X} to {wr}, with some default format. */

void interval_gen_print(FILE *wr, interval_t *X, char *fmt, char *lp, char *sep, char *rp);
  /* Writes {X} to {wr}, bracketed by {lp,rp} and with the two endpoints 
    separated by {sep}.  When NULL, they default to "%24.16e", "[", " ", and "]",
    respectively. */

#endif
