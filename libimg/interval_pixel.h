#ifndef interval_pixel_H
#define interval_pixel_H

/* Interval-valued image pixels. */
/* Last edited on 2007-07-11 22:25:45 by stolfi */ 

#include <r2.h>
#include <r2x2.h>
#include <bool.h>
#include <interval.h>

void ivpix_accum_pixel(int chns, interval_t vs[], double wt, interval_t v[], double *wtotP);
  /* Adds {wt*vs[k]} to {v[k]} for {k} in {0..chns-1}. Also adds {wt} to {*wtotP}. */

void ivpix_debug_itv_pixel(char *label, double x, double y, int chns, interval_t v[], char *tail);
  /* Prints the point {(x,y)} and the interval-valued samples {v[0..chns-1]} to {stderr}. */

void ivpix_make_pixel_undef(int chns, interval_t v[]);
  /* Sets {v[0..chns-1]} to the undefined interval {[-INF _ +INF]}. */

void ivpix_scale_pixel(int chns, double s, interval_t v[]);
  /* Scales {v[0..chns-1]} by {s}, assumed positive. */

float ivpix_floatize_interval(interval_t *v);
  /* Picks the midpoint of {v}, safely. */

void ivpix_print_interval(FILE *wr, interval_t *v, int width, int prec);
  /* Prints the interval {v} to {wr}, with format "%{width}.{prec}f". */

void ivpix_print_bound(FILE *wr, double v, int width, int prec);
  /* Prints the interval bound {v} to {wr}, with format "%{width}.{prec}f". */

#endif
