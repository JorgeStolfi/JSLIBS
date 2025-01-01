/* Miscellaneous general-purpose definitions */
/* Last edited on 2024-12-31 16:42:18 by stolfi */

#ifndef cpk_basic_H
#define cpk_basic_H

#include <stdint.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>

#include <bool.h>
#include <vec.h>

/* HANDY BASIC TYPES */

#define INF INFINITY
  /* Plus infinity as a {double} value. */

typedef uint32_t nat;

typedef struct ui2_t { uint32_t c[2]; } ui2_t;
  /* A pair of unsigned integers, e. g. vertex indices. */

vec_typedef(ui2_vec_t, ui2_vec, ui2_t);

/* COORDINATES OF POINTS */

/* Get X and Y coordinates from an {r2_t}, {i2+t}, or {ui2.t}: */
#define X(p) ((p).c[0])
#define Y(p) ((p).c[1])

/* TIME/DATE UNTILITIES */

char *today(void);
  /* Returns a newly allocated string containing 
    today's date in the format "yy-mm-dd hh:mm:ss" */

double cpk_cpu_time_1(void);
double cpk_cpu_time_2(void);
  /* Both procedures return the current user time of the process, in
    microseconds. They differ in implementation and resolution:
    {cpk_now_1} is based on the older {times} routine {cpk_now_2} on
    the POSIX {clock_gettime} routine. */

/* MATH UTILITIES */

int32_t intcmp(int32_t x, int32_t y);
int32_t dblcmp(double x, double y);
  /* Returns -1, 0, or +1 depending on whether {x} is less than,
    equal to, or greater than {y}. */

#endif
