/* See r4_extra.h */
/* Last edited on 2014-01-12 12:48:18 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <r4.h>
#include <r4_extra.h>

#include <affirm.h>

#define N 4

void r4_pick_ortho (r4_t *u, r4_t *r)
  {
    double t0 = u->c[0]; r->c[0] = -(u->c[1]); r->c[1] = t0; 
    double t2 = u->c[2]; r->c[2] = -(u->c[3]); r->c[3] = t2; 
  }
