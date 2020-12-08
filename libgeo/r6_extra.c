/* See r6_extra.h */
/* Last edited on 2014-01-12 15:12:13 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <r6.h>
#include <r6_extra.h>

#include <affirm.h>

#define N 6

void r6_pick_ortho (r6_t *u, r6_t *r)
  {
    double t0 = u->c[0]; r->c[0] = -(u->c[1]); r->c[1] = t0; 
    double t2 = u->c[2]; r->c[2] = -(u->c[3]); r->c[3] = t2; 
    double t4 = u->c[4]; r->c[4] = -(u->c[5]); r->c[5] = t4; 
  }
