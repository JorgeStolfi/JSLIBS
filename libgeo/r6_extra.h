/* r6_extra.h --- additional operations on points and vectors of R^6 */
/* Last edited on 2014-01-12 15:11:26 by stolfilocal */

#ifndef r6_extra_H
#define r6_extra_H

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>

#include <r6.h>

void r6_pick_ortho (r6_t *u, r6_t *r);
  /* Sets {r} to an arbitrary vector orthogonal to {u},
    with the same L-infinity and Euclidean norms as {u}.
    In particular, if {u} is {(0,0,0,0,0,0)}, sets {r} to {(0,0,0,0,0,0)}.
    The result {r} is a continuous function of {u}.  */

#endif
