/* r4_extra.h --- additional operations on points and vectors of R^4 */
/* Last edited on 2014-01-12 15:10:56 by stolfilocal */

#ifndef r4_extra_H
#define r4_extra_H

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>

#include <r4.h>

void r4_pick_ortho (r4_t *u, r4_t *r);
  /* Sets {r} to an arbitrary vector orthogonal to {u},
    with the same L-infinity and Euclidean norms as {u}.
    In particular, if {u} is {(0,0,0,0)}, sets {r} to {(0,0,0,0)}.
    The result {r} is a continuous function of {u}.  */

#endif
