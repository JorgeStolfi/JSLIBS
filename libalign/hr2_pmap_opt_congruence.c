/* See {hr2_pmap_opt_translation.h}. */
/* Last edited on 2023-10-21 12:07:29 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>
#include <math.h>
 
#include <bool.h>
#include <r2.h>
#include <affirm.h>

#include <hr2_pmap_opt_congruence.h>

void hr2_pmap_opt_congruence_encode(hr2_pmap_t *M, double rad[], double y[])
  { 
    y[0] = M->dir.c[0][1] / rad[0];
    y[1] = M->dir.c[0][2] / rad[1];
    r2_t u = (r2_t){{ +M->dir.c[1][1], +M->dir.c[1][2] }};
    r2_t v = (r2_t){{ +M->dir.c[2][2], -M->dir.c[2][1] }};
    double angu = atan2(u.c[1], u.c[0]);
    double angv = atan2(v.c[1], v.c[0]);
    if (angv > angu + M_PI) { angv = angv - 2*M_PI; }
    if (angv <= angu - M_PI) { angv = angv + 2*M_PI; }
    double ang = (angu + angv)/2;
    y[2] = ang / rad[2];
  }

void hr2_pmap_opt_congruence_decode(double y[], double rad[], hr2_pmap_t *M)
  { M->dir.c[0][1] = y[0] * rad[0];
    M->dir.c[0][2] = y[1] * rad[1];
    double ang = y[2] * rad[2];
    double c = cos(ang), s = sin(ang);
    M->dir.c[1][1] = +c;
    M->dir.c[1][2] = +s;
    M->dir.c[2][1] = -s;
    M->dir.c[2][2] = +c;

    M->inv.c[0][1] = - c*M->dir.c[0][1] - s*M->dir.c[0][2];
    M->inv.c[0][2] = + s*M->dir.c[0][1] - c*M->dir.c[0][2];
    M->inv.c[1][1] = +c;
    M->inv.c[1][2] = -s;
    M->inv.c[2][1] = +s;
    M->inv.c[2][2] = +c;
  }
