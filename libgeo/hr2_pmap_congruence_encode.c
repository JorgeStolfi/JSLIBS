/* See {hr2_pmap_translation_encode.h}. */
/* Last edited on 2024-09-17 12:20:51 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>
#include <math.h>
 
#include <bool.h>
#include <sign.h>
#include <r2.h>
#include <affirm.h>

#include <hr2_pmap_congruence_encode.h>

void hr2_pmap_congruence_encode(hr2_pmap_t *M, double y[])
  { 
    double A00 = M->dir.c[0][0];
    demand(A00 > 1.0e-200, "map is not a congruence");
    
    y[0] = M->dir.c[0][1] / A00;
    y[1] = M->dir.c[0][2] / A00;
    
    double A11 = M->dir.c[1][1] / A00;
    double A12 = M->dir.c[1][2] / A00;
    double A21 = M->dir.c[2][1] / A00;
    double A22 = M->dir.c[2][2] / A00;

    if (A11*A22 - A21*A12 > 0)
      { y[2] = atan2(A12, A11); }
    else
      { y[2] = atan2(A22, A21); }

    /* Paranoia check: */
    double R1 = hypot(A11, A12);
    demand(fabs(R1 - 1) <= 1.0e-14, "map is not an isometry");
    double R2 = hypot(A21, A22);
    demand(fabs(R2 - 1) <= 1.0e-14, "map is not an isometry");
  }

void hr2_pmap_congruence_decode(double y[], hr2_pmap_t *M)
  { 
    r3x3_ident(&M->dir);
    M->dir.c[0][1] = y[0];
    M->dir.c[0][2] = y[1];
    
    double c = cos(y[2]), s = sin(y[2]);
    
    M->dir.c[1][1] = + c;
    M->dir.c[1][2] = + s;

    M->dir.c[2][1] = - s;
    M->dir.c[2][2] = + c;
      
    r3x3_inv(&(M->dir), &(M->inv));

    /* Just in case: */
    M->inv.c[1][0] = 0.0;
    M->inv.c[2][0] = 0.0;
  }
