/* See {hr2_pmap_translation_encode.h}. */
/* Last edited on 2024-11-03 06:41:15 by stolfi */

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
    r3x3_t *A = &(M->dir);

    bool_t paranoia = TRUE;
     if (paranoia)
      { /* Paranoia checks: */
        demand(hr2_pmap_is_congruence(M, 1.0e-14), "map is not a congruence");
        demand(r3x3_det(A) > 0, "map is is not a positive congruence");
      }
   
    double A00 = A->c[0][0];
    demand(A00 > 1.0e-200, "map is not affine");
    
    /* Rotation: */
    double A11 = A->c[1][1] / A00;
    double A12 = A->c[1][2] / A00;
    y[0] = atan2(A12, A11);
    /* Hack to handle the DAMN minus zero: */
    if (y[0] == -M_PI) { y[0] = +M_PI; }
    
    y[1] = A->c[0][1] / A00;
    y[2] = A->c[0][2] / A00;
  }

void hr2_pmap_congruence_decode(double y[], hr2_pmap_t *M)
  { 
    r3x3_t *A = &(M->dir);
    r3x3_ident(&M->dir);
    
    /* Rotation: */
    demand(isfinite(y[0]), "invalid rotation angle {y[0]}");
    double c = cos(y[0]), s = sin(y[0]);
    A->c[1][1] = + c;
    A->c[1][2] = + s;
    A->c[2][1] = - s;
    A->c[2][2] = + c;

    /* Translation: */
    demand(isfinite(y[1]), "invalid X translation {y[1]}");
    demand(isfinite(y[2]), "invalid Y translation {y[2]}");
    A->c[0][1] = y[1];
    A->c[0][2] = y[2];
      
    r3x3_inv(&(M->dir), &(M->inv));

    /* Just in case: */
    M->inv.c[1][0] = 0.0;
    M->inv.c[2][0] = 0.0;
  }
