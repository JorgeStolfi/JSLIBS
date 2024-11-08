/* See {hr2_pmap_similarity_encode.h}. */
/* Last edited on 2024-11-03 12:12:05 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>
#include <math.h>
 
#include <bool.h>
#include <sign.h>
#include <r2.h>
#include <affirm.h>

#include <hr2_pmap_similarity_encode.h>

void hr2_pmap_similarity_encode(hr2_pmap_t *M, double y[])
  { 
    r3x3_t *A = &(M->dir);
    
    bool_t paranoia = TRUE;
    if (paranoia)
      { /* Paranoia checks: */
        demand(hr2_pmap_is_similarity(M, NAN, 1.0e-14), "map is not a similarity");
        demand(r3x3_det(A) > 0, "map is is not a positive similarity");
      }
  
    double A00 = A->c[0][0];
    demand(A00 > 1.0e-200, "map is not affine");
   
    double A11 = A->c[1][1] / A00;
    double A12 = A->c[1][2] / A00;
    double A21 = A->c[2][1] / A00;
    double A22 = A->c[2][2] / A00;
      
    /* Scaling: */
    double f = sqrt(hypot(A11, A12)*hypot(A21, A22));
    demand(f > 1.0e-200, "scale factor of map is zero");
    y[0] = 0.5*(f - 1/f);

    /* Rotation: */
    y[1] = atan2(A12, A11);
    /* Hack to handle the DAMN minus zero: */
    if (y[1] == -M_PI) { y[1] = +M_PI; }
    
    y[2] = A->c[0][1] / A00;
    y[3] = A->c[0][2] / A00;
  }

void hr2_pmap_similarity_decode(double y[], hr2_pmap_t *M)
  { 
    r3x3_t *A = &(M->dir);
    A->c[0][0] = 1.00;
    A->c[1][0] = 0.00;
    A->c[2][0] = 0.0;
    
    /* Scaling: */
    demand(isfinite(y[0]), "invalid scaling parameter {y[0]}");
    double f = hypot(y[0], 1) + y[0];
    demand(isfinite(f), "scaling parameter {y[0]} too large");
    
    /* Rotation: */
    demand(isfinite(y[1]), "invalid rotation angle {y[1]}");
    double fc = f*cos(y[1]), fs = f*sin(y[1]);
    
    A->c[1][1] = + fc;
    A->c[1][2] = + fs;

    A->c[2][1] = - fs;
    A->c[2][2] = + fc;
      
    /* Translation: */
    demand(isfinite(y[2]), "invalid X translation {y[2]}");
    demand(isfinite(y[3]), "invalid Y translation {y[3]}");
    A->c[0][1] = y[2];
    A->c[0][2] = y[3];
    
    r3x3_inv(A, &(M->inv));

    /* Just in case: */
    M->inv.c[1][0] = 0.0;
    M->inv.c[2][0] = 0.0;
  }
