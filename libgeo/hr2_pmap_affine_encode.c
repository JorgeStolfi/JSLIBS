/* See {hr2_pmap_affine_encode.h}. */
/* Last edited on 2024-09-17 16:30:22 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>
#include <math.h>
 
#include <bool.h>
#include <r2.h>
#include <hr2.h>
#include <hr2_pmap.h>
#include <affirm.h>

#include <hr2_pmap_affine_encode.h>

void hr2_pmap_affine_encode(hr2_pmap_t *M, double y[])
  { r3x3_t *A = &(M->dir);
    double A00 = M->dir.c[0][0];
    demand(A00 >= 1.0e-200, "map is not affine");
  
    y[0] = A->c[0][1] / A00;
    y[1] = A->c[0][2] / A00;
    
    r2_t u = (r2_t){{ (+ A->c[1][1] + A->c[2][1])/A00,  (+ A->c[1][2] + A->c[2][2]) / A00 }};
    r2_t v = (r2_t){{ (- A->c[1][1] + A->c[2][1])/A00,  (- A->c[1][2] + A->c[2][2]) / A00 }};
    
    double u2 = r2_dot(&u, &u);
    demand(u2 >= 1.0e-200, "map is singular");
    double um = sqrt(u2);
    
    y[2] = atan2(u.c[1], u.c[0]);
    y[3] = um - 1/um;
    
    double va = r2_dot(&u, &v)/u2;
    double vb = r2_det(&u, &v)/u2;
    if (vb < 0) { va = -va; vb = -vb; }
    y[4] = va;
    y[5] = vb - 1/vb;
  }

void hr2_pmap_affine_decode(double y[], hr2_pmap_t *M)
  { r3x3_t *A = &(M->dir);
    r3x3_ident(A);
    A->c[0][1] = y[0];
    A->c[0][2] = y[1];
    
    double ang = y[2];
    double R = y[3]/2 + hypot(y[3]/2, 1);
    r2_t u = (r2_t){{ R*cos(ang), R*sin(ang) }};  /* Image of vector {(+1,+1)}. */
    r2_t w = (r2_t){{ - u.c[1], u.c[0] }};        /* {u} rotated 90 deg CCW. */
    
    double va = y[4];
    double vb = y[5]/2 + hypot(y[5]/2, 1);
    
    r2_t v; r2_mix(va, &u, vb, &w, &v);  /* Image of vector {(-1,+1)}. */
    
    A->c[1][1] = 0.5*(u.c[0] - v.c[0]);
    A->c[1][2] = 0.5*(u.c[1] - v.c[1]);

    A->c[2][1] = 0.5*(u.c[0] + v.c[0]);
    A->c[2][2] = 0.5*(u.c[1] + v.c[1]);
    
    r3x3_inv(&(M->dir), &(M->inv));

    /* Just in case: */
    M->inv.c[1][0] = 0.0;
    M->inv.c[2][0] = 0.0;
  }
