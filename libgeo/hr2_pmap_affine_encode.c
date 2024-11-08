/* See {hr2_pmap_affine_encode.h}. */
/* Last edited on 2024-11-03 12:25:02 by stolfi */

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
  { 
    r3x3_t *A = &(M->dir);

    bool_t paranoia = TRUE;
    if (paranoia)
      { /* Paranoia checks: */
        demand(hr2_pmap_is_affine(M, 1.0e-14), "map is not affine");
        demand(r3x3_det(A) > 0, "map is is not positive");
      }
 
    double A00 = A->c[0][0];
    demand(A00 >= 1.0e-200, "map is not affine");
   
    double A11 = A->c[1][1] / A00;
    double A12 = A->c[1][2] / A00;
    double A21 = A->c[2][1] / A00;
    double A22 = A->c[2][2] / A00;
    
    r2_t u = (r2_t){{ A11,  A12 }};
    r2_t v = (r2_t){{ A21,  A22 }};
    
    double u2 = r2_dot(&u, &u);
    demand(u2 >= 1.0e-200, "map is singular");
    
    double vc = r2_dot(&u, &v)/u2;
    double vs = fabs(r2_det(&u, &v)/u2); /* Ignore handedness. */
    demand(vs > 0, "map is singular");
    
    y[0] = vc;
    y[1] = 0.5*(vs - 1/vs);
    
    double f = sqrt(u2);
    y[2] = 0.5*(f - 1/f);
    y[3] = atan2(u.c[1], u.c[0]);
    /* Hack to handle the DAMN minus zero: */
    if (y[3] == -M_PI) { y[3] = +M_PI; }
  
    y[4] = A->c[0][1] / A00;
    y[5] = A->c[0][2] / A00;
  }

void hr2_pmap_affine_decode(double y[], hr2_pmap_t *M)
  { 
    r3x3_t *A = &(M->dir);
    r3x3_ident(A);
    
    /* Shearing: */
    demand(isfinite(y[0]), "invalid {X} shearing parameter {y[0]}");
    double vc = y[0];
    demand(isfinite(y[1]), "invalid {Y} pre scaling parameter {y[1]}");
    double vs = y[1]+ hypot(y[1], 1);

    /* Scaling: */
    demand(isfinite(y[2]), "invalid scaling parameter {y[2]}");
    double f = y[2] + hypot(y[2], 1);
    demand(isfinite(f), "scaling parameter {y[2]} too large");
    
    /* Rotation: */
    demand(isfinite(y[3]),  "invalid rotation angle {y[3]}");
    double ang = y[3];
    double fuc = f*cos(ang);
    double fus = f*sin(ang);
    
    A->c[1][1] = + fuc;
    A->c[1][2] = + fus;

    A->c[2][1] = vc*fuc - vs*fus;
    A->c[2][2] = vc*fus + vs*fuc;

    A->c[0][1] = y[4];
    A->c[0][2] = y[5];

    r3x3_inv(&(M->dir), &(M->inv));

    /* Just in case: */
    M->inv.c[1][0] = 0.0;
    M->inv.c[2][0] = 0.0;
  }
