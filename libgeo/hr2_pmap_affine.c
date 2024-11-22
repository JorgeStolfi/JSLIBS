/* See hr2_pmap_affine.h */
/* Last edited on 2024-11-21 20:24:00 by stolfi */ 

#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <sign_get.h>
#include <sign.h>
#include <affirm.h>
#include <r2x2.h>
#include <r3x3.h>

#include <bool.h>
#include <r2.h>
#include <hr2.h>
#include <hr2_pmap.h>

#include <hr2_pmap_affine.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

hr2_pmap_t hr2_pmap_affine_from_mat_and_disp(r2x2_t *E, r2_t *d)
  {
    hr2_pmap_t M;

    M.dir.c[0][0] = 1.0;
    M.dir.c[0][1] = d->c[0];
    M.dir.c[0][2] = d->c[1];
    
    M.dir.c[1][0] = 0.0;
    M.dir.c[1][1] = E->c[0][0];
    M.dir.c[1][2] = E->c[0][1];
    
    M.dir.c[2][0] = 0.0;
    M.dir.c[2][1] = E->c[1][0];
    M.dir.c[2][2] = E->c[1][1];
    
    r2x2_t F;
    r2x2_inv(E, &(F));
    r2_t psid;
    r2x2_map_row(d, &(F), &psid);
    
    M.inv.c[0][0] = 1.0;
    M.inv.c[0][1] = - psid.c[0];
    M.inv.c[0][2] = - psid.c[1];
    
    M.inv.c[1][0] = 0.0;
    M.inv.c[1][1] = F.c[0][0];
    M.inv.c[1][2] = F.c[0][1];
    
    M.inv.c[2][0] = 0.0;
    M.inv.c[2][1] = F.c[1][0];
    M.inv.c[2][2] = F.c[1][1];
    
    return M;
  }

hr2_pmap_t hr2_pmap_affine_from_three_points(r2_t *o, r2_t *p, r2_t *q)
  {
    hr2_pmap_t M; 
    
    M.dir.c[0][0] = 1.0;
    M.dir.c[0][1] = o->c[0];
    M.dir.c[0][2] = o->c[1];

    M.dir.c[1][0] = 0.0;
    M.dir.c[1][1] = p->c[0] - o->c[0];
    M.dir.c[1][2] = p->c[1] - o->c[1];

    M.dir.c[2][0] = 0.0;
    M.dir.c[2][1] = q->c[0] - o->c[0];
    M.dir.c[2][2] = q->c[1] - o->c[1];

    r3x3_inv(&(M.dir), &(M.inv));

    /* Just in case: */
    M.inv.c[1][0] = 0.0;
    M.inv.c[2][0] = 0.0;
    
    return M;
  }

double hr2_pmap_affine_discr_sqr(hr2_pmap_t *M, hr2_pmap_t *N)
  {
    demand((M->dir.c[1][0] == 0) && (M->dir.c[2][0] == 0), "{M} is not affine");
    demand(M->dir.c[0][0] > 0, "map {M} does not preserve side");
    r3x3_t A; double wA = M->dir.c[0][0]; r3x3_scale(1/wA, &(M->dir), &A);
   
    demand((N->dir.c[1][0] == 0) && (N->dir.c[2][0] == 0), "{N} is not affine");
    demand(M->dir.c[0][0] > 0, "map {N} does not preserve side");
    r3x3_t B; double wB = N->dir.c[0][0]; r3x3_scale(1/wB, &(N->dir), &B);
   
    /* Hope the math is right: */
    r3x3_t E, H;
    r3x3_sub(&A, &B, &E);
    r3x3_mul_tr(&E, &E, &H);
    double h2 = (H.c[1][1] + H.c[2][2])/2;
    double d2 = H.c[0][0];
    return h2 + d2;
  }
