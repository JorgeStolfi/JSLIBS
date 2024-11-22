/* See hr2_pmap_congruence.h */
/* Last edited on 2024-11-15 12:58:31 by stolfi */ 

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

#include <hr2_pmap_congruence.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

hr2_pmap_t hr2_pmap_congruence_from_point_and_dir(r2_t *p, r2_t *u, sign_t sgn)
  { 
    demand((sgn == -1) || (sgn == +1), "invalid {sgn}");
    /* Normalize {u} to unit length to get the image of the {(1,0)} vector: */
    r2_t du;  
    double mu = r2_dir(u, &du);
    demand(mu != 0, "invalid direction {u}");

    /* Get the unit vector {dv} that is to be the image of the {(0,1)} vector: */
    r2_t dv;
    if (sgn < 0)
      { dv = (r2_t){{ +du.c[1], -du.c[0] }}; }
    else
      { dv = (r2_t){{ -du.c[1], +du.c[0] }}; }
      
    /* Assemble the matrix: */
    hr2_pmap_t M; /* The resulting map. */
    M.dir.c[0][0] = 1.0;     
    M.dir.c[0][1] = p->c[0]; 
    M.dir.c[0][2] = p->c[1]; 

    M.dir.c[1][0] = 0.0;     
    M.dir.c[1][1] = du.c[0]; 
    M.dir.c[1][2] = du.c[1]; 

    M.dir.c[2][0] = 0.0;     
    M.dir.c[2][1] = dv.c[0]; 
    M.dir.c[2][2] = dv.c[1]; 
    
    r3x3_inv(&(M.dir), &(M.inv));

    return M;
  }
