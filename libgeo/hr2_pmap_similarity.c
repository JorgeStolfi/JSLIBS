/* See hr2_pmap_similarity.h */
/* Last edited on 2024-11-15 12:58:19 by stolfi */ 

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

#include <hr2_pmap_similarity.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */


hr2_pmap_t hr2_pmap_similarity_from_two_points(r2_t *p, r2_t *q, sign_t sgn)
  {
    demand((sgn == -1) || (sgn == +1), "invalid {sgn}");

    /* Compute the vector {u} that is the image of the {(1,0)} vector: */
    r2_t u; r2_sub(q, p, &u);

    /* Get the squared length {d2} of {u}: */
    double d2 = r2_norm_sqr(&u);
    demand(d2 != 0, "points {p,q} coincide");

    /* Get the vector {v} that is to be the image of the {(0,1)} vector: */
    r2_t v;
    if (sgn < 0)
      { v = (r2_t){{ +u.c[1], -u.c[0] }}; }
    else
      { v = (r2_t){{ -u.c[1], +u.c[0] }}; }
      
    /* Assemble the matrix: */
    hr2_pmap_t M; /* The resulting map. */
    M.dir.c[0][0] = 1.0;        
    M.dir.c[0][1] = p->c[0];    
    M.dir.c[0][2] = p->c[1];    

    M.dir.c[1][0] = 0.0;         
    M.dir.c[1][1] = +u.c[0];    
    M.dir.c[1][2] = +u.c[1];    

    M.dir.c[2][0] = 0.0;         
    M.dir.c[2][1] = +v.c[0];    
    M.dir.c[2][2] = +v.c[1];    
    
    r3x3_inv(&(M.dir), &(M.inv));

    return M;
  }
