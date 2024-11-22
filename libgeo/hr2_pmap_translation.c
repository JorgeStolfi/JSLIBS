/* See hr2_pmap_translation.h */
/* Last edited on 2024-11-09 14:03:19 by stolfi */ 

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

#include <hr2_pmap_translation.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

hr2_pmap_t hr2_pmap_translation_from_disp(r2_t *disp)
  { hr2_pmap_t M;
    r3x3_ident(&(M.dir));

    M.dir.c[0][1] = +disp->c[0];
    M.dir.c[0][2] = +disp->c[1];

    r3x3_ident(&(M.inv));
    M.inv.c[0][1] = -disp->c[0];
    M.inv.c[0][2] = -disp->c[1];
    return M;
  }
 
