/* See {hr2_pmap_encode.h}. */
/* Last edited on 2024-09-17 15:38:36 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <hr2.h>
#include <r3x3.h>
#include <jsrandom.h>
#include <affirm.h>

#include <hr2_pmap_encode.h>
#include <hr2_pmap_throw_type.h>

hr2_pmap_t hr2_pmap_throw_type(hr2_pmap_type_t type, sign_t sgn)
  {
    bool_t debug = TRUE;
    if (debug) 
      { char *xtype = hr2_pmap_type_to_string(type);
        fprintf(stderr, "    > --- %s type = %s sgn = %+d ---\n", __FUNCTION__, xtype, sgn);
      }
    
    int32_t ny = hr2_pmap_encode_num_parameters(type);
    assert (ny <= 9);
    double y[ny];
    for (int32_t ky = 0; ky < ny; ky++) { y[ky] = 10*(drandom() - 0.5); }
    
    hr2_pmap_t M;
    hr2_pmap_decode(ny, y, type, sgn, &M);
    if (debug) {
      fprintf(stderr, "  thrown:\n");
      hr2_pmap_gen_print(stderr, &M, "%12.7f", "    ", "[ ","  "," ]\n    ", "[ "," "," ]", "\n");
    }
    demand(sgn*r3x3_det(&(M.dir)) > 0, "{hr2_pmap_decode} returns wrong handedness");
    
    /* Check that map has the requested type apart from {XY} swap: */
    hr2_pmap_t N = M;
    hr2_pmap_set_sign(&N, sgn);
    double eqtol = 1.0e-14;
    demand(hr2_pmap_is_type(&N, type, sgn, eqtol), "wrong map type or sign");

    r3x3_normalize(&(M.dir));
    r3x3_normalize(&(M.inv));
    if (debug) { fprintf(stderr, "    < --- %s ---\n", __FUNCTION__); }
    return M;
  }
