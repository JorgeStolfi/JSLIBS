/* See {hr2_pmap_encode.h}. */
/* Last edited on 2024-11-20 12:02:38 by stolfi */

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
#include <hr2_pmap_throw_by_type.h>

hr2_pmap_t hr2_pmap_throw_by_type(hr2_pmap_type_t type, sign_t sgn)
  {
    demand((sgn == -1) || (sgn == +1), "invalid sgn");
    bool_t debug = FALSE;
    if (debug) 
      { char *xtype = hr2_pmap_type_to_string(type);
        fprintf(stderr, "      > --- %s type = %s sgn = %+d ---\n", __FUNCTION__, xtype, sgn);
      }
    hr2_pmap_t M;
    if (type == hr2_pmap_type_GENERIC)
      { hr2_pmap_throw(&M);
        hr2_pmap_set_sign(&M, sgn);
      }
    else
      { uint32_t ny = hr2_pmap_encode_num_parameters(type);
        assert (ny <= 9);
        double y[ny];
        for (int32_t ky = 0; ky < ny; ky++) { y[ky] = 10*(drandom() - 0.5); }

        hr2_pmap_decode(ny, y, type, sgn, &M);
        if (debug)
            { fprintf(stderr, "        thrown:\n");
              hr2_pmap_gen_print(stderr, &M, "%12.7f", "          ", "[ ","  "," ]\n          ", "[ "," "," ]", "\n");
            }
        double eqtol = 1.0e-14;
        demand(hr2_pmap_is_type(&M, type, sgn, eqtol), "wrong map type or sign");
      }
    /* The chances of the matrix being all zeros are nil: */
    r3x3_normalize(&(M.dir));
    r3x3_normalize(&(M.inv));
    if (debug) { fprintf(stderr, "      < --- %s ---\n", __FUNCTION__); }
    return M;
  }
