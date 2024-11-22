 /* Last edited on 2024-11-18 09:23:49 by stolfi */
/* See {interp_spline_O.h}. */

/* !!! Implement the other cases !!! */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <interp_spline_O.h>

uint32_t interp_spline_O_compute_num_samples(int32_t ord)
  {
    demand(ord >= -1, "invalid {ord}");
    fprintf(stderr, "%s(%d) not implemented yet\n", __FUNCTION__, ord);
    return 0;
    // /* !!! Rethink !!! */
    // if (ord == -1)
    //   { return 1; }
    // else if (ord == 0)
    //   { return 2; }
    // else if (ord == 1)
    //   { return 4; }
    // else
    //   { demand(FALSE, "invalid interpolation ord"); }
  }
 
void interp_spline_O_get_weights(double z, int32_t ord, uint32_t nw, double w[])
  {
    bool_t debug = FALSE;

    demand(ord >= -1, "invalid {ord}");

    /* Adjust {z} to be relative to start of first pixel: */
    z = z - 0.5*(nw - 1);
    
    /* Get the raw index of the first pixel: */
    int32_t iz = (int32_t)floor(z);
    
    /* Get the fraction {fz}: */
    double fz = z - iz;
    assert(fz >= 0.0);
    assert(fz < 1.0);
    
    /* Compute the fraction's complement {gz}: */
    double gz = 1.0 - fz;

    demand(FALSE, "not implemented yet");
    
    if (debug)
      { fprintf(stderr, "z = %9.6f  fz = %9.6f  gz = %9.6f", z, fz, gz);
        for (int32_t k = 0; k < nw; k++) { fprintf(stderr, "  w[%d] = %10.7f\n", k, w[k]); }
        fprintf(stderr, "\n");
      }
  }

// # O-spline with 4 pieces of degree 1, hopefully of order 0
//   D0O0x0(z) = (-0.000000)*z+0.000000
//   D0O0x1(z) = (+1.000000)*z+0.000000
//   D0O0x2(z) = (-1.000000)*z+1.000000
//   D0O0x3(z) = (+0.000000)*z+0.000000
//   D0O0R(x) = ((x < -2.0)||(x > +2.0) ? 0 : (x < +0.0 ? (x < -1.0 ? D0O0x0(x+2.0) : D0O0x1(x+1.0)) : (x < +1.0 ? D0O0x2(x-0.0) : D0O0x3(x-1.0))))
//   D0O0(x) = D0O0R(x)
//   O0(x) = D0O0(x)
// # O-spline with 6 pieces of degree 3, hopefully of order 1
//   D0O1x0(z) = (((-0.083333)*z+0.083333)*z-0.000000)*z+0.000000
//   D0O1x1(z) = (((+0.583333)*z-0.500000)*z-0.083333)*z+0.000000
//   D0O1x2(z) = (((-1.333333)*z+1.666667)*z+0.666667)*z-0.000000
//   D0O1x3(z) = (((+1.333333)*z-2.333333)*z+0.000000)*z+1.000000
//   D0O1x4(z) = (((-0.583333)*z+1.250000)*z-0.666667)*z-0.000000
//   D0O1x5(z) = (((+0.083333)*z-0.166667)*z+0.083333)*z+0.000000
//   D0O1R(x) = ((x < -3.0)||(x > +3.0) ? 0 : (x < +0.0 ? (x < -1.0 ? (x < -2.0 ? D0O1x0(x+3.0) : D0O1x1(x+2.0)) : D0O1x2(x+1.0)) : (x < +2.0 ? (x < +1.0 ? D0O1x3(x-0.0) : D0O1x4(x-1.0)) : D0O1x5(x-2.0))))
//   D0O1(x) = D0O1R(x)
//   O1(x) = D0O1(x)
// # O-spline with 8 pieces of degree 4, hopefully of order 2
//   D0O2x0(z) = ((((+0.017361)*z-0.017361)*z-0.000000)*z+0.000000)*z+0.000000
//   D0O2x1(z) = ((((-0.065972)*z-0.003472)*z+0.052083)*z+0.017361)*z-0.000000
//   D0O2x2(z) = ((((+0.093750)*z+0.413194)*z-0.354167)*z-0.152778)*z-0.000000
//   D0O2x3(z) = ((((-0.045139)*z-1.156250)*z+1.447917)*z+0.753472)*z+0.000000
//   D0O2x4(z) = ((((-0.045139)*z+1.336806)*z-2.291667)*z+0.000000)*z+1.000000
//   D0O2x5(z) = ((((+0.093750)*z-0.788194)*z+1.447917)*z-0.753472)*z+0.000000
//   D0O2x6(z) = ((((-0.065972)*z+0.267361)*z-0.354167)*z+0.152778)*z-0.000000
//   D0O2x7(z) = ((((+0.017361)*z-0.052083)*z+0.052083)*z-0.017361)*z-0.000000
//   D0O2R(x) = ((x < -4.0)||(x > +4.0) ? 0 : (x < +0.0 ? (x < -2.0 ? (x < -3.0 ? D0O2x0(x+4.0) : D0O2x1(x+3.0)) : (x < -1.0 ? D0O2x2(x+2.0) : D0O2x3(x+1.0))) : (x < +2.0 ? (x < +1.0 ? D0O2x4(x-0.0) : D0O2x5(x-1.0)) : (x < +3.0 ? D0O2x6(x-2.0) : D0O2x7(x-3.0)))))
//   D0O2(x) = D0O2R(x)
//   O2(x) = D0O2(x)
