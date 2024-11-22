/* Last edited on 2024-11-18 09:22:25 by stolfi */
/* See {interp_spline_I.h}. */

/* !!! Implement the other cases !!! */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <interp_spline_I.h>

uint32_t interp_spline_I_compute_num_samples(int32_t ord)
  {
    demand(ord >= -1, "invalid {ord}");
    switch(ord)
      { 
        case -1:
          /* Nearest-sample replication: */
          return 1;
        case 0:
          /* Linear interpolation of two samples: */
          return 2;
        case 1:
          /* Quadratic C1 interpolation of 4 samples: */
          return 4;
        default:
          fprintf(stderr, "%s(%d) not implemented yet\n", __FUNCTION__, ord);
          return 0;
      }
  }
 
void interp_spline_I_get_weights(double z, int32_t ord, uint32_t nw, double wt[])
  {
    demand(ord >= -1, "invalid {ord}");
    bool_t debug = FALSE;
    if (ord == -1)
      { /* Just replicate the nearest sample: */
        assert(nw == 1);
        wt[0] = 1.0;
      }
    else
      { /* Polynomial interpolation of nearest {nw} samples: */
        /* Adjust {z} to be relative to first filter tap: */
        z = z - 0.5*(nw - 1);

        /* Get the raw index of the nearest lower tap: */
        int32_t iz = (int32_t)floor(z);

        /* Get the fraction {fz} and its complement {gz}: */
        double fz = z - iz;
        assert(fz >= 0.0);
        assert(fz < 1.0);
        double gz = 1.0 - fz;
        
        if (debug) { fprintf(stderr, "z = %9.6f  fz = %9.6f  gz = %9.6f", z, fz, gz); }

        switch(ord)
          { 
            case 0:
              /* Linear interpolation: */
              assert(nw == 2);
              wt[0] = gz;
              wt[1] = fz;
              break;
              
            case 1:
              assert(nw == 4);
              if (fz < 0.5)
                { wt[0] = ((+0.75)*fz - 0.5)*fz;
                  wt[1] = (-1.75)*fz*fz + 1;
                  wt[2] = ((+1.25)*fz + 0.5)*fz;
                  wt[3] = (-0.25)*fz*fz;
                }
              else
                { wt[0] = (-0.25)*gz*gz;
                  wt[1] = ((+1.25)*gz + 0.5)*gz;
                  wt[2] = (-1.75)*gz*gz + 1;
                  wt[3] = ((+0.75)*gz - 0.5)*gz;
                }
              break;
              
            default: 
              demand(FALSE, "not implemented yet");  
          }

        if (debug)
          { int32_t k;
            for (k = 0; k < nw; k++) { fprintf(stderr, "  wt[%d] = %10.7f\n", k, wt[k]); }
            fprintf(stderr, "\n");
          }
     }
  }

// # I-spline with 12 pieces of degree 3, hopefully of order 2
//   D0I2x0(z) = (((+0.003472)*z+0.000000)*z-0.000000)*z+0.000000
//   D0I2x1(z) = (((-0.024306)*z+0.010417)*z+0.010417)*z+0.003472
//   D0I2x2(z) = (((+0.031250)*z-0.062500)*z-0.041667)*z-0.000000
//   D0I2x3(z) = (((+0.114583)*z+0.031250)*z-0.072917)*z-0.072917
//   D0I2x4(z) = (((-0.138889)*z+0.375000)*z+0.333333)*z-0.000000
//   D0I2x5(z) = (((-0.194444)*z-0.041667)*z+0.666667)*z+0.569444
//   D0I2x6(z) = (((+0.194444)*z-0.625000)*z-0.000000)*z+1.000000
//   D0I2x7(z) = (((+0.138889)*z-0.041667)*z-0.666667)*z+0.569444
//   D0I2x8(z) = (((-0.114583)*z+0.375000)*z-0.333333)*z-0.000000
//   D0I2x9(z) = (((-0.031250)*z+0.031250)*z+0.072917)*z-0.072917
//   D0I2x10(z) = (((+0.024306)*z-0.062500)*z+0.041667)*z-0.000000
//   D0I2x11(z) = (((-0.003472)*z+0.010417)*z-0.010417)*z+0.003472
//   D0I2R(x) = ((x < -6.0)||(x > +6.0) ? 0 : (x < +0.0 ? (x < -3.0 ? (x < -4.0 ? (x < -5.0 ? D0I2x0(x+6.0) : D0I2x1(x+5.0)) : D0I2x2(x+4.0)) : (x < -1.0 ? (x < -2.0 ? D0I2x3(x+3.0) : D0I2x4(x+2.0)) : D0I2x5(x+1.0))) : (x < +3.0 ? (x < +2.0 ? (x < +1.0 ? D0I2x6(x-0.0) : D0I2x7(x-1.0)) : D0I2x8(x-2.0)) : (x < +5.0 ? (x < +4.0 ? D0I2x9(x-3.0) : D0I2x10(x-4.0)) : D0I2x11(x-5.0)))))
//   D0I2(x) = D0I2R(x/0.5)
//   I2(x) = D0I2(x)
