/* See hermite3.h */
/* Last edited on 2014-07-27 16:58:14 by stolfilocal */

#define hermite3_C_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"
  
#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>

#include <hermite3.h>  

double hermite3_interp(double v0, double d0, double v1, double d1, double a, double b, double t)
  { /* Map to the range {[0 _ 1]}: */
    double h = b - a;
    double r = t - a;
    double s = b - t;
    if (h != 1.0) { r /= h; s /= h; d0 *= h; d1 *= h; }
    double r2 = r*r;
    double s2 = s*s;
    double P0 = (2*r + 1)*s2;
    double P1 = (2*s + 1)*r2;
    return v0*P0 + (d0*s - d1*r)*r*s + v1*P1;
  }

void hermite3_subsample(int nx, double x[], double dx[], int ns, int ny, double y[])
  {
    demand(nx > 0, "invalid {ns}");
    demand(ns > 0, "invalid {ns}");
    demand(ny == ns*(nx-1) + 1, "inconsistent lengths");
    bool_t smooth = (nx >= 3); /* True means should use C1 cubic interpolation. */
    bool_t self = smooth && (dx == NULL); /* {dx} is self-allocated. */
    if (self)
      { /* Allocate and estimate the derivatives: */
        dx = notnull(malloc(nx*sizeof(double)), "no mem");
        hermite3_estimate_derivs(nx, x, dx);
      }
    /* Copy the first sample: */
    y[0] = x[0];
    if (ns == 1)
      { /* Just copy the rest: */
        int ix;
        for (ix = 1; ix < nx; ix++) { y[ix] = x[ix]; }
      }
    else 
      { /* Actual interpolation: */
        /* Initialize values {u2,u1,v,v1} and (if smooth) derivatives {d0,d1}: */
        double v1 = x[0]; /* Value at end of previous interval. */
        double d1 = (smooth ? dx[0] : NAN); /* Derivative at end of previous interval. */
        int iy = 0; /* Last new sample set. */
        int ix;
        for (ix = 1; ix < nx; ix++)
          { /* Get the values {v0,v1} at {ix-1,ix}: */
            double v0 = v1;
            v1 = x[ix];
            /* Get the derivatives {d0,d1} at {ix-1,ix}, if needed, else {NAN,NAN}: */
            double d0 = d1;  
            d1 = (smooth ? dx[ix] : NAN);
            /* Interpolate new samples: */ 
            int jn;
            for (jn = 1; jn < ns; jn++)
              { /* Interpolation at fraction {jn/ns}: */
                double r0 = ((double)jn)/((double)ns);
                double r1 = 1 - r0;
                iy++;
                y[iy] = (smooth ? hermite3_interp(v0, d0, v1, d1, 0.0, 1.0, r0) : r1*v0 + r0*v1);
              }
            /* Copy the non-interpolated point: */
            iy++;
            y[iy] = x[ix];
          }
      }
    /* Cleanup: */
    if (self) { free(dx); }
  }
  
void hermite3_estimate_derivs(int nx, double x[], double dx[])
  {
    demand(nx >= 0, "invalid {nx}");
    if (nx == 0)
      { /* Nothing to do: */ }
    else if (nx == 1)
      { /* Assume constant: */
        dx[0] = 0.0;
      }
    else if (nx == 2)
      { /* Straight line: */
        double d = x[1] - x[0];
        dx[0] = d; dx[1] = d;
      }
    else if (nx == 3)
      { /* Parabola: */
        dx[0] = hermite3_estimate_deriv_0_1_2(x[0], x[1], x[2]);
        dx[1] = hermite3_estimate_deriv_1_0_1(x[0], x[2]);
        dx[2] = - hermite3_estimate_deriv_0_1_2(x[2], x[1], x[0]);
      }
    else
      { /* Cubics ({nx >= 4}): */
        dx[0] = hermite3_estimate_deriv_0_1_3(x[0], x[1], x[2], x[3]);
        dx[1] = hermite3_estimate_deriv_1_1_2(x[0], x[1], x[2], x[3]);
        int i;
        for (i = 2; i < nx-2; i++) { dx[i] = hermite3_estimate_deriv_2_0_2(x[i-2], x[i-1], x[i+1], x[i+2]); }
        dx[nx-2] = - hermite3_estimate_deriv_1_1_2(x[nx-1], x[nx-2], x[nx-3], x[nx-4]);
        dx[nx-1] = - hermite3_estimate_deriv_0_1_3(x[nx-1], x[nx-2], x[nx-3], x[nx-4]);
      }
  }
  
double hermite3_estimate_deriv_0_1_2(double v, double vp, double vpp)
  { return 0.5*(-3*v + 4*vp - vpp); }

double hermite3_estimate_deriv_1_0_1(double vm, double vp)
  { return 0.5*(vp - vm); }

double hermite3_estimate_deriv_2_0_2(double vmm, double vm, double vp, double vpp)
  { return 2*(vp - vm)/3 - (vpp - vmm)/12; }

double hermite3_estimate_deriv_1_1_2(double vm, double v, double vp, double vpp)
  { return (-2*vm -3*v + 6*vp - vpp)/6; }

double hermite3_estimate_deriv_0_1_3(double v, double vp, double vpp, double vppp)
  { return (-11*v + 18*vp - 9*vpp + 2*vppp)/6; }
