// See minu_gen.h
// Last edited on 2002-07-16 10:35:42 by stolfi

#include <math.h>
#include <minu_gen.h>

void ComputeDerivatives
  ( double u, double fu, double v, double fv, double x, double fx,
    double *Dfx, double *Dq, double *DDfx, double *DDq )
{
  double uv = v - u;
  double ux = x - u;
  double xv = v - x;
  double du = xv * (fx - fu);
  double dv = ux * (fv - fx);
  double q = uv * ux * xv;
  
  *Dfx = du*xv + dv*ux;
  *DDfx = 2.0 * (dv - du);
  if (q < 0.0) { *Dfx = -*Dfx; *DDfx = -*DDfx; q = -q; }
  *Dq = q;
  *DDq = q;
}

