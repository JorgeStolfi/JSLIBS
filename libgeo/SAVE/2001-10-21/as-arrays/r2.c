/* See r2.h */
/* Last edited on 2001-10-21 21:39:52 by stolfi */

#include <math.h>
#include <stdio.h>
#include <ioprotos.h>
#include "r2.h"
#include <js.h>

#define N 2

void r2_zero (r2_t res)
{
  res[0] = 0.0;
  res[1] = 0.0;
}
  
void r2_axis (int i, r2_t res)
{
  assert((i >= 0) && (i < N), "r2_axis: bad index");
  res[0] = 0.0;
  res[1] = 0.0;

  res[i] = 1.0;
}

void r2_add (r2_t a, r2_t b, r2_t res)
{
  res[0] = a[0] + b[0];
  res[1] = a[1] + b[1];
}

void r2_sub (r2_t a, r2_t b, r2_t res)
{
  res[0] = a[0] - b[0];
  res[1] = a[1] - b[1];
}

void r2_scale (double s, r2_t a, r2_t res)
{
  res[0] = s * a[0];
  res[1] = s * a[1];
}

void r2_mix_in (double s, r2_t a, r2_t res)
{
  res[0] += s * a[0];
  res[1] += s * a[1];
}

double r2_dist (r2_t a, r2_t b)
{
  double d0, d1;
  d0 = (a[0] - b[0]);
  d1 = (a[1] - b[1]);
  return sqrt(d0*d0 + d1*d1);
}

double r2_orthize (r2_t a, r2_t u, r2_t res)
{
  double sau = a[0]*u[0] + a[1]*u[1];
  double suu = u[0]*u[0] + u[1]*u[1];
  double c = sau / suu;
  res[0] = a[0] - c * u[0];
  res[1] = a[1] - c * u[1];
  return (c);
}

double r2_normalize_inf (r2_t a)
{
  double d = 0.0;
  double a0 = fabs(a[0]);
  double a1 = fabs(a[1]);
  if (a0 > d) d = a0;
  if (a1 > d) d = a1;
  a[0] /= d;
  a[1] /= d;
  return (d);
}

double r2_normalize (r2_t a)
{
  double d = sqrt(a[0]*a[0] + a[1]*a[1]);
  a[0] /= d;
  a[1] /= d;
  return (d);
}

double r2_dot (r2_t a, r2_t b)
{
  return a[0]*b[0] + a[1]*b[1];
}

double r2_cross (r2_t a, r2_t b)
{
  return a[0]*b[1] - a[1]*b[0];
}

void r2_print (FILE *f, r2_t a)
{
  fprintf(f, "(%16.8e %16.8e)", a[0], a[1]);
}

void r2_throw_cube (r2_t res)
{
  res[0] = 2.0 * frandom() - 1.0;
  res[1] = 2.0 * frandom() - 1.0;
}
  
void r2_throw_ball (r2_t res)
{
  double x, y;
  do
    { x = 2.0 * frandom() - 1.0; res[0] = x;
      y = 2.0 * frandom() - 1.0; res[1] = y;
    }
  while (x*x + y*y > 1.0);
}
 
