#include <math.h>
#include <stdio.h>
#include <ioprotos.h>
#include "r4.h"

void r4_add (r4_t a, r4_t b, r4_t res)
{
  res[0] = a[0] + b[0];
  res[1] = a[1] + b[1];
  res[2] = a[2] + b[2];
  res[3] = a[3] + b[3];
}

void r4_sub (r4_t a, r4_t b, r4_t res)
{
  res[0] = a[0] - b[0];
  res[1] = a[1] - b[1];
  res[2] = a[2] - b[2];
  res[3] = a[3] - b[3];
}

void r4_scale (double s, r4_t a, r4_t res)
{
  res[0] = s * a[0];
  res[1] = s * a[1];
  res[2] = s * a[2];
  res[3] = s * a[3];
}

double r4_dist (r4_t a, r4_t b)
{
  double d;
  r4_t o;
  o[0] = (a[0] - b[0]);
  o[1] = (a[1] - b[1]);
  o[2] = (a[2] - b[2]);
  o[3] = (a[3] - b[3]);
  d = sqrt(o[0]*o[0] + o[1]*o[1] + o[2]*o[2] + o[3]*o[3]);
  return (d);
}

double r4_orthize (r4_t a, r4_t u, r4_t res)
{
  double sau = a[0]*u[0] + a[1]*u[1] + a[2]*u[2] + a[3]*u[3];
  double suu = u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
  double c = sau / suu;
  res[0] = a[0] - c * u[0];
  res[1] = a[1] - c * u[1];
  res[2] = a[2] - c * u[2];
  res[3] = a[3] - c * u[3];
  return (c);
}

double r4_normalize_inf (r4_t a)
{
  double d = 0.0;
  double a0 = fabs(a[0]);
  double a1 = fabs(a[1]);
  double a2 = fabs(a[2]);
  double a3 = fabs(a[3]);
  if (a0 > d) d = a0;
  if (a1 > d) d = a1;
  if (a2 > d) d = a2;
  if (a3 > d) d = a3;
  a[0] /= d;
  a[1] /= d;
  a[2] /= d;
  a[3] /= d;
  return (d);
}

double r4_normalize (r4_t a)
{
  double d;

  d = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3]);
  a[0] /= d;
  a[1] /= d;
  a[2] /= d;
  a[3] /= d;
  return (d);
}

double r4_dot (r4_t a, r4_t b)
{
  double d;

  d = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
  return (d);
}

void r4_cross (r4_t a, r4_t b, r4_t c, r4_t res)
{
  double d01 = a[0]*b[1] - a[1]*b[0];
  double d02 = a[0]*b[2] - a[2]*b[0];
  double d12 = a[1]*b[2] - a[2]*b[1];
  double d03 = a[0]*b[3] - a[3]*b[0];
  double d13 = a[1]*b[3] - a[3]*b[1];
  double d23 = a[2]*b[3] - a[3]*b[2];

  res[0] = - d12*c[3] + d13*c[2] - d23*c[1];
  res[1] = + d02*c[3] - d03*c[2] + d23*c[0];
  res[2] = - d01*c[3] + d03*c[1] - d13*c[0];
  res[3] = + d01*c[2] - d02*c[1] + d12*c[0];
}

void r4_print (FILE *f, r4_t a)
{
  fprintf(f, 
    "(%16.8e %16.8e %16.8e %16.8e)", 
    a[0], a[1], a[2], a[3]
  );
}
    







