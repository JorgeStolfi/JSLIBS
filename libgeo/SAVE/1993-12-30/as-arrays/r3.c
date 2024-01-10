#include <math.h>
#include <stdio.h>
#include <ioprotos.h>
#include "r3.h"

void r3_add (r3_t a, r3_t b, r3_t res)
{
  res[0] = a[0] + b[0];
  res[1] = a[1] + b[1];
  res[2] = a[2] + b[2];
}

void r3_sub (r3_t a, r3_t b, r3_t res)
{
  res[0] = a[0] - b[0];
  res[1] = a[1] - b[1];
  res[2] = a[2] - b[2];
}

void r3_scale (double s, r3_t a, r3_t res)
{
  res[0] = s * a[0];
  res[1] = s * a[1];
  res[2] = s * a[2];
}

double r3_dist (r3_t a, r3_t b)
{
  double d;
  r3_t o;
  o[0] = (a[0] - b[0]);
  o[1] = (a[1] - b[1]);
  o[2] = (a[2] - b[2]);
  d = sqrt(o[0]*o[0] + o[1]*o[1] + o[2]*o[2]);
  return (d);
}

double r3_orthize (r3_t a, r3_t u, r3_t res)
{
  double sau = a[0]*u[0] + a[1]*u[1] + a[2]*u[2];
  double suu = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
  double c = sau / suu;
  res[0] = a[0] - c * u[0];
  res[1] = a[1] - c * u[1];
  res[2] = a[2] - c * u[2];
  return (c);
}

double r3_normalize_inf (r3_t a)
{
  double d = 0.0;
  double a0 = fabs(a[0]);
  double a1 = fabs(a[1]);
  double a2 = fabs(a[2]);
  if (a0 > d) d = a0;
  if (a1 > d) d = a1;
  if (a2 > d) d = a2;
  a[0] /= d;
  a[1] /= d;
  a[2] /= d;
  return (d);
}

double r3_normalize (r3_t a)
{
  double d;

  d = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  a[0] /= d;
  a[1] /= d;
  a[2] /= d;
  return (d);
}

double r3_dot (r3_t a, r3_t b)
{
  double d;

  d = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return (d);
}

void r3_cross (r3_t a, r3_t b, r3_t res)
{
  res[0] = (a[1] * b[2]) - (a[2] * b[1]);
  res[1] = (a[2] * b[0]) - (a[0] * b[2]);
  res[2] = (a[0] * b[1]) - (a[1] * b[0]);
}

void r3_print (FILE *f, r3_t a)
{
  fprintf(f, 
    "(%16.8e %16.8e %16.8e)",
    a[0], a[1], a[2]
  );
}
    







