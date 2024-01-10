#include <math.h>
#include <stdio.h>
#include <ioprotos.h>
#include "r4.h"
#include <js.h>

#define N 4

void r4_zero (r4_t *res)
{
  res->c[0] = 0.0;
  res->c[1] = 0.0;
  res->c[2] = 0.0;
  res->c[3] = 0.0;
}
  
void r4_axis (int i, r4_t *res)
{
  assert((i >= 0) && (i < N), "r4_axis: bad index");
  res->c[0] = 0.0;
  res->c[1] = 0.0;
  res->c[2] = 0.0;
  res->c[3] = 0.0;

  res->c[i] = 1.0;
}

void r4_add (r4_t *a, r4_t *b, r4_t *res)
{
  res->c[0] = a->c[0] + b->c[0];
  res->c[1] = a->c[1] + b->c[1];
  res->c[2] = a->c[2] + b->c[2];
  res->c[3] = a->c[3] + b->c[3];
}

void r4_sub (r4_t *a, r4_t *b, r4_t *res)
{
  res->c[0] = a->c[0] - b->c[0];
  res->c[1] = a->c[1] - b->c[1];
  res->c[2] = a->c[2] - b->c[2];
  res->c[3] = a->c[3] - b->c[3];
}

void r4_scale (double s, r4_t *a, r4_t *res)
{
  res->c[0] = s * a->c[0];
  res->c[1] = s * a->c[1];
  res->c[2] = s * a->c[2];
  res->c[3] = s * a->c[3];
}

void r4_mix_in (double s, r4_t *a, r4_t *res)
{
  res->c[0] += s * a->c[0];
  res->c[1] += s * a->c[1];
  res->c[2] += s * a->c[2];
  res->c[3] += s * a->c[3];
}

double r4_dist (r4_t *a, r4_t *b)
{
  double d;
  r4_t o;
  o.c[0] = (a->c[0] - b->c[0]);
  o.c[1] = (a->c[1] - b->c[1]);
  o.c[2] = (a->c[2] - b->c[2]);
  o.c[3] = (a->c[3] - b->c[3]);
  d = sqrt(o.c[0]*o.c[0] + o.c[1]*o.c[1] + o.c[2]*o.c[2] + o.c[3]*o.c[3]);
  return (d);
}

double r4_orthize (r4_t *a, r4_t *u, r4_t *res)
{
  double sau = a->c[0]*u->c[0] + a->c[1]*u->c[1] + a->c[2]*u->c[2] + a->c[3]*u->c[3];
  double suu = u->c[0]*u->c[0] + u->c[1]*u->c[1] + u->c[2]*u->c[2] + u->c[3]*u->c[3];
  double c = sau / suu;
  res->c[0] = a->c[0] - c * u->c[0];
  res->c[1] = a->c[1] - c * u->c[1];
  res->c[2] = a->c[2] - c * u->c[2];
  res->c[3] = a->c[3] - c * u->c[3];
  return (c);
}

double r4_normalize_inf (r4_t *a)
{
  double d = 0.0;
  double a0 = fabs(a->c[0]);
  double a1 = fabs(a->c[1]);
  double a2 = fabs(a->c[2]);
  double a3 = fabs(a->c[3]);
  if (a0 > d) d = a0;
  if (a1 > d) d = a1;
  if (a2 > d) d = a2;
  if (a3 > d) d = a3;
  a->c[0] /= d;
  a->c[1] /= d;
  a->c[2] /= d;
  a->c[3] /= d;
  return (d);
}

double r4_normalize (r4_t *a)
{
  double d;

  d = sqrt(a->c[0]*a->c[0] + a->c[1]*a->c[1] + a->c[2]*a->c[2] + a->c[3]*a->c[3]);
  a->c[0] /= d;
  a->c[1] /= d;
  a->c[2] /= d;
  a->c[3] /= d;
  return (d);
}

double r4_dot (r4_t *a, r4_t *b)
{
  double d;

  d = a->c[0]*b->c[0] + a->c[1]*b->c[1] + a->c[2]*b->c[2] + a->c[3]*b->c[3];
  return (d);
}

void r4_cross (r4_t *a, r4_t *b, r4_t *c, r4_t *res)
{
  double d01 = a->c[0]*b->c[1] - a->c[1]*b->c[0];
  double d02 = a->c[0]*b->c[2] - a->c[2]*b->c[0];
  double d12 = a->c[1]*b->c[2] - a->c[2]*b->c[1];
  double d03 = a->c[0]*b->c[3] - a->c[3]*b->c[0];
  double d13 = a->c[1]*b->c[3] - a->c[3]*b->c[1];
  double d23 = a->c[2]*b->c[3] - a->c[3]*b->c[2];

  res->c[0] = - d12*c->c[3] + d13*c->c[2] - d23*c->c[1];
  res->c[1] = + d02*c->c[3] - d03*c->c[2] + d23*c->c[0];
  res->c[2] = - d01*c->c[3] + d03*c->c[1] - d13*c->c[0];
  res->c[3] = + d01*c->c[2] - d02*c->c[1] + d12*c->c[0];
}

void r4_print (FILE *f, r4_t *a)
{
  fprintf(f, 
    "(%16.8e %16.8e %16.8e %16.8e)", 
    a->c[0], a->c[1], a->c[2], a->c[3]
  );
}
    
void r4_throw_cube (r4_t *res)
{
  int i;
  for (i=0; i<N; i++)
    { res->c[i] = 2.0 * frandom() - 1.0; }
}
  
void r4_throw_ball (r4_t *res)
{
  int i;
  double c, d;
  do
    { d = 0.0;
      for (i=0; i<N; i++)
        { c = 2.0 * frandom() - 1.0; 
          res->c[i] = c;
          d += c*c;
        }
    }
  while (d > 1.0);
}








