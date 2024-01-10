#include <math.h>
#include <stdio.h>
#include <ioprotos.h>
#include "r3.h"
#include <js.h>

#define N 3

void r3_zero (r3_t *res)
{
  res->c[0] = 0.0;
  res->c[1] = 0.0;
  res->c[2] = 0.0;
}
  
void r3_axis (int i, r3_t *res)
{
  assert((i >= 0) && (i < N), "r3_axis: bad index");
  res->c[0] = 0.0;
  res->c[1] = 0.0;
  res->c[2] = 0.0;

  res->c[i] = 1.0;
}

void r3_add (r3_t *a, r3_t *b, r3_t *res)
{
  res->c[0] = a->c[0] + b->c[0];
  res->c[1] = a->c[1] + b->c[1];
  res->c[2] = a->c[2] + b->c[2];
}

void r3_sub (r3_t *a, r3_t *b, r3_t *res)
{
  res->c[0] = a->c[0] - b->c[0];
  res->c[1] = a->c[1] - b->c[1];
  res->c[2] = a->c[2] - b->c[2];
}

void r3_scale (double s, r3_t *a, r3_t *res)
{
  res->c[0] = s * a->c[0];
  res->c[1] = s * a->c[1];
  res->c[2] = s * a->c[2];
}

void r3_mix_in (double s, r3_t *a, r3_t *res)
{
  res->c[0] += s * a->c[0];
  res->c[1] += s * a->c[1];
  res->c[2] += s * a->c[2];
}

double r3_dist (r3_t *a, r3_t *b)
{
  double d;
  r3_t o;
  o.c[0] = (a->c[0] - b->c[0]);
  o.c[1] = (a->c[1] - b->c[1]);
  o.c[2] = (a->c[2] - b->c[2]);
  d = sqrt(o.c[0]*o.c[0] + o.c[1]*o.c[1] + o.c[2]*o.c[2]);
  return (d);
}

double r3_orthize (r3_t *a, r3_t *u, r3_t *res)
{
  double sau = a->c[0]*u->c[0] + a->c[1]*u->c[1] + a->c[2]*u->c[2];
  double suu = u->c[0]*u->c[0] + u->c[1]*u->c[1] + u->c[2]*u->c[2];
  double c = sau / suu;
  res->c[0] = a->c[0] - c * u->c[0];
  res->c[1] = a->c[1] - c * u->c[1];
  res->c[2] = a->c[2] - c * u->c[2];
  return (c);
}

double r3_normalize_inf (r3_t *a)
{
  double d = 0.0;
  double a0 = fabs(a->c[0]);
  double a1 = fabs(a->c[1]);
  double a2 = fabs(a->c[2]);
  if (a0 > d) d = a0;
  if (a1 > d) d = a1;
  if (a2 > d) d = a2;
  a->c[0] /= d;
  a->c[1] /= d;
  a->c[2] /= d;
  return (d);
}

double r3_normalize (r3_t *a)
{
  double d;

  d = sqrt(a->c[0]*a->c[0] + a->c[1]*a->c[1] + a->c[2]*a->c[2]);
  a->c[0] /= d;
  a->c[1] /= d;
  a->c[2] /= d;
  return (d);
}

double r3_dot (r3_t *a, r3_t *b)
{
  double d;

  d = a->c[0]*b->c[0] + a->c[1]*b->c[1] + a->c[2]*b->c[2];
  return (d);
}

void r3_cross (r3_t *a, r3_t *b, r3_t *res)
{
  res->c[0] = (a->c[1] * b->c[2]) - (a->c[2] * b->c[1]);
  res->c[1] = (a->c[2] * b->c[0]) - (a->c[0] * b->c[2]);
  res->c[2] = (a->c[0] * b->c[1]) - (a->c[1] * b->c[0]);
}

void r3_print (FILE *f, r3_t *a)
{
  fprintf(f, 
    "(%16.8e %16.8e %16.8e)",
    a->c[0], a->c[1], a->c[2]
  );
}
    
void r3_throw_cube (r3_t *res)
{
  int i;
  for (i=0; i<N; i++)
    { res->c[i] = 2.0 * frandom() - 1.0; }
}
  
void r3_throw_ball (r3_t *res)
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








