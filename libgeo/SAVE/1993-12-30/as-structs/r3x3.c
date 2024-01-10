#include <math.h>
#include <stdio.h>
#include <ioprotos.h>
#include "r3.h"
#include "r3x3.h"

#define N 3

void r3x3_map_row (r3_t *a, r3x3_t *m, r3_t *res)
{
  int j, k;
  for (j=0; j<N; j++)
    { double s = 0.0;
      for (k=0; k<N; k++) s += a->c[k] * m->c[k,j];
      res->c[j] = s;
    }
}

void r3x3_map_col (r3x3_t *m, r3_t *a, r3_t *res)
{
  int i, k;
  for (i=0; i<N; i++)
    { double s = 0.0;
      for (k=0; k<N; k++) s += m->c[i,k] * a->c[k];
      res->c[i] = s;
    }
}

void r3x3_mul (r3x3_t *m, r3x3_t *n, r3x3_t *res)
{
  int i, j, k;
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      { double s = 0.0;
        for (k=0; k<N; k++)  s += m->c[i,k]*n->c[k,j];
        res->c[i,j] = s;
      }
}

double r3x3_det (r3x3_t *m)
{
  double d0 = m->c[1,1]*m->c[2,2] - m->c[1,2]*m->c[2,1];
  double d1 = m->c[1,0]*m->c[2,2] - m->c[1,2]*m->c[2,0];
  double d2 = m->c[1,0]*m->c[2,1] - m->c[1,1]*m->c[2,0];
  double det = d0 * m->c[0,0] - d1 * m->c[0,1] + d2 * m->c[0,2];
  return (det);
}

double r3x3_cof (r3x3_t *m, int ix, int jx)
{
  int i0, i1, j0, j1;
  double d;

  /* Select columns of minor:*/
  if      (ix == 0) {i0 = 1; i1 = 2; }
  else if (ix == 1) {i0 = 0; i1 = 2; }
  else if (ix == 2) {i0 = 0; i1 = 1; }
  if      (jx == 0) {j0 = 1; j1 = 2; }
  else if (jx == 1) {j0 = 0; j1 = 2; }
  else if (jx == 2) {j0 = 0; j1 = 1; }

  /* Compute determinant: */
  d = m->c[i0,j0]*m->c[i1,j1] - m->c[i0,j1]*m->c[i1,j0];
  if ((ix + jx) % 2) d = -d;
  return (d);
}

void r3x3_adj (r3x3_t *m, r3x3_t *res)
{
  int i, j;
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      res->c[i,j] = r3x3_cof (m, j, i);
}

void r3x3_print (FILE *f, r3x3_t *m)
{
  int i,j;
  for (i=0; i<N; i++)
    {
      for (j=0; j<N; j++) fprintf(f, "%16.8e", m->c[i,j]);
      fputc('\n', f);
    }
  fflush(f);
}

