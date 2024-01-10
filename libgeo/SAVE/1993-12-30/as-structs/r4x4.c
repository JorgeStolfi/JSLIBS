#include <math.h>
#include <stdio.h>
#include <ioprotos.h>
#include "r4x4.h"
#include "r4.h"
#include "r3x3.h"

#define N 4

void r4x4_map_row (r4_t *a, r4x4_t *m, r4_t *res)
{
  int j, k;
  for (j=0; j<N; j++)
    { double s = 0.0;
      for (k=0; k<N; k++) s += a->c[k] * m->c[k,j];
      res->c[j] = s;
    }
}

void r4x4_map_col (r4x4_t *m, r4_t *a, r4_t *res)
{
  int i, k;
  for (i=0; i<N; i++)
    { double s = 0.0;
      for (k=0; k<N; k++) s += m->c[i,k] * a->c[k];
      res->c[i] = s;
    }
}

void r4x4_mul (r4x4_t *m, r4x4_t *n, r4x4_t *res)
{
  int i, j, k;
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      { double s = 0.0;
        for (k=0; k<N; k++)  s += m->c[i,k]*n->c[k,j];
        res->c[i,j] = s;
      }
}

double r4x4_cof (r4x4_t *m, int ix, int jx)
{
  r3x3_t t;
  int i, j, ii, jj;
  double d;

  /* Extract minor: */
  ii = 0;
  for (i=0; i<N; i++)
    if (i != ix)
      { jj = 0;
        for (j=0; j<N; j++) 
          if (j != jx) { t.c[ii,jj] = m->c[i,j]; jj++; }
        ii++;
      }
  /* Compute its determinant: */
  d = r3x3_det(&t);
  if ((ix + jx) % 2) d = -d;
  return (d);
}

void r4x4_adj (r4x4_t *m, r4x4_t *res)
{
  int i, j;
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      res->c[i,j] = r4x4_cof (m, j, i);
}

void r4x4_print (FILE *f, r4x4_t *m)
{
  int i,j;
  for (i=0; i<N; i++)
    {
      for (j=0; j<N; j++) fprintf(f, "%16.8e", m->c[i,j]);
      fputc('\n', f);
    }
  fflush(f);
}

