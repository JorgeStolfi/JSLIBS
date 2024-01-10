/* See r3x3.h */
/* Last edited on 2001-09-30 21:57:40 by stolfi */

#include <math.h>
#include <stdio.h>
#include <ioprotos.h>
#include "r3x3.h"
#include "r4.h"
#include "r3x3.h"

#define N 3

void r3x3_map_row (r3_t a, r3x3_t m, r3_t res)
{
  int j, k;
  for (j=0; j<N; j++)
    { double s = 0.0;
      for (k=0; k<N; k++) s += a[k] * m[k][j];
      res[j] = s;
    }
}

void r3x3_map_col (r3x3_t m, r3_t a, r3_t res)
{
  int i, k;
  for (i=0; i<N; i++)
    { double s = 0.0;
      for (k=0; k<N; k++) s += m[i][k] * a[k];
      res[i] = s;
    }
}

void r3x3_mul (r3x3_t m, r3x3_t n, r3x3_t res)
{
  int i, j, k;
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      { double s = 0.0;
        for (k=0; k<N; k++)  s += m[i][k]*n[k][j];
        res[i][j] = s;
      }
}

double r3x3_det (r3x3_t m)
{
  double d0 = m[1][1]*m[2][2] - m[1][2]*m[2][1];
  double d1 = m[1][0]*m[2][2] - m[1][2]*m[2][0];
  double d2 = m[1][0]*m[2][1] - m[1][1]*m[2][0];
  double d = d0 * m[0][0] - d1 * m[0][1] + d2 * m[0][2];
  return (d);
}

double r3x3_cof (r3x3_t m, int ix, int jx)
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
  d = m[i0][j0]*m[i1][j1] - m[i0][j1]*m[i1][j0];
  if ((ix + jx) % 2) d = -d;
  return (d);
}

void r3x3_adj (r3x3_t m, r3x3_t res)
{
  int i, j;
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      res[i][j] = r3x3_cof (m, j, i);
}

void r3x3_print (FILE *f, r3x3_t m)
{
  int i,j;
  for (i=0; i<N; i++)
    {
      for (j=0; j<N; j++) fprintf(f, "%16.8e", m[i][j]);
      fputc('\n', f);
    }
  fflush(f);
}

