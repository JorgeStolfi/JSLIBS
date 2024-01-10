/* See r2x2.h. */
/* Last edited on 2001-10-21 21:15:36 by stolfi */

#include <math.h>
#include <stdio.h>
#include <ioprotos.h>
#include "r2.h"
#include "r2x2.h"

#define N 2

void r2x2_map_row (r2_t a, r2x2_t m, r2_t res)
{
  res[0] = a[0] * m[0][0] + a[1] * m[1][0];
  res[1] = a[0] * m[0][1] + a[1] * m[1][1];
}

void r2x2_map_col (r2x2_t m, r2_t a, r2_t res)
{
  res[0] = m[0][0] * a[0] + m[0][1] * a[1];
  res[1] = m[1][0] * a[0] + m[1][1] * a[1];
}

void r2x2_mul (r2x2_t m, r2x2_t n, r2x2_t res)
{
  int i, j, k;
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      { double s = 0.0;
        for (k=0; k<N; k++)  s += m[i][k]*n[k][j];
        res[i][j] = s;
      }
}

double r2x2_det (r2x2_t m)
{
  return m[0][0] * m[1][1] - m[1][0] * m[0][1];
}

double r2x2_cof (r2x2_t m, int ix, int jx)
{
  int ic = 1-ix, jc = 1-jx;
  double d = m[ic][jc];
  if (ix != jx) { d = -d; }
  return d;
}

void r2x2_adj (r2x2_t m, r2x2_t res)
{
  res[0][0] =   m[1][1];
  res[0][1] = - m[0][1];
  res[1][0] = - m[1][0];
  res[1][1] =   m[0][0];
}

void r2x2_print (FILE *f, r2x2_t m)
{
  int i,j;
  for (i=0; i<N; i++)
    {
      for (j=0; j<N; j++) fprintf(f, "%16.8e", m[i][j]);
      fputc('\n', f);
    }
  fflush(f);
}

