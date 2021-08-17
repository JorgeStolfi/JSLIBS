/* See rf4x4.h. */
/* Last edited on 2021-08-17 08:50:48 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <r4x4.h>
#include <rmxn.h>
#include <rfmxn.h>

#include <rf4.h>
#include <rf4x4.h>

#define N 4

rf4x4_t rf4x4_zero(void)
  { rf4x4_t M;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M.c[i][j] = 0.0; }
    return M;
  }

rf4x4_t rf4x4_ident(void)
  { rf4x4_t M;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M.c[i][j] = (i == j ? 1.0 : 0.0); }
    return M;
  }

rf4x4_t rf4x4_transp (rf4x4_t *A)
  {
    rf4x4_t M;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M.c[i][j] = A->c[j][i]; }
    return M;
  }

rf4_t rf4x4_get_row(rf4x4_t *A, int32_t i)
  { rf4_t x;
    assert((i >= 0) && (i < N));
    float *v = &(A->c[i][0]);
    x.c[0] = v[0];
    x.c[1] = v[1];
    x.c[2] = v[2];
    x.c[3] = v[3];
    return x;
  }
  
void rf4x4_set_row(rf4x4_t *A, int32_t i, rf4_t *x)
  { assert((i >= 0) && (i < N));
    float *v = &(A->c[i][0]);
    v[0] = x->c[0];
    v[1] = x->c[1];
    v[2] = x->c[2];
    v[3] = x->c[3];
  }

rf4_t rf4x4_get_col(rf4x4_t *A, int32_t j)
  { rf4_t x;
    assert((j >= 0) && (j < N));
    x.c[0] = A->c[0][j];
    x.c[1] = A->c[1][j];
    x.c[2] = A->c[2][j];
    x.c[3] = A->c[3][j];
    return x;
  }
  
void rf4x4_set_col(rf4x4_t *A, int32_t j, rf4_t *x)
  { assert((j >= 0) && (j < N));
    A->c[0][j] = x->c[0];
    A->c[1][j] = x->c[1];
    A->c[2][j] = x->c[2];
    A->c[3][j] = x->c[3];
  }

rf4_t rf4x4_map_row (rf4_t *x, rf4x4_t *A)
  { rf4_t r;
    for (int32_t j = 0; j < N; j++)
      { double s = 0.0;
        for (int32_t k = 0; k < N; k++) s += (double)(x->c[k])*A->c[k][j];
        r.c[j] = (float)s;
      }
    return r;
  }

rf4_t rf4x4_map_col (rf4x4_t *A, rf4_t *x)
  { rf4_t r;
    for (int32_t i = 0; i < N; i++)
      { double s = 0.0;
        for (int32_t k = 0; k < N; k++) s += (double)(A->c[i][k])*x->c[k];
        r.c[i] = (float)s;
      }
    return r;
  }

rf4x4_t rf4x4_scale (double s, rf4x4_t *A)  
  { rf4x4_t M;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M.c[i][j] = (float)(s * A->c[i][j]); }
    return M;
  }

rf4x4_t rf4x4_mul (rf4x4_t *A, rf4x4_t *B)
  { rf4x4_t M;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { double s = 0.0;
          for (int32_t k = 0; k < N; k++)  s += (double)(A->c[i][k])*B->c[k][j];
          M.c[i][j] = (float)s;
        }
    return M;
  }

rf4x4_t rf4x4_mul_tr (rf4x4_t *A, rf4x4_t *B)
  {
    rf4x4_t M;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { double s = 0.0;
          for (int32_t k = 0; k < N; k++)  s += (double)(A->c[i][k])*B->c[j][k];
          M.c[i][j] = (float)s;
        }
    return M;
  }

double rf4x4_det (rf4x4_t *A)
  { double A00 = A->c[0][0];
    double A01 = A->c[0][1];
    double A02 = A->c[0][2];
    double A03 = A->c[0][3];

    double A10 = A->c[1][0];
    double A11 = A->c[1][1];
    double A12 = A->c[1][2];
    double A13 = A->c[1][3];

    double A01_01 = A00 * A11 - A01 * A10;
    double A01_02 = A00 * A12 - A02 * A10;
    double A01_12 = A01 * A12 - A02 * A11;
    double A01_03 = A00 * A13 - A03 * A10;
    double A01_13 = A01 * A13 - A03 * A11;
    double A01_23 = A02 * A13 - A03 * A12;

    double A20 = A->c[2][0];
    double A21 = A->c[2][1];
    double A22 = A->c[2][2];
    double A23 = A->c[2][3];

    double A30 = A->c[3][0];
    double A31 = A->c[3][1];
    double A32 = A->c[3][2];
    double A33 = A->c[3][3];

    double A23_01 = A20 * A31 - A21 * A30;
    double A23_02 = A20 * A32 - A22 * A30;
    double A23_12 = A21 * A32 - A22 * A31;
    double A23_03 = A20 * A33 - A23 * A30;
    double A23_13 = A21 * A33 - A23 * A31;
    double A23_23 = A22 * A33 - A23 * A32;

    double d = 
        A01_01 * A23_23 - A01_02 * A23_13 + A01_12 * A23_03
      + A01_03 * A23_12 - A01_13 * A23_02 + A01_23 * A23_01;

    return d;
  }

rf4x4_t rf4x4_adj (rf4x4_t *A)
  { r4x4_t C;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { C.c[i][j] = (double)A->c[i][j]; }
    r4x4_adj(&C, &C);
    rf4x4_t M;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M.c[i][j] = (float)C.c[i][j]; }
    return M;
  }

rf4x4_t rf4x4_inv (rf4x4_t *A)
  { r4x4_t C;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { C.c[i][j] = (double)A->c[i][j]; }
    r4x4_inv(&C, &C);
    rf4x4_t M;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M.c[i][j] = (float)C.c[i][j]; }
    return M;
  }

double rf4x4_norm_sqr(rf4x4_t* A)
  {
    double A00 = A->c[0][0];
    double A01 = A->c[0][1];
    double A02 = A->c[0][2];
    double A03 = A->c[0][3];

    double A10 = A->c[1][0];
    double A11 = A->c[1][1];
    double A12 = A->c[1][2];
    double A13 = A->c[1][3];

    double A20 = A->c[2][0];
    double A21 = A->c[2][1];
    double A22 = A->c[2][2];
    double A23 = A->c[2][3];

    double A30 = A->c[3][0];
    double A31 = A->c[3][1];
    double A32 = A->c[3][2];
    double A33 = A->c[3][3];

    return 
      A00*A00 + A01*A01 + A02*A02 + A03*A03 + 
      A10*A10 + A11*A11 + A12*A12 + A13*A13 + 
      A20*A20 + A21*A21 + A22*A22 + A23*A23 +
      A30*A30 + A31*A31 + A32*A32 + A33*A33; 
  }

double rf4x4_norm(rf4x4_t* A)
  {
    double A00 = A->c[0][0];
    double A01 = A->c[0][1];
    double A02 = A->c[0][2];
    double A03 = A->c[0][3];

    double A10 = A->c[1][0];
    double A11 = A->c[1][1];
    double A12 = A->c[1][2];
    double A13 = A->c[1][3];

    double A20 = A->c[2][0];
    double A21 = A->c[2][1];
    double A22 = A->c[2][2];
    double A23 = A->c[2][3];

    double A30 = A->c[3][0];
    double A31 = A->c[3][1];
    double A32 = A->c[3][2];
    double A33 = A->c[3][3];

    return sqrt
      ( A00*A00 + A01*A01 + A02*A02 + A03*A03 + 
        A10*A10 + A11*A11 + A12*A12 + A13*A13 + 
        A20*A20 + A21*A21 + A22*A22 + A23*A23 +
        A30*A30 + A31*A31 + A32*A32 + A33*A33
      ); 
  }

double rf4x4_mod_norm_sqr(rf4x4_t* A)
  {
    double D00 = A->c[0][0] - 1;
    double D01 = A->c[0][1];
    double D02 = A->c[0][2];
    double D03 = A->c[0][3];

    double D10 = A->c[1][0];
    double D11 = A->c[1][1] - 1;
    double D12 = A->c[1][2];
    double D13 = A->c[1][3];

    double D20 = A->c[2][0];
    double D21 = A->c[2][1];
    double D22 = A->c[2][2] - 1;
    double D23 = A->c[2][3];

    double D30 = A->c[3][0];
    double D31 = A->c[3][1];
    double D32 = A->c[3][2];
    double D33 = A->c[3][3] - 1;

    return 
      D00*D00 + D01*D01 + D02*D02 + D03*D03 + 
      D10*D10 + D11*D11 + D12*D12 + D13*D13 + 
      D20*D20 + D21*D21 + D22*D22 + D23*D23 +
      D30*D30 + D31*D31 + D32*D32 + D33*D33; 
  }

bool_t rf4x4_is_unif_scaling(rf4x4_t *M, double s)
  {
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { if (M->c[i][j] != (i == j ? s : 0.0)) { return FALSE; } }
    return TRUE;
  }

rf4x4_t rf4x4_from_rows(rf4_t *a, rf4_t *b, rf4_t *c, rf4_t *d)
  { rf4x4_t M;
    for (int32_t j = 0; j < N; j++)
      { M.c[0][j] = a->c[j];
        M.c[1][j] = b->c[j];
        M.c[2][j] = c->c[j];
        M.c[3][j] = d->c[j];
      }
    return M;
  }

rf4x4_t rf4x4_from_cols(rf4_t *a, rf4_t *b, rf4_t *c, rf4_t *d)
  { rf4x4_t M;
    for (int32_t j = 0; j < N; j++)
      { M.c[j][0] = a->c[j];
        M.c[j][1] = b->c[j];
        M.c[j][2] = c->c[j];
        M.c[j][3] = d->c[j];
      }
    return M;
  }

void rf4x4_print (FILE *f, rf4x4_t *A)
  { rf4x4_gen_print(f, A, NULL, NULL, NULL, NULL, NULL, NULL, NULL); }

void rf4x4_gen_print
  ( FILE *f, rf4x4_t *A,
    char *fmt, 
    char *olp, char *osep, char *orp,
    char *ilp, char *isep, char *irp
  )
  { rfmxn_gen_print(f, N, N, &(A->c[0][0]), fmt, olp, osep, orp, ilp, isep, irp); }

