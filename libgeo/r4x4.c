/* See r4x4.h. */
/* Last edited on 2024-11-20 16:21:11 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <r4.h>
#include <rmxn.h>
#include <jsrandom.h>

#include <r4x4.h>

#define N 4

void r4x4_zero(r4x4_t *M)
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = 0.0; }
  }

void r4x4_ident(r4x4_t *M)
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = (i == j ? 1.0 : 0.0); }
  }

void r4x4_throw(r4x4_t *M, sign_t sgn)
  { assert((sgn >= -1) && (sgn <= +1));
    while (TRUE)
      { for (int32_t i = 0; i < N; i++)
          { for (int32_t j = 0; j < N; j++) 
              { M->c[i][j] = 2*drandom() - 1; }
          }
        if (sgn == 0) { break; }
        double det = r4x4_det(M);
        if (det == 0) { continue; }
        if (det*sgn < 0) 
          { /* Negate first row: */
            for (int32_t j = 0; j < N; j++) { M->c[0][j] = - M->c[0][j]; }
          }
        /* At this point, {sgn*det} must be positive: */
        break;
      }
  }

void r4x4_transp(r4x4_t *A, r4x4_t *M)
  {
    double a, b;
    /* Copy diagonal elements: */
    M->c[0][0] = A->c[0][0];
    M->c[1][1] = A->c[1][1];
    M->c[2][2] = A->c[2][2];
    M->c[3][3] = A->c[3][3];
    /* Copy remaining elements, transposed. */
    a = A->c[0][1]; b = A->c[1][0]; M->c[0][1] = b; M->c[1][0] = a;
    a = A->c[0][2]; b = A->c[2][0]; M->c[0][2] = b; M->c[2][0] = a;
    a = A->c[0][3]; b = A->c[3][0]; M->c[0][3] = b; M->c[3][0] = a;
    a = A->c[1][2]; b = A->c[2][1]; M->c[1][2] = b; M->c[2][1] = a;
    a = A->c[1][3]; b = A->c[3][1]; M->c[1][3] = b; M->c[3][1] = a;
    a = A->c[2][3]; b = A->c[3][2]; M->c[2][3] = b; M->c[3][2] = a;
  }

void r4x4_get_row(r4x4_t *A, uint32_t i, r4_t *x)
  { assert(i < N);
    double *v = &(A->c[i][0]);
    x->c[0] = v[0];
    x->c[1] = v[1];
    x->c[2] = v[2];
    x->c[3] = v[3];
  }
  
void r4x4_set_row(r4x4_t *A, uint32_t i, r4_t *x)
  { assert(i < N);
    double *v = &(A->c[i][0]);
    v[0] = x->c[0];
    v[1] = x->c[1];
    v[2] = x->c[2];
    v[3] = x->c[3];
  }

void r4x4_get_col(r4x4_t *A, uint32_t j, r4_t *x)
  { assert(j < N);
    x->c[0] = A->c[0][j];
    x->c[1] = A->c[1][j];
    x->c[2] = A->c[2][j];
    x->c[3] = A->c[3][j];
  }
  
void r4x4_set_col(r4x4_t *A, uint32_t j, r4_t *x)
  { assert(j < N);
    A->c[0][j] = x->c[0];
    A->c[1][j] = x->c[1];
    A->c[2][j] = x->c[2];
    A->c[3][j] = x->c[3];
  }

void r4x4_map_row(r4_t *x, r4x4_t *A, r4_t *r)
  { r4_t rr;
    for (int32_t j = 0; j < N; j++)
      { double s = 0.0;
        for (int32_t k = 0; k < N; k++) s += x->c[k] * A->c[k][j];
        rr.c[j] = s;
      }
    (*r) = rr;
  }

void r4x4_map_col(r4x4_t *A, r4_t *x, r4_t *r)
  { r4_t rr;
    for (int32_t i = 0; i < N; i++)
      { double s = 0.0;
        for (int32_t k = 0; k < N; k++) s += A->c[i][k] * x->c[k];
        rr.c[i] = s;
      }
    (*r) = rr;
  }

void r4x4_add(r4x4_t *A, r4x4_t *B, r4x4_t *M)
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = A->c[i][j] + B->c[i][j]; }
  }

void r4x4_sub(r4x4_t *A, r4x4_t *B, r4x4_t *M)
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = A->c[i][j] - B->c[i][j]; }
  }

void r4x4_neg(r4x4_t *A, r4x4_t *M)
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = - A->c[i][j]; }
  }

void r4x4_scale(double s, r4x4_t *A, r4x4_t *M)  
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = s * A->c[i][j]; }
  }

void r4x4_mul(r4x4_t *A, r4x4_t *B, r4x4_t *M)
  { r4x4_t RR;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { double s = 0.0;
          for (int32_t k = 0; k < N; k++)  s += A->c[i][k]*B->c[k][j];
          RR.c[i][j] = s;
        }
    (*M) = RR;
  }

void r4x4_mul_tr(r4x4_t *A, r4x4_t *B, r4x4_t *M)
  {
    r4x4_t RR;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { double s = 0.0;
          for (int32_t k = 0; k < N; k++)  s += A->c[i][k]*B->c[j][k];
          RR.c[i][j] = s;
        }
    (*M) = RR;
  }

double r4x4_det(r4x4_t *A)
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

void r4x4_adj(r4x4_t *A, r4x4_t *M)
  { double A00 = A->c[0][0];
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

    double A01_01 = A00 * A11 - A01 * A10;
    double A01_02 = A00 * A12 - A02 * A10;
    double A01_12 = A01 * A12 - A02 * A11;
    double A01_03 = A00 * A13 - A03 * A10;
    double A01_13 = A01 * A13 - A03 * A11;
    double A01_23 = A02 * A13 - A03 * A12;

    double A23_01 = A20 * A31 - A21 * A30;
    double A23_02 = A20 * A32 - A22 * A30;
    double A23_12 = A21 * A32 - A22 * A31;
    double A23_03 = A20 * A33 - A23 * A30;
    double A23_13 = A21 * A33 - A23 * A31;
    double A23_23 = A22 * A33 - A23 * A32;

    double A012_012 = A01_01 * A22 - A01_02 * A21 + A01_12 * A20;
    double A012_013 = A01_01 * A23 - A01_03 * A21 + A01_13 * A20;
    double A012_023 = A01_02 * A23 - A01_03 * A22 + A01_23 * A20;
    double A012_123 = A01_12 * A23 - A01_13 * A22 + A01_23 * A21;

    double A013_012 = A01_01 * A32 - A01_02 * A31 + A01_12 * A30;
    double A013_013 = A01_01 * A33 - A01_03 * A31 + A01_13 * A30;
    double A013_023 = A01_02 * A33 - A01_03 * A32 + A01_23 * A30;
    double A013_123 = A01_12 * A33 - A01_13 * A32 + A01_23 * A31;

    double A023_012 = A00 * A23_12 - A01 * A23_02 + A02 * A23_01;
    double A023_013 = A00 * A23_13 - A01 * A23_03 + A03 * A23_01;
    double A023_023 = A00 * A23_23 - A02 * A23_03 + A03 * A23_02;
    double A023_123 = A01 * A23_23 - A02 * A23_13 + A03 * A23_12;

    double A123_012 = A10 * A23_12 - A11 * A23_02 + A12 * A23_01;
    double A123_013 = A10 * A23_13 - A11 * A23_03 + A13 * A23_01;
    double A123_023 = A10 * A23_23 - A12 * A23_03 + A13 * A23_02;
    double A123_123 = A11 * A23_23 - A12 * A23_13 + A13 * A23_12;

    M->c[0][0] =  A123_123; 
    M->c[0][1] = -A023_123; 
    M->c[0][2] =  A013_123; 
    M->c[0][3] = -A012_123;

    M->c[1][0] = -A123_023; 
    M->c[1][1] =  A023_023; 
    M->c[1][2] = -A013_023; 
    M->c[1][3] =  A012_023;

    M->c[2][0] =  A123_013; 
    M->c[2][1] = -A023_013; 
    M->c[2][2] =  A013_013; 
    M->c[2][3] = -A012_013;

    M->c[3][0] = -A123_012; 
    M->c[3][1] =  A023_012; 
    M->c[3][2] = -A013_012; 
    M->c[3][3] =  A012_012;
  }

void r4x4_inv(r4x4_t *A, r4x4_t *M)
  { r4x4_t RR;
    r4x4_adj(A, &RR);
    double d = 
        A->c[0][0]*RR.c[0][0]
      + A->c[0][1]*RR.c[1][0] 
      + A->c[0][2]*RR.c[2][0] 
      + A->c[0][3]*RR.c[3][0];
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = RR.c[i][j]/d; }
  }

double r4x4_norm_sqr(r4x4_t *A)
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

double r4x4_norm(r4x4_t *A)
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

double r4x4_normalize(r4x4_t *A)
  { 
    double w = r4x4_norm(A);
    if (w != 0)
      { for (int32_t i = 0; i < N; i++)
          { for (int32_t j = 0; j < N; j++)
             { A->c[i][j] /= w; }
          }
      }
    else
      { for (int32_t i = 0; i < N; i++)
          { for (int32_t j = 0; j < N; j++)
             { A->c[i][j] = NAN; }
          }
      }
    return w;
  }

double r4x4_mod_norm_sqr(r4x4_t *A)
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

bool_t r4x4_is_unif_scaling(r4x4_t *M, double s)
  {
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { if (M->c[i][j] != (i == j ? s : 0.0)) { return FALSE; } }
    return TRUE;
  }

void r4x4_from_rows(r4_t *a, r4_t *b, r4_t *c, r4_t *d, r4x4_t *M)
  { for (int32_t j = 0; j < N; j++)
      { M->c[0][j] = a->c[j];
        M->c[1][j] = b->c[j];
        M->c[2][j] = c->c[j];
        M->c[3][j] = d->c[j];
      }
  }

void r4x4_from_cols(r4_t *a, r4_t *b, r4_t *c, r4_t *d, r4x4_t *M)
  { for (int32_t j = 0; j < N; j++)
      { M->c[j][0] = a->c[j];
        M->c[j][1] = b->c[j];
        M->c[j][2] = c->c[j];
        M->c[j][3] = d->c[j];
      }
  }

void r4x4_print(FILE *f, r4x4_t *A)
  { r4x4_gen_print(f, A, NULL, NULL, NULL, NULL, NULL, NULL, NULL); }

void r4x4_gen_print
  ( FILE *f, r4x4_t *A,
    char *fmt, 
    char *olp, char *osep, char *orp,
    char *ilp, char *isep, char *irp
  )
  { 
    rmxn_gen_print(f, N, N, &(A->c[0][0]), fmt, olp, osep, orp, ilp, isep, irp);
  }

