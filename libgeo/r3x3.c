/* See r3x3.h. */
/* Last edited on 2024-11-07 23:50:05 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <r3.h>
#include <affirm.h>
#include <rmxn.h>
#include <sign.h>

#include <r3x3.h>

#define N 3

void r3x3_zero(r3x3_t *M)
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = 0.0; }
  }

void r3x3_ident(r3x3_t *M)
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = (i == j ? 1.0 : 0.0); }
  }

void r3x3_throw(r3x3_t *A, sign_t sgn)
  { demand((sgn >= -1) && (sgn <= +1), "invalid {sgn}");
    while (TRUE)
      { for (int32_t i = 0; i < N; i++)
          { for (int32_t j = 0; j < N; j++) 
              { A->c[i][j] += 2*drand48() - 1; }
          }
        if (sgn == 0) { break; }
        double det = r3x3_det(A);
        if (det == 0) { continue; }
        if (det*sgn < 0) 
          { /* Negate first row: */
            for (int32_t j = 0; j < N; j++) { A->c[0][j] = - A->c[0][j]; }
          }
        /* At this point, {sgn*det} must be positive: */
        break;
      }
  }

void r3x3_transp(r3x3_t *A, r3x3_t *M)
  {
    double a, b;
    /* Copy diagonal elements: */
    M->c[0][0] = A->c[0][0];
    M->c[1][1] = A->c[1][1];
    M->c[2][2] = A->c[2][2];
    /* Copy remaining elements, transposed. */
    a = A->c[0][1]; b = A->c[1][0]; M->c[0][1] = b; M->c[1][0] = a;
    a = A->c[0][2]; b = A->c[2][0]; M->c[0][2] = b; M->c[2][0] = a;
    a = A->c[1][2]; b = A->c[2][1]; M->c[1][2] = b; M->c[2][1] = a;
  }

void r3x3_get_row(r3x3_t *A, int32_t i, r3_t *x)
  { assert((i >= 0) && (i < N));
    double *v = &(A->c[i][0]);
    x->c[0] = v[0];
    x->c[1] = v[1];
    x->c[2] = v[2];
  }
  
void r3x3_set_row(r3x3_t *A, int32_t i, r3_t *x)
  { assert((i >= 0) && (i < N));
    double *v = &(A->c[i][0]);
    v[0] = x->c[0];
    v[1] = x->c[1];
    v[2] = x->c[2];
  }

void r3x3_get_col(r3x3_t *A, int32_t j, r3_t *x)
  { assert((j >= 0) && (j < N));
    x->c[0] = A->c[0][j];
    x->c[1] = A->c[1][j];
    x->c[2] = A->c[2][j];
  }
  
void r3x3_set_col(r3x3_t *A, int32_t j, r3_t *x)
  { assert((j >= 0) && (j < N));
    A->c[0][j] = x->c[0];
    A->c[1][j] = x->c[1];
    A->c[2][j] = x->c[2];
  }

void r3x3_map_row(r3_t *x, r3x3_t *A, r3_t *r)
  { r3_t rr;
    for (int32_t j = 0; j < N; j++)
      { double s = 0.0;
        for (int32_t k = 0; k < N; k++) s += x->c[k] * A->c[k][j];
        rr.c[j] = s;
      }
    for (int32_t j = 0; j < N; j++) { r->c[j] = rr.c[j]; }
  }

void r3x3_map_col(r3x3_t *A, r3_t *x, r3_t *r)
  { r3_t rr;
    for (int32_t i = 0; i < N; i++)
      { double s = 0.0;
        for (int32_t k = 0; k < N; k++) s += A->c[i][k] * x->c[k];
        rr.c[i] = s;
      }
    for (int32_t i = 0; i < N; i++) { r->c[i] = rr.c[i]; }
  }

void r3x3_add(r3x3_t *A, r3x3_t *B, r3x3_t *M)
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = A->c[i][j] + B->c[i][j]; }
  }

void r3x3_sub(r3x3_t *A, r3x3_t *B, r3x3_t *M)
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = A->c[i][j] - B->c[i][j]; }
  }

void r3x3_neg(r3x3_t *A, r3x3_t *M)
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = - A->c[i][j]; }
  }

void r3x3_scale(double s, r3x3_t *A, r3x3_t *M)  
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = s * A->c[i][j]; }
  }

void r3x3_mix(double s, r3x3_t *A, double t, r3x3_t *B, r3x3_t *M)
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { M->c[i][j] = s * A->c[i][j] + t * B->c[i][j]; }
  }

void r3x3_mul(r3x3_t *A, r3x3_t *B, r3x3_t *M)
  { r3x3_t RR;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { double s = 0.0;
          for (int32_t k = 0; k < N; k++)  s += A->c[i][k]*B->c[k][j];
          RR.c[i][j] = s;
        }
    (*M) = RR;
  }

void r3x3_mul_tr(r3x3_t *A, r3x3_t *B, r3x3_t *M)
  {
    r3x3_t RR;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { double s = 0.0;
          for (int32_t k = 0; k < N; k++)  s += A->c[i][k]*B->c[j][k];
          RR.c[i][j] = s;
        }
    (*M) = RR;
  }
  
void r3x3_tr_mul(r3x3_t *A, r3x3_t *B, r3x3_t *M)
  {
    r3x3_t RR;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { double s = 0.0;
          for (int32_t k = 0; k < N; k++)  s += A->c[k][i]*B->c[k][j];
          RR.c[i][j] = s;
        }
    (*M) = RR;
  }

double r3x3_det(r3x3_t *A)
  { double d0 = A->c[1][1]*A->c[2][2] - A->c[1][2]*A->c[2][1];
    double d1 = A->c[1][0]*A->c[2][2] - A->c[1][2]*A->c[2][0];
    double d2 = A->c[1][0]*A->c[2][1] - A->c[1][1]*A->c[2][0];
    double det = d0 * A->c[0][0] - d1 * A->c[0][1] + d2 * A->c[0][2];
    return det;
  }

void r3x3_adj(r3x3_t *A, r3x3_t *M)
  { double A00 = A->c[0][0];
    double A01 = A->c[0][1];
    double A02 = A->c[0][2];

    double A10 = A->c[1][0];
    double A11 = A->c[1][1];
    double A12 = A->c[1][2];

    double A20 = A->c[2][0];
    double A21 = A->c[2][1];
    double A22 = A->c[2][2];

    double A01_01 = A00 * A11 - A01 * A10;
    double A01_02 = A00 * A12 - A02 * A10;
    double A01_12 = A01 * A12 - A02 * A11;

    double A02_01 = A00 * A21 - A01 * A20;
    double A02_02 = A00 * A22 - A02 * A20;
    double A02_12 = A01 * A22 - A02 * A21;

    double A12_01 = A10 * A21 - A11 * A20;
    double A12_02 = A10 * A22 - A12 * A20;
    double A12_12 = A11 * A22 - A12 * A21;

    M->c[0][0] =  A12_12; 
    M->c[0][1] = -A02_12; 
    M->c[0][2] =  A01_12;

    M->c[1][0] = -A12_02; 
    M->c[1][1] =  A02_02; 
    M->c[1][2] = -A01_02;

    M->c[2][0] =  A12_01; 
    M->c[2][1] = -A02_01; 
    M->c[2][2] =  A01_01;
  }

void r3x3_inv(r3x3_t *A, r3x3_t *M)
  { r3x3_t RR;
    r3x3_adj(A, &RR);
    double d = 
        A->c[0][0]*RR.c[0][0]
      + A->c[0][1]*RR.c[1][0]
      + A->c[0][2]*RR.c[2][0];
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        M->c[i][j] = RR.c[i][j]/d;
  }

double r3x3_norm_sqr(r3x3_t *A)
  {
    double A00 = A->c[0][0];
    double A01 = A->c[0][1];
    double A02 = A->c[0][2];

    double A10 = A->c[1][0];
    double A11 = A->c[1][1];
    double A12 = A->c[1][2];

    double A20 = A->c[2][0];
    double A21 = A->c[2][1];
    double A22 = A->c[2][2];

    return A00*A00 + A01*A01 + A02*A02 + A10*A10 + A11*A11 + A12*A12 + A20*A20 + A21*A21 + A22*A22; 
  }

double r3x3_norm(r3x3_t *A)
  {
    double A00 = A->c[0][0];
    double A01 = A->c[0][1];
    double A02 = A->c[0][2];

    double A10 = A->c[1][0];
    double A11 = A->c[1][1];
    double A12 = A->c[1][2];

    double A20 = A->c[2][0];
    double A21 = A->c[2][1];
    double A22 = A->c[2][2];

    return sqrt( A00*A00 + A01*A01 + A02*A02 + A10*A10 + A11*A11 + A12*A12 + A20*A20 + A21*A21 + A22*A22 ); 
  }

double r3x3_normalize(r3x3_t *A)
  { 
    double w = r3x3_norm(A);
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

double r3x3_mod_norm_sqr(r3x3_t *A)
  {
    double D00 = A->c[0][0] - 1;
    double D01 = A->c[0][1];
    double D02 = A->c[0][2];

    double D10 = A->c[1][0];
    double D11 = A->c[1][1] - 1;
    double D12 = A->c[1][2];

    double D20 = A->c[2][0];
    double D21 = A->c[2][1];
    double D22 = A->c[2][2] - 1;

    return 
      D00*D00 + D01*D01 + D02*D02 + 
      D10*D10 + D11*D11 + D12*D12 + 
      D20*D20 + D21*D21 + D22*D22; 
  }

double r3x3_L_inf_norm(r3x3_t *A)
  {
    double mv = 0.0;
    mv = fmax(mv, fabs(A->c[0][0]));
    mv = fmax(mv, fabs(A->c[0][1]));
    mv = fmax(mv, fabs(A->c[0][2]));
    
    mv = fmax(mv, fabs(A->c[1][0]));
    mv = fmax(mv, fabs(A->c[1][1]));
    mv = fmax(mv, fabs(A->c[1][2]));
    
    mv = fmax(mv, fabs(A->c[2][0]));
    mv = fmax(mv, fabs(A->c[2][1]));
    mv = fmax(mv, fabs(A->c[2][2]));
    
    return mv;
  }

double r3x3_L_inf_normalize(r3x3_t *A)
  { 
    double mv = r3x3_L_inf_norm(A);
    A->c[0][0] /= mv;
    A->c[0][1] /= mv;
    A->c[0][2] /= mv;
              
    A->c[1][0] /= mv;
    A->c[1][1] /= mv;
    A->c[1][2] /= mv;
              
    A->c[2][0] /= mv;
    A->c[2][1] /= mv;
    A->c[2][2] /= mv;
    
    return mv;
  }

void r3x3_diff_sqr(r3x3_t *A, r3x3_t *B, r3x3_t *R, double *dabs2P, double *drel2P)
  {
    double dabs2 = 0.0;
    double drel2 = 0.0;
    for (int32_t i = 0; i < N; i++) 
      { for (int32_t j = 0; j < N; j++) 
          { double *Rp = &(R->c[i][j]);
            if ((*Rp) != 0.0)
              { double Re = (*Rp);
                double *Ap = &(A->c[i][j]); double Ae = (*Ap);
                double *Bp = &(B->c[i][j]); double Be = (*Bp);
                double d = Ae - Be;
                dabs2 += d*d;
                drel2 += (d/Re)*(d/Re);
              }
          }
      }
    if (dabs2P != NULL) { (*dabs2P) = dabs2; }
    if (drel2P != NULL) { (*drel2P) = drel2; }
  }

bool_t r3x3_is_unif_scaling(r3x3_t *M, double s)
  { for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { if (M->c[i][j] != (i == j ? s : 0.0)) { return FALSE; } }
    return TRUE;
  }

void r3x3_from_rows(r3_t *a, r3_t *b, r3_t *c, r3x3_t *M)
  { for (int32_t j = 0; j < N; j++)
      { M->c[0][j] = a->c[j];
        M->c[1][j] = b->c[j];
        M->c[2][j] = c->c[j];
      }
  }

void r3x3_from_cols(r3_t *a, r3_t *b, r3_t *c, r3x3_t *M)
  { for (int32_t j = 0; j < N; j++)
      { M->c[j][0] = a->c[j];
        M->c[j][1] = b->c[j];
        M->c[j][2] = c->c[j];
      }
  }

void r3x3_u_v_rotation(r3_t *u, r3_t *v, r3x3_t *M)
  { r3_t e;
    double cost, sint;
    cost = r3_dot(u, v);
    r3_cross(u,v, &e);
    sint = r3_norm(&e);
    r3_dir(&e, &e);
    /* We should have {cost*cost + sint*sint == 1.0} or nearly so. */
    if ((fabs(cost + 1.0) < 1.0e-14) || (fabs(sint) < 1.0e-14))
      { /* Vectors are practically collinear but opposite. */
        /* Pick any {e} orthogonal to {u}. */
        { unsigned k = 0; /* will be index of smallest {u} coordinate. */
          if (fabs(u->c[1]) < fabs(u->c[k])) { k = 1; }
          if (fabs(u->c[2]) < fabs(u->c[k])) { k = 2; }
          e = (r3_t){{0.0, 0.0, 0.0}};
          e.c[k] = 1.0;;
        }
        r3_t ue; r3_cross(u,&e,&ue); r3_dir(&ue, &e);
        cost = -1; sint = 0.0;
      }
    double x = e.c[0], y = e.c[1], z = e.c[2];
    double xx = x*x, yy = y*y, zz = z*z;
    double xy = x*y, yz = y*z, zx = z*x;
    
    M->c[0][0] = xx + cost*(yy + zz);
    M->c[0][1] = (1.0 - cost)*xy + sint*z;
    M->c[0][2] = (1.0 - cost)*zx - sint*y;
    M->c[1][0] = (1.0 - cost)*xy - sint*z;
    M->c[1][1] = yy + cost*(xx + zz);
    M->c[1][2] = (1.0 - cost)*yz + sint*x;
    M->c[2][0] = (1.0 - cost)*zx + sint*y;
    M->c[2][1] = (1.0 - cost)*yz - sint*x;
    M->c[2][2] = zz + cost*(xx + yy);
  }

void r3x3_print(FILE *f, r3x3_t *A)
  { r3x3_gen_print(f, A, NULL, NULL, NULL, NULL, NULL, NULL, NULL); }

void r3x3_gen_print 
  ( FILE *f, r3x3_t *A,
    char *fmt, 
    char *olp, char *osep, char *orp,
    char *ilp, char *isep, char *irp
  )
  { rmxn_gen_print(f, N, N, &(A->c[0][0]), fmt, olp, osep, orp, ilp, isep, irp); }

