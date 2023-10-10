/* See r6x6.h. */
/* Last edited on 2023-10-09 09:02:44 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <r6.h>
#include <rmxn.h>

#include <r6x6.h>

#define N 6

void r6x6_zero(r6x6_t *M)
  { int32_t i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { M->c[i][j] = 0.0; }
  }

void r6x6_ident(r6x6_t *M)
  { int32_t i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { M->c[i][j] = (i == j ? 1.0 : 0.0); }
  }

void r6x6_transp (r6x6_t *A, r6x6_t *M)
  { 
    int32_t i, j;
    for (i = 0; i < N; i++)
      { M->c[i][i] = A->c[i][i]; 
        for (j = 0; j < i; j++)
          { double a = A->c[i][j];
            double b = A->c[j][i];
            M->c[i][j] = b;
            M->c[j][i] = a;
          }
      }
  }

void r6x6_get_row(r6x6_t *A, int32_t i, r6_t *x)
  { assert((i >= 0) && (i < N));
    double *v = &(A->c[i][0]);
    int32_t j;
    for (j = 0; j < N; j++) { x->c[j] = v[j]; }
  }
  
void r6x6_set_row(r6x6_t *A, int32_t i, r6_t *x)
  { assert((i >= 0) && (i < N));
    double *v = &(A->c[i][0]);
    int32_t j;
    for (j = 0; j < N; j++) { v[j] = x->c[j]; }
  }

void r6x6_get_col(r6x6_t *A, int32_t j, r6_t *x)
  { assert((j >= 0) && (j < N));
    int32_t i;
    for (i = 0; i < N; i++) { x->c[j] = A->c[i][j]; }
  }
  
void r6x6_set_col(r6x6_t *A, int32_t j, r6_t *x)
  { assert((j >= 0) && (j < N));
    int32_t i;
    for (i = 0; i < N; i++) { A->c[i][j] = x->c[j]; }
  }

void r6x6_map_row (r6_t *x, r6x6_t *A, r6_t *r)
  { r6_t rr;
    int32_t j, k;
    for (j = 0; j < N; j++)
      { double s = 0.0;
        for (k = 0; k < N; k++) s += x->c[k] * A->c[k][j];
        rr.c[j] = s;
      }
    (*r) = rr;
  }

void r6x6_map_col (r6x6_t *A, r6_t *x, r6_t *r)
  { r6_t rr;
    int32_t i, k;
    for (i = 0; i < N; i++)
      { double s = 0.0;
        for (k = 0; k < N; k++) s += A->c[i][k] * x->c[k];
        rr.c[i] = s;
      }
    (*r) = rr;
  }

void r6x6_scale (double s, r6x6_t *A, r6x6_t *M)  
  { int32_t i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { M->c[i][j] = s * A->c[i][j]; }
  }

void r6x6_mul (r6x6_t *A, r6x6_t *B, r6x6_t *M)
  { r6x6_t RR;
    int32_t i, j, k;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { double s = 0.0;
          for (k = 0; k < N; k++)  s += A->c[i][k]*B->c[k][j];
          RR.c[i][j] = s;
        }
    (*M) = RR;
  }

void r6x6_mul_tr (r6x6_t *A, r6x6_t *B, r6x6_t *M)
  {
    r6x6_t RR;
    int32_t i, j, k;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { double s = 0.0;
          for (k = 0; k < N; k++)  s += A->c[i][k]*B->c[j][k];
          RR.c[i][j] = s;
        }
    (*M) = RR;
  }

double r6x6_det (r6x6_t *A)
  { return rmxn_det(N, &(A->c[0][0])); }

/* !!! Uncomment when we get {rmxn_adj}: !!!
void r6x6_adj (r6x6_t *A, r6x6_t *M)
  { (void)rmxn_adj(N, &(A->c[0][0]), &(M->c[0][0])); }
*/

void r6x6_inv (r6x6_t *A, r6x6_t *M)
  { (void)rmxn_inv(N, &(A->c[0][0]), &(M->c[0][0])); }

double r6x6_norm_sqr(r6x6_t* A)
  {
    double s = 0.0;
    int32_t i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { double Aij = A->c[i][j];
          s += Aij*Aij;
        }
    return s; 
  }

double r6x6_norm(r6x6_t* A)
  {
    double s = 0.0;
    int32_t i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { double Aij = A->c[i][j];
          s += Aij*Aij;
        }
    return sqrt(s); 
  }

double r6x6_normalize(r6x6_t *A)
  { 
    double w = r6x6_norm(A);
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

double r6x6_mod_norm_sqr(r6x6_t* A)
  {
    double s = 0.0;
    int32_t i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { double Aij = A->c[i][j];
          double Dij = Aij - (i == j ? 1 : 0);
          s += Dij*Dij;
        }
    return s; 
  }

bool_t r6x6_is_unif_scaling(r6x6_t *M, double s)
  {
    int32_t i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { if (M->c[i][j] != (i == j ? s : 0.0)) { return FALSE; } }
    return TRUE;
  }

void r6x6_from_rows(r6_t *a, r6_t *b, r6_t *c, r6_t *d, r6_t *e, r6_t *f, r6x6_t *M)
  { int32_t j;
    for (j = 0; j < N; j++)
      { M->c[0][j] = a->c[j];
        M->c[1][j] = b->c[j];
        M->c[2][j] = c->c[j];
        M->c[3][j] = d->c[j];
        M->c[4][j] = e->c[j];
        M->c[5][j] = f->c[j];
      }
  }

void r6x6_from_cols(r6_t *a, r6_t *b, r6_t *c, r6_t *d, r6_t *e, r6_t *f, r6x6_t *M)
  { int32_t j;
    for (j = 0; j < N; j++)
      { M->c[j][0] = a->c[j];
        M->c[j][1] = b->c[j];
        M->c[j][2] = c->c[j];
        M->c[j][3] = d->c[j];
        M->c[j][4] = e->c[j];
        M->c[j][5] = f->c[j];
      }
  }

void r6x6_print (FILE *f, r6x6_t *A)
  { r6x6_gen_print(f, A, NULL, NULL, NULL, NULL, NULL, NULL, NULL); }

void r6x6_gen_print 
  ( FILE *f, r6x6_t *A,
    char *fmt, 
    char *olp, char *osep, char *orp,
    char *ilp, char *isep, char *irp
  )
  { rmxn_gen_print(f, N, N, &(A->c[0][0]), fmt, olp, osep, orp, ilp, isep, irp); }

