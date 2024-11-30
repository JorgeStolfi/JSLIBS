/* See r2x2.h. */
/* Last edited on 2024-11-20 12:52:03 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <r2.h>
#include <affirm.h>
#include <rmxn.h>

#include <r2x2.h>

#define N 2

void r2x2_zero(r2x2_t *M)
  { M->c[0][0] = 0.0;
    M->c[0][1] = 0.0;
    M->c[1][0] = 0.0;
    M->c[1][1] = 0.0;
  }

void r2x2_ident(r2x2_t *M)
  { M->c[0][0] = 1.0;
    M->c[0][1] = 0.0;
    M->c[1][0] = 0.0;
    M->c[1][1] = 1.0;
  }
    
void r2x2_transp(r2x2_t *A, r2x2_t *M)
  {
    double a, b;
    /* Copy diagonal elements: */
    M->c[0][0] = A->c[0][0];
    M->c[1][1] = A->c[1][1];
    /* Copy remaining elements, transposed. */
    a = A->c[0][1]; b = A->c[1][0]; M->c[0][1] = b; M->c[1][0] = a;
  }

void r2x2_get_row(r2x2_t *A, uint32_t i, r2_t *x)
  { assert((i >= 0) && (i < N));
    double *v = &(A->c[i][0]);
    x->c[0] = v[0];
    x->c[1] = v[1];
  }
  
void r2x2_set_row(r2x2_t *A, uint32_t i, r2_t *x)
  { assert((i >= 0) && (i < N));
    double *v = &(A->c[i][0]);
    v[0] = x->c[0];
    v[1] = x->c[1];
  }

void r2x2_get_col(r2x2_t *A, uint32_t j, r2_t *x)
  { assert((j >= 0) && (j < N));
    x->c[0] = A->c[0][j];
    x->c[1] = A->c[1][j];
  }
  
void r2x2_set_col(r2x2_t *A, uint32_t j, r2_t *x)
  { assert((j >= 0) && (j < N));
    A->c[0][j] = x->c[0];
    A->c[1][j] = x->c[1];
  }

void r2x2_map_row(r2_t *x, r2x2_t *A, r2_t *r)
  { r2_t rr;
    rr.c[0] = x->c[0] * A->c[0][0] + x->c[1] * A->c[1][0];
    rr.c[1] = x->c[0] * A->c[0][1] + x->c[1] * A->c[1][1];
    
    r->c[0] = rr.c[0];
    r->c[1] = rr.c[1];
  }

void r2x2_map_col(r2x2_t *A, r2_t *x, r2_t *r)
  { r2_t rr;
    rr.c[0] = A->c[0][0] * x->c[0] + A->c[0][1] * x->c[1];
    rr.c[1] = A->c[1][0] * x->c[0] + A->c[1][1] * x->c[1];
    
    r->c[0] = rr.c[0];
    r->c[1] = rr.c[1];
  }
  
void r2x2_add(r2x2_t *A, r2x2_t *B, r2x2_t *M)
  { M->c[0][0] = A->c[0][0] + B->c[0][0];
    M->c[0][1] = A->c[0][1] + B->c[0][1];
    M->c[1][0] = A->c[1][0] + B->c[1][0];
    M->c[1][1] = A->c[1][1] + B->c[1][1];
  }

void r2x2_sub(r2x2_t *A, r2x2_t *B, r2x2_t *M)
  { M->c[0][0] = A->c[0][0] - B->c[0][0];
    M->c[0][1] = A->c[0][1] - B->c[0][1];
    M->c[1][0] = A->c[1][0] - B->c[1][0];
    M->c[1][1] = A->c[1][1] - B->c[1][1];
  }

void r2x2_neg(r2x2_t *A, r2x2_t *M)
  { M->c[0][0] = - A->c[0][0];
    M->c[0][1] = - A->c[0][1];
    M->c[1][0] = - A->c[1][0];
    M->c[1][1] = - A->c[1][1];
  }

void r2x2_scale(double s, r2x2_t *A, r2x2_t *M)  
  { M->c[0][0] = s * A->c[0][0];
    M->c[0][1] = s * A->c[0][1];
    M->c[1][0] = s * A->c[1][0];
    M->c[1][1] = s * A->c[1][1];
  }

void r2x2_mix(double s, r2x2_t *A, double t, r2x2_t *B, r2x2_t *M)
  { M->c[0][0] = s * A->c[0][0] + t * B->c[0][0];
    M->c[0][1] = s * A->c[0][1] + t * B->c[0][1];
    M->c[1][0] = s * A->c[1][0] + t * B->c[1][0];
    M->c[1][1] = s * A->c[1][1] + t * B->c[1][1];
  }

void r2x2_mul(r2x2_t *A, r2x2_t *B, r2x2_t *M)
  { r2x2_t RR;
    RR.c[0][0] = A->c[0][0]*B->c[0][0] + A->c[0][1]*B->c[1][0];
    RR.c[0][1] = A->c[0][0]*B->c[0][1] + A->c[0][1]*B->c[1][1];
    RR.c[1][0] = A->c[1][0]*B->c[0][0] + A->c[1][1]*B->c[1][0];
    RR.c[1][1] = A->c[1][0]*B->c[0][1] + A->c[1][1]*B->c[1][1];
    (*M) = RR;
  }

void r2x2_mul_tr(r2x2_t *A, r2x2_t *B, r2x2_t *M)
  {
    r2x2_t RR;
    RR.c[0][0] = A->c[0][0]*B->c[0][0] + A->c[0][1]*B->c[0][1];
    RR.c[0][1] = A->c[0][0]*B->c[1][0] + A->c[0][1]*B->c[1][1];
    RR.c[1][0] = A->c[1][0]*B->c[0][0] + A->c[1][1]*B->c[0][1];
    RR.c[1][1] = A->c[1][0]*B->c[1][0] + A->c[1][1]*B->c[1][1];
    (*M) = RR;
  }

void r2x2_tr_mul(r2x2_t *A, r2x2_t *B, r2x2_t *M)
  {
    r2x2_t RR;
    RR.c[0][0] = A->c[0][0]*B->c[0][0] + A->c[1][0]*B->c[1][0];
    RR.c[0][1] = A->c[0][0]*B->c[0][1] + A->c[1][0]*B->c[1][1];
    RR.c[1][0] = A->c[0][1]*B->c[0][0] + A->c[1][1]*B->c[1][0];
    RR.c[1][1] = A->c[0][1]*B->c[0][1] + A->c[1][1]*B->c[1][1];
    (*M) = RR;
  }

double r2x2_det(r2x2_t *A)
  { return A->c[0][0] * A->c[1][1] - A->c[1][0] * A->c[0][1]; }

void r2x2_adj(r2x2_t *A, r2x2_t *M)
  { r2x2_t RR;

    RR.c[0][0] =   A->c[1][1];
    RR.c[0][1] = - A->c[0][1];
    RR.c[1][0] = - A->c[1][0];
    RR.c[1][1] =   A->c[0][0];
    
    (*M) = RR;
  }

void r2x2_inv(r2x2_t *A, r2x2_t *M)
  { double d = A->c[0][0] * A->c[1][1] - A->c[1][0] * A->c[0][1];
    r2x2_t RR;

    RR.c[0][0] =   A->c[1][1]/d;
    RR.c[0][1] = - A->c[0][1]/d;
    RR.c[1][0] = - A->c[1][0]/d;
    RR.c[1][1] =   A->c[0][0]/d;
    
    (*M) = RR;
  }

void r2x2_sqrt(r2x2_t *A, r2x2_t *M)
 {
    double A00 = A->c[0][0];
    double A01 = A->c[0][1];
    double A10 = A->c[1][0];
    double A11 = A->c[1][1];
    double mx = fmax(fmax(fabs(A00),fabs(A01)),fmax(fabs(A10),fabs(A11)));
    double tr = A00 + A11; /* Trace. */
    double dt = A00*A11 - A01*A10;
    demand (dt >= 0.0, "determinant is negative\n");
    double sd = sqrt(dt);
    double R2 = tr + 2*sd;
    demand (R2 > 1.0e-15*mx, "denominator is not real\n");
    double R = sqrt(R2);
    M->c[0][0] = (A00 + sd)/R;
    M->c[0][1] = A01/R;
    M->c[1][0] = A10/R;
    M->c[1][1] = (A11 + sd)/R;
 }

double r2x2_norm_sqr(r2x2_t *A)
  {
    double A00 = A->c[0][0];
    double A01 = A->c[0][1];

    double A10 = A->c[1][0];
    double A11 = A->c[1][1];

    return A00*A00 + A01*A01 + A10*A10 + A11*A11; 
  }

double r2x2_norm(r2x2_t *A)
  {
    double A00 = A->c[0][0];
    double A01 = A->c[0][1];

    double A10 = A->c[1][0];
    double A11 = A->c[1][1];

    return sqrt( A00*A00 + A01*A01 + A10*A10 + A11*A11 ); 
  }


double r2x2_normalize(r2x2_t *A)
  { 
    double w = r2x2_norm(A);
    if (w != 0)
      { A->c[0][0] /= w; 
        A->c[0][1] /= w; 
        A->c[1][0] /= w; 
        A->c[1][1] /= w;
      }
    else
      { A->c[0][0] = NAN; 
        A->c[0][1] = NAN; 
        A->c[1][0] = NAN; 
        A->c[1][1] = NAN;
      } 
    return w;
  }

double r2x2_mod_norm_sqr(r2x2_t *A)
  {
    double D00 = A->c[0][0] - 1;
    double D01 = A->c[0][1];

    double D10 = A->c[1][0];
    double D11 = A->c[1][1] - 1;

    return D00*D00 + D01*D01 + D10*D10 + D11*D11; 
  }

double r2x2_L_inf_norm(r2x2_t *A)
  {
    double mv = 0.0;
    mv = fmax(mv, fabs(A->c[0][0]));
    mv = fmax(mv, fabs(A->c[0][1]));
    
    mv = fmax(mv, fabs(A->c[1][0]));
    mv = fmax(mv, fabs(A->c[1][1]));
    
    return mv;
  }

double r2x2_L_inf_normalize(r2x2_t *A)
  {
    double mv = r2x2_L_inf_norm(A);
    A->c[0][0] /= mv;
    A->c[0][1] /= mv;
              
    A->c[1][0] /= mv;
    A->c[1][1] /= mv;
    
    return mv;
  }

void r2x2_diff_sqr(r2x2_t *A, r2x2_t *B, r2x2_t *R, double *dabs2P, double *drel2P)
  {
    double dabs2 = 0.0;
    double drel2 = 0.0;
    for (uint32_t i = 0;  i < N; i++) 
      { for (uint32_t j = 0;  j < N; j++) 
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

void r2x2_from_rows(r2_t *a, r2_t *b, r2x2_t *M)
  { M->c[0][0] = a->c[0];
    M->c[1][0] = b->c[0];
    
    M->c[0][1] = a->c[1];
    M->c[1][1] = b->c[1];
  }

void r2x2_from_cols(r2_t *a, r2_t *b, r2x2_t *M)
  { M->c[0][0] = a->c[0];
    M->c[1][0] = a->c[1];
           
    M->c[0][1] = b->c[0];
    M->c[1][1] = b->c[1];
  }

bool_t r2x2_is_unif_scaling(r2x2_t *M, double s)
  {
    for (uint32_t i = 0;  i < N; i++)
      for (uint32_t j = 0;  j < N; j++)
        { if (M->c[i][j] != (i == j ? s : 0.0)) { return FALSE; } }
    return TRUE;
  }

void r2x2_rot90(r2x2_t *M)
  { M->c[0][0] = 00.0;
    M->c[0][1] = +1.0;
    M->c[1][0] = -1.0;
    M->c[1][1] = 00.0;
  }

void r2x2_sym_eigen(r2x2_t *A, r2_t *e, r2x2_t *M)
  {
    /* Grab the elements of {A}: */
    double A00 = A->c[0][0];
    double A01 = A->c[0][1];
    double A10 = A->c[1][0];
    double A11 = A->c[1][1];
    /* Symmetrize the matrix: */
    double Y = 0.5*A01 + 0.5*A10; /* Off-diagonal element. */
    if (Y == 0)
      { /* Matrix is diagonal -- eigenvalues are {A00} and {A11}: */
        if (A00 >= A11)
          { /* Largest eigenvalue is {A00}, or the two are equal: */
            e->c[0] = A00; e->c[1] = A11;
            /* The corresponding eigenvectors are {(+1,00)} and {(00,+1)}: */
            if (M != NULL) { r2x2_ident(M); }
          }
        else
          { /* Largest eigenvalue is {A11}: */
            e->c[0] = A11; e->c[1] = A00;
            /* The corresponding eigenvectors are {(00,+1)} and {(-1,00)}: */
            if (M != NULL) { r2x2_rot90(M); }
          }
      }
    else
      { /* Get intermediate quantities {T,X}: */
        double T = 0.5*A00 + 0.5*A11; /* Half-trace of matrix. */
        double X = 0.5*A00 - 0.5*A11; /* Half-unbalance of matrix. */
        /* Eigenvalues are {T ± sqrt(X^2 + Y^2)}: */
        double Q = hypot(X, Y); /* Square root of discriminant. */
        e->c[0] = T + Q; /* Max eigenvalue. */
        e->c[1] = T - Q; /* Min eigenvalue. */
        /* Compute the eigenvectors if so requested: */
        if (M != NULL)
          { /* assert(Q != 0); */
            /* The eigenvectors are the square roots of {±(X/Q + Y/Q*IMAG)}. */
            double x = X/Q, y = Y/Q;
            double c, s;  /* The eigenvector of eigenvalue {T+Q} is {±(c,s)}. */
            /* Split into two cases for numerical stability: */
            if (x >= 0)
              { c = sqrt(0.5*(1 + x)); s = 0.5*y/c; }
            else
              { s = sqrt(0.5*(1 - x)); c = 0.5*y/s; }
            /* Store into array with proper signs: */
            /* Fix the sign so that the arg is in {-PI} (exc.) to {+PI} (inc.): */
            if ((c < 0) || ((c == 0) && (s < 0))) { c = -c; s = -s; }
            M->c[0][0] = +c;
            M->c[0][1] = +s;
            M->c[1][0] = -s;
            M->c[1][1] = +c;
          }
      }
  }

void r2x2_rot_and_scale(r2_t *p, r2x2_t *M)
  {
    M->c[0][0] = p->c[0];
    M->c[0][1] = p->c[1];
    M->c[1][0] = - p->c[1];
    M->c[1][1] = p->c[0];
  }

void r2x2_from_point_pairs(r2_vec_t *p1, r2_t *bar1, r2_vec_t *p2, r2_t *bar2, r2x2_t *M)
  {
    uint32_t np = p1->ne;
    assert(np == p2->ne);
    double debug = FALSE;
    
    if (debug) { fprintf(stderr, "--- computing the linear matrix ---\n"); }
    
    r2x2_ident(M);
    /* Compute mean 2×2 linear transformation {S} from {p1-bar1} to {p2-bar2} by least squares: */
    r2x2_t A; r2x2_zero(&A); /* Moment matrix. */
    r2x2_t B; r2x2_zero(&B); /* Projection matrix. */
    for (uint32_t k = 0;  k < np; k++)
      { 
        /* Reduce points {p1.e[k],p2.e[k]} relative to barycenter: */
        r2_t q1k = p1->e[k];
        r2_t q2k = p2->e[k];
        if (bar1 != NULL) { r2_sub(&q1k, bar1, &q1k); }
        if (bar2 != NULL) { r2_sub(&q2k, bar2, &q2k); }
        /* Accumulate moments and projections: */
        for (int32_t i = 0; i < 2; i ++)
          { for (uint32_t j = 0;  j < 2; j++)
              { A.c[i][j] += q1k.c[i]*q1k.c[j];
                B.c[i][j] += q1k.c[i]*q2k.c[j];
              }
          }
      }
    r2x2_t Z; r2x2_inv(&A, &Z);
    r2x2_mul(&Z, &B, M);

    if (debug) 
      { fprintf(stderr, "  linear matrix:\n");
        r2x2_gen_print(stderr, M, "%13.6e", "", "\n", "\n", "    [ ", " ", " ]");
      }
  }

void r2x2_print(FILE *f, r2x2_t *A)
  { r2x2_gen_print(f, A, NULL, NULL, NULL, NULL, NULL, NULL, NULL); }

void r2x2_gen_print 
  ( FILE *f, r2x2_t *A,
    char *fmt, 
    char *olp, char *osep, char *orp,
    char *ilp, char *isep, char *irp
  )
  { rmxn_gen_print(f, N, N, &(A->c[0][0]), fmt, olp, osep, orp, ilp, isep, irp); }  
