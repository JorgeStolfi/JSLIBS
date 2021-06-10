/* See r2_aff_map.h */
/* Last edited on 2021-06-09 19:49:49 by jstolfi */ 

/* Based on R2_AFF_MAP.m3, created 1994-05-04 by J. Stolfi. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <r2.h>
#include <r2x2.h>
#include <affirm.h>
#include <sign.h>
#include <sign_get.h>

#include <r2_aff_map.h>

#define NC 2
  /* Number of Cartesian coordinates in a point. */

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

double *r2_aff_map_elem_addr(r2_aff_map_t *A, int32_t s)
  { demand((s >= 0) && (s < 6), "invalid elem index");
    if (s < 4) 
      { int32_t i = s / 2, j = s % 2;
        return &(A->mat.c[i][j]);
      }
    else
      { int32_t j = s - 4;
        return &(A->disp.c[j]);
      }
  }

void r2_aff_map_apply(r2_t *p, r2_aff_map_t *A, r2_t *q)
  {
    r2x2_map_row(p, &(A->mat), q);
    r2_add(q, &(A->disp), q);
  }

void r2_aff_map_invert(r2_aff_map_t *A, r2_aff_map_t *B)
  {
    r2x2_inv(&(A->mat), &(B->mat));
    r2_t td;
    r2x2_map_row(&(A->disp), &(B->mat), &td);
    r2_neg(&td, &(B->disp));
  }

void r2_aff_map_compose(r2_aff_map_t *A, r2_aff_map_t *B, r2_aff_map_t *C)
  {
    r2_t td;
    r2x2_map_row(&(A->disp), &(B->mat), &td);
    r2_add(&td, &(B->disp), &(C->disp));
    r2x2_mul(&(A->mat), &(B->mat), &(C->mat));
  }

r2_aff_map_t r2_aff_map_translation(r2_t *vec)
  { r2_aff_map_t A;
    r2x2_ident(&(A.mat));

    A.disp.c[0] = vec->c[0];
    A.disp.c[1] = vec->c[1];
    
    return A;
  }

r2_aff_map_t r2_aff_map_rot_scale(double ang, double scale)
  {
    double c = cos(ang);
    double s = sin(ang);

    r2_aff_map_t A;
    A.mat.c[0][0] = +c*scale;
    A.mat.c[0][1] = +s*scale;
    A.mat.c[1][0] = -s*scale;
    A.mat.c[1][1] = +c*scale;

    A.disp.c[0] = 0.0;
    A.disp.c[1] = 0.0;
    
    return A;
  }

void r2_aff_map_disp_sqr
  ( r2_aff_map_t *A, 
    r2_aff_map_t *B, 
    r2_aff_map_t *R,
    double *dabs2P,
    double *drel2P
  )
  {
    double dabs2 = 0.0;
    double drel2 = 0.0;
    for (int32_t s = 0; s < 6; s++) 
      { double *Rp = r2_aff_map_elem_addr(R, s);
        if ((*Rp) != 0.0)
          { double Re = (*Rp);
            double *Ap = r2_aff_map_elem_addr(A, s); double Ae = (*Ap);
            double *Bp = r2_aff_map_elem_addr(B, s); double Be = (*Bp);
            double d = Ae - Be;
            dabs2 += d*d;
            drel2 += (d/Re)*(d/Re);
          }
      }
    if (dabs2P != NULL) { (*dabs2P) = dabs2; }
    if (drel2P != NULL) { (*drel2P) = drel2; }
  }
    
double r2_aff_map_mismatch_sqr(r2_aff_map_t *A, r2_aff_map_t *B)
  {
    /* Hope the math is right: */
    r2x2_t M, H; r2_t v;
    r2x2_sub(&(A->mat), &(B->mat), &M);
    r2_sub(&(A->disp), &(B->disp), &v);
    r2x2_mul_tr(&M, &M, &H);
    double h2 = (H.c[0][0] + H.c[1][1])/2;
    double d2 = r2_dot(&v, &v);
    return h2 + d2;
  }

r2_aff_map_t r2_aff_map_from_points(r2_t *o, r2_t *p, r2_t *q)
  {
    r2_aff_map_t A; /* The resulting map. */
    
    A.disp = (*o);
    A.mat.c[0][0] = p->c[0] - o->c[0];
    A.mat.c[0][1] = p->c[1] - o->c[1];
    A.mat.c[1][0] = q->c[0] - o->c[0];
    A.mat.c[1][1] = q->c[1] - o->c[1];

    return A;
  }

r2_aff_map_t r2_aff_map_from_point_pairs(int32_t np, r2_t p1[], r2_t p2[], double w[])
  {
    int32_t k;
    double debug = FALSE;
    
    if (debug) { fprintf(stderr, "--- computing the affine matrix ---\n"); }
    
    r2_aff_map_t A;
    if (np == 0)
      { /* Identity map: */
        r2x2_ident(&(A.mat));
        A.disp = (r2_t){{ 0,0 }};
      }
    else
      { /* Compute the barycenters of {p1} and {p2}: */
        r2_t bar1; r2_barycenter(np, p1, w, &bar1);
        if (debug) { r2_gen_print(stderr, &bar1, "%10.4f", "  bar1 = ( ", " ", " )\n"); }
        demand(r2_is_finite(&bar1), "barycenter of {p1} is undefined");

        r2_t bar2; r2_barycenter(np, p2, w, &bar2);
        if (debug) { r2_gen_print(stderr, &bar2, "%10.4f", "  bar2 = ( ", " ", " )\n"); }
        demand(r2_is_finite(&bar2), "barycenter of {p2} is undefined");

        /* Determine the linear map matrix {A.mat}: */
        if (np == 1)
          { /* Translation: */
            r2x2_ident(&(A.mat));
            r2_sub(&(bar2), &bar1, &(A.disp));
          }
        else if (np == 2)
          { /* Compute rotation & scale matrix: */
            r2_t q1; r2_sub(&(p1[0]), &bar1, &q1);
            r2x2_t R1; r2x2_rot_and_scale(&q1, &R1);

            r2_t q2; r2_sub(&(p2[0]), &bar2, &q2);
            r2x2_t R2; r2x2_rot_and_scale(&q2, &R2);

            r2x2_inv(&R1, &R1); r2x2_mul(&R1, &R2, &(A.mat));
          }
        else
          { /* Compute general 2×2 linear transformation by least squares: */

            r2x2_t M; r2x2_zero(&M); /* Moment matrix. */
            r2x2_t P; r2x2_zero(&P); /* Projection matrix. */
            for (k = 0; k < np; k++)
              { 
                double wk = (w == NULL ? 1.0 : w[k]);
                /* Reduce points relative to barycenter: */
                r2_t q1k, q2k;
                r2_sub(&(p1[k]), &bar1, &q1k);
                r2_sub(&(p2[k]), &bar2, &q2k);
                /* Accumulate moments and projections: */
                int32_t i,j;
                for(i = 0; i < 2; i ++)
                  { for (j = 0; j < 2; j++)
                      { M.c[i][j] += wk*q1k.c[i]*q1k.c[j];
                        P.c[i][j] += wk*q1k.c[i]*q2k.c[j];
                      }
                  }
              }
            r2x2_t Z; r2x2_inv(&M, &Z);
            r2x2_mul(&Z, &P, &(A.mat));
          }
        if (debug) 
          { fprintf(stderr, "  linear matrix:\n");
            r2x2_gen_print(stderr, &(A.mat), "%13.6e", "", "\n", "\n", "    [ ", " ", " ]");
          }

        /* Compute the displacement vector {A.disp}: */
        r2_t v; r2x2_map_row(&bar1, &(A.mat), &v);
        r2_sub(&bar2, &v, &(A.disp));
      }

    return A;
  }

void r2_aff_map_print (FILE *f, r2_aff_map_t *A)
  { r2_aff_map_gen_print(f, A, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); }

void r2_aff_map_gen_print 
  ( FILE *f, r2_aff_map_t *A,
    char *mfmt, char *dfmt, 
    char *olp, char *osep, char *orp,
    char *ilp, char *isep, char *irp,
    char *dsep
  )
  { if (dsep == NULL) { dsep = "\n"; }
    r2x2_gen_print(f, &(A->mat), mfmt, olp, osep, orp, ilp, isep, irp);
    fputs(dsep, f);
    r2_gen_print(f, &(A->disp), dfmt, ilp, isep, irp);
  }  
        
