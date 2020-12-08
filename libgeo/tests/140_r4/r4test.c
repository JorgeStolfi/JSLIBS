/* r4test --- test program for r4.h, r4x4.h  */
/* Last edited on 2016-04-03 13:09:22 by stolfilocal */

#include <r4.h>
#include <r4_extra.h>
#include <r4x4.h>
#include <rn_test_tools.h>

#include <affirm.h>
#include <jsrandom.h>
#include <flt.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 4
#define NO NULL

/* Internal prototypes */

int main (int argc, char **argv);
void test_r4(int verbose);
void test_r4x4(int verbose);
void throw_matrix(r4x4_t *m);

int main (int argc, char **argv)
  { int i;
    srand(1993);
    srandom(1993);

    for (i = 0; i < 100; i++) test_r4(i < 3);
    for (i = 0; i < 100; i++) test_r4x4(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_r4(int verbose)
  { r4_t a, b, c, d, e, para, perp;
    double r, s, t;
    double rr, ss, tt;
    double mag;
    int i, j, k;

    if (verbose)
      { fprintf(stderr,
          "sizeof(r4_t) = %lud  %d*sizeof(double) = %lud\n",
          sizeof(r4_t), N, N*sizeof(double)
        );
      }

    if (verbose) { fprintf(stderr, "--- r4_zero ---\n"); }
    r4_zero(&a);
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i],0.0, NO, NO, "r4_zero error"); }

    if (verbose) { fprintf(stderr, "--- r4_all ---\n"); }
    r4_all(3.14, &a);
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i],3.14, NO, NO, "r4_all error"); }

    if (verbose) { fprintf(stderr, "--- r4_axis ---\n"); }
    for (k = 0; k < N; k++)
      { r4_axis(k, &a);
        for (i = 0; i < N; i++)
          { rn_check_eq(a.c[i],(i == k ? 1.0 : 0.0), NO, NO, "r4_axis error"); }
      }

    if (verbose) { fprintf(stderr, "--- r4_throw_cube ---\n"); }
    r4_throw_cube(&a);
    for (i = 0; i < N; i++)
      { affirm(a.c[i] != a.c[(i+1)%N], "r4_throw probable error(1)"); 
        /* Check whether there are more than 8 nonzero bits: */
        double vv = a.c[i]*256.0;
        affirm(vv != floor(vv), "r4_throw error(3)"); 
        affirm((a.c[i] > -1.0) && (a.c[i] < 1.0), "r4_throw error(2)"); 
      }

    if (verbose) { fprintf(stderr, "--- r4_throw_dir ---\n"); }
    /* Should check uniformity... */
    r4_throw_dir(&a);
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "r4_throw_dir error(1)"); }
    /* Check whether the norm is 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    rn_check_eps(1,rr,0.000000001 * rr, NO, NO, "r4_throw_dir error (2)");

    if (verbose) { fprintf(stderr, "--- r4_throw_ball ---\n"); }
    /* Should check uniformity... */
    r4_throw_ball(&a);
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "r4_throw_ball error(1)"); }
    /* Check whether the norm is at most 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    demand(rr <= 1 + 0.000000001*rr, "r4_throw_ball error (2)");

    if (verbose) { fprintf(stderr, "--- r4_add ---\n"); }
    r4_throw_cube(&a);
    r4_throw_cube(&b);
    r4_add(&a, &b, &d);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] + b.c[i], NO, NO, "r4_add error"); }

    if (verbose) { fprintf(stderr, "--- r4_sub ---\n"); }
    r4_throw_cube(&a);
    r4_throw_cube(&b);
    r4_sub(&a, &b, &d);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] - b.c[i], NO, NO, "r4_sub error"); }

    if (verbose) { fprintf(stderr, "--- r4_neg ---\n"); }
    r4_throw_cube(&a);
    r4_neg(&a, &d);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],- a.c[i], NO, NO, "r4_neg error"); }

    if (verbose) { fprintf(stderr, "--- r4_scale ---\n"); }
    s = drandom();
    r4_throw_cube(&a);
    r4_scale(s, &a, &d);
    for (i = 0; i < N; i++)
      { double zi = s*a.c[i];
        rn_check_eq(d.c[i],zi, NO, NO, "r4_scale error(1)");
      }

    if (verbose) { fprintf(stderr, "--- r4_mix ---\n"); }
    s = drandom();
    t = drandom();
    r4_throw_cube(&a);
    r4_throw_cube(&b);
    r4_mix(s, &a, t, &b, &d);
    for (i = 0; i < N; i++)
      { double ddi = s * a.c[i] + t * b.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r4_mix error");
      }

    if (verbose) { fprintf(stderr, "--- r4_mix_in ---\n"); }
    s = drandom();
    r4_throw_cube(&a);
    r4_throw_cube(&b);
    d = b;
    r4_mix_in(s, &a, &d);
    for (i = 0; i < N; i++)
      { double ddi = b.c[i] + s * a.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r4_mix_in error");
      }

    if (verbose) { fprintf(stderr, "--- r4_weigh ---\n"); }
    r4_throw_cube(&a);
    r4_throw_cube(&b);
    r4_weigh(&a, &b, &d);
    for (i = 0; i < N; i++)
      { double ddi = a.c[i] * b.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r4_weigh error");
      }

    if (verbose) { fprintf(stderr, "--- r4_norm, r4_norm_sqr, r4_L_inf_norm ---\n"); }
    r4_throw_cube(&a);
    r = r4_norm(&a);
    s = r4_norm_sqr(&a);
    t = r4_L_inf_norm(&a);
    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < N; i++)
      { double ai = fabs(a.c[i]);
        ss += ai*ai; 
        if (ai > tt) { tt = ai; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.000000001 * rr, NO, NO, "r4_norm error");
    rn_check_eps(s,ss,0.000000001 * ss, NO, NO, "r4_norm_sqr error");
    rn_check_eq(t,tt, NO, NO, "r4_L_inf_norm error");

    if (verbose) { fprintf(stderr, "--- r4_dist, r4_dist_sqr, r4_L_inf_dist ---\n"); }
    r4_throw_cube(&a);
    r4_throw_cube(&b);
    r = r4_dist(&a, &b);
    s = r4_dist_sqr(&a, &b);
    t = r4_L_inf_dist(&a, &b);
    ss = 0.0;
    tt = 0.0;

    for (i = 0; i < N; i++)
      { double di = fabs(a.c[i] - b.c[i]);
        ss += di*di; 
        if (di > tt) { tt = di; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.000000001 * rr, NO, NO, "r4_dist error");
    rn_check_eps(s,ss,0.000000001 * ss, NO, NO, "r4_dist_sqr error");
    rn_check_eq(t,tt, NO, NO, "r4_L_inf_dist error");

    if (verbose) { fprintf(stderr, "--- r4_dir, r4_L_inf_dir ---\n"); }
    r4_throw_cube(&a);
    r = r4_dir(&a, &b);
    s = r4_L_inf_dir(&a, &d);
    ss = r4_norm(&a);
    tt = r4_L_inf_norm(&a);
    for (i = 0; i < N; i++)
      { rn_check_eps(b.c[i],a.c[i]/ss,0.000000001 * ss, NO, NO, "r4_dir error");
        rn_check_eps(d.c[i],a.c[i]/tt,0.000000001 * tt, NO, NO, "r4_L_inf_dir error");
      }

    if (verbose) { fprintf(stderr, "--- r4_dot, r4_cos, r4_sin, r4_angle ---\n"); }
    r4_throw_cube(&a);
    r4_throw_cube(&b);
    r = r4_dot(&a, &b);
    double S = r4_sin(&a, &b);
    double C = r4_cos(&a, &b);
    double A = r4_angle(&a, &b);
    mag = sqrt(r4_dot(&a,&a)*r4_dot(&b,&b));
    rr = 0.0;
    for (i = 0; i < N; i++) { rr += a.c[i]*b.c[i]; }
    double CC = rr/(r4_norm(&a)*r4_norm(&b));
    rn_check_eps(r,rr,0.000000001 * mag, NO, NO, "r4_dot error(1)");
    rn_check_eps(C,CC,0.000000001, NO, NO, "r4_cos error(1)");
    d = a;
    r4_mix_in(-rr/r4_norm_sqr(&b), &b, &d);
    double SS = r4_norm(&d)/r4_norm(&a);
    rn_check_eps(S,SS,0.000000001, NO, NO, "r4_sin error(1)");
    double AA = atan2(SS, CC);
    rn_check_eps(A,AA,0.000000001, NO, NO, "r4_angle error(1)");
    for (i = 0; i < N; i++)
      { r4_axis(i, &a);
        for (j = 0; j < N; j++)
          { r4_axis(j, &b);
            r = r4_dot(&a, &b);
            s = r4_sin(&a, &b);
            t = r4_cos(&a, &b);
            rr = (i == j ? 1.0 : 0.0);
            rn_check_eq(r,rr, NO, NO, "r4_dot error(2)");
            rn_check_eq(t,rr, NO, NO, "r4_dot error(2)");
            rn_check_eq(s,1.0 - rr, NO, NO, "r4_dot error(2)");
          }
      }

    if (verbose) { fprintf(stderr, "--- r4_cross ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int i0 = (i + 0) % N;
        int i1 = (i + 1) % N;
        int i2 = (i + 2) % N;
        int i3 = (i + 3) % N;
        double sign = ((i % 2) == 0 ? 1.0 : -1.0);
        int p;
        r4_axis(i0, &a);
        r4_axis(i1, &b);
        r4_axis(i2, &c);
        r4_cross(&a, &b, &c, &d);
        r4_axis(i3, &e);
        for (p = 0; p < N; p++)
          { double ep = sign*e.c[p];
            rn_check_eq(d.c[p],ep, NO, NO, "r4_cross error(x)");
          }
      }
    /* Test on random vectors: */
    r4_throw_cube(&a);
    r4_throw_cube(&b);
    r4_throw_cube(&c);
    r4_cross(&a, &b, &c, &d);
    mag = r4_norm(&a)*r4_norm(&b)*r4_norm(&c);
    r = r4_dot(&a, &d);
    rn_check_eps(r,0.0,0.00000001 * mag*r4_norm(&a), NO, NO, "r4_cross error(1)");
    r = r4_dot(&b, &d);
    rn_check_eps(r,0.0,0.00000001 * mag*r4_norm(&b), NO, NO, "r4_cross error(2)");
    r = r4_dot(&c, &d);
    rn_check_eps(r,0.0,0.00000001 * mag*r4_norm(&c), NO, NO, "r4_cross error(3)");
 
    if (verbose) { fprintf(stderr, "--- r4_pick_ortho ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int i0 = (i + 0) % N;
        r4_axis(i0, &a);
        double ma = r4_L_inf_norm(&a);
        r4_pick_ortho(&a, &d);
        double mo = r4_L_inf_norm(&d);
        rn_check_eq(ma, mo, NO, NO, "r4_pick_ortho error(0)");
        r = r4_dot(&a, &d);
        rn_check_eps(r, 0.0, 0.00000001*ma, NO, NO, "r4_pick_ortho error(1)");
      }
    /* Test on random vectors: */
    { r4_throw_cube(&a);
      double ma = r4_L_inf_norm(&a);
      r4_pick_ortho(&a, &d);
      double mo = r4_L_inf_norm(&d); 
      rn_check_eq(ma, mo, NO, NO, "r4_pick_ortho error(2)");
      r = r4_dot(&a, &d);
      rn_check_eps(r, 0.0, 0.00000001*ma, NO, NO, "r4_pick_ortho error(3)");
    }
    
    if (verbose) { fprintf(stderr, "--- r4_det ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int i0 = (i + 0) % N;
        int i1 = (i + 1) % N;
        int i2 = (i + 2) % N;
        int i3 = (i + 3) % N;
        double sign = ((i % 2) == 0 ? 1.0 : -1.0);
        r4_axis(i0, &a);
        r4_axis(i1, &b);
        r4_axis(i2, &c);
        r4_axis(i3, &d);
        r = r4_det(&a, &b, &c, &d);
        rn_check_eq(r,sign, NO, NO, "r4_det error(2)");
      }
    /* Test on random vectors: */
    r4_throw_cube(&a);
    r4_throw_cube(&b);
    r4_throw_cube(&c);
    r4_throw_cube(&d);
    r = r4_det(&a, &b, &c, &d);
    r4_cross(&a, &b, &c, &e);
    rr = r4_dot(&d, &e);
    mag = r4_norm(&a)*r4_norm(&b)*r4_norm(&c)*r4_norm(&d);
    rn_check_eps(r,rr,0.00000001 * mag, NO, NO, "r4_det error(1)");

    if (verbose) { fprintf(stderr, "--- r4_decomp ---\n"); }
    r4_throw_cube(&a);
    r4_throw_cube(&b);
    r = r4_decomp(&a, &b, &para, &perp);
    rr = r4_dot(&a, &b)/r4_norm_sqr(&b);  
    rn_check_eps(r,rr,0.000000001 * (fabs(r) + fabs(rr)), NO, NO, "r4_decomp error(1)");
    r4_add(&para, &perp, &c);
    s = r4_dist(&a, &c);
    affirm (s <= 0.000000001 * r4_norm(&a), "r4_decomp error(2)");
    s = r4_dot(&perp, &b);
    rn_check_eps(s,0.0,0.000000001 * r4_norm(&b), NO, NO, "r4_decomp error(3)");
    t = r4_dot(&para, &perp);
    rn_check_eps(t,0.0,0.000000001 * r4_norm(&a), NO, NO, "r4_decomp error(4)");

    if (verbose) { fprintf(stderr, "--- r4_print ---\n"); }
    if (verbose)
      { r4_throw_cube (&a);
        fprintf(stderr, "a = ");
        r4_print(stderr, &a);
        fputc('\n', stderr);
      }
  }

void test_r4x4(int verbose)
  { r4x4_t A, B, C;
    r4_t a, b, c, bb, cc;
    double r, s, t;
    double rr, ss, tt;
    double mag;
    int i, j, k;

    if (verbose) { fprintf(stderr, "--- Size and allocation ---\n"); }
    if (verbose)
      { fprintf(stderr,
          "sizeof(r4x4_t) = %lud  %d*%d*sizeof(double) = %lud\n",
          sizeof(r4x4_t), N, N, N*N*sizeof(double)
        );
        fprintf(stderr, "&B = %016lx\n", (long unsigned)&B);
        fprintf(stderr, "&A-&B = %lud\n", ((long unsigned)(&A))-((long unsigned)(&B)));
        fprintf(stderr, "&B-&C = %lud\n", ((long unsigned)(&B))-((long unsigned)(&C)));
        fprintf(stderr, "&(B.c) = %016lx\n", (long unsigned)&(B.c));
        fprintf(stderr, "B.c = %016lx\n", (long unsigned)(B.c));
        fprintf(stderr, "&(B.c[0]) = %016lx\n", (long unsigned)&(B.c[0]));
        fprintf(stderr, "B.c[0] = %016lx\n", (long unsigned)(B.c[0]));
        fprintf(stderr, "&(B.c[0][0]) = %016lx\n", (long unsigned)&(B.c[0][0]));
      }

    if (verbose) { fprintf(stderr, "--- Indexing and addressing ---\n"); }
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { double *Aij = &(A.c[i][j]); 
          affirm(Aij = ((double *)&A)+(N*i)+j, "r4x4_t indexing error");
        }

    if (verbose) { fprintf(stderr, "--- r4x4_zero, r4x4_ident ---\n"); }
    r4x4_zero(&A);
    r4x4_ident(&B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(A.c[i][j],0.0, NO, NO, "r4x4_zero error"); 
            rn_check_eq(B.c[i][j],(i == j ? 1.0 : 0.0), NO, NO, "r4x4_ident error");
          }
      }

    if (verbose) { fprintf(stderr, "--- r4x4_get_row, r4x4_set_row, r4x4_get_col, r4x4_set_col ---\n"); }
    throw_matrix(&A);
    int dir; /* 0 for row, 1 for col. */
    for (dir = 0; dir < 2; dir++)
      { for (i = 0; i < N; i++)
          { /* Check {r4x4_get_row,r4x4_get_col}: */
            r4_throw_cube(&a);
            if (dir == 0) { r4x4_get_row(&A, i, &a); } else { r4x4_get_col(&A, i, &a); }
            for (j = 0; j < N; j++)
              { double vj = (dir == 0 ? A.c[i][j] : A.c[j][i]);
                affirm(vj = a.c[j], "r4x4_get_row/r4x4_get_col error");
              }
            /* Check {r4x4_set_row,r4x4_set_col}: */
            r4_throw_cube(&a);
            if (dir == 0) { r4x4_set_row(&A, i, &a); } else { r4x4_set_col(&A, i, &a); }
            for (j = 0; j < N; j++)
              { double vj = (dir == 0 ? A.c[i][j] : A.c[j][i]);
                affirm(vj = a.c[j], "r4x4_set_row/r4x4_set_col error");
              }
          }
      }

    if (verbose) { fprintf(stderr, "--- r4x4_map_row, r4x4_map_col ---\n"); }
    throw_matrix(&A);
    r4_throw_cube(&a);
    r4x4_map_row(&a, &A, &b);
    r4x4_map_col(&A, &a, &c);
    r4_zero(&bb);
    r4_zero(&cc);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { bb.c[j] += a.c[i] * A.c[i][j];
            cc.c[i] += A.c[i][j] * a.c[j];
          }
      }
    r = r4_dist(&b, &bb);
    affirm(r < 0.000000001 * r4_norm(&bb), "r4_map_row error");
    s = r4_dist(&c, &cc);
    affirm(s < 0.000000001 * r4_norm(&cc), "r4_map_col error");

    if (verbose) { fprintf(stderr, "--- r4x4_scale ---\n"); }
    throw_matrix(&A);
    r = drandom();
    r4x4_scale(r, &A, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sel = r * A.c[i][j];
            rn_check_eps(C.c[i][j],sel,0.000000001 * fabs(sel), NO, NO,
              "r4x4_scale error"
            );
          }
      }

    if (verbose) { fprintf(stderr, "--- r4x4_mul ---\n"); }
    throw_matrix(&A);
    throw_matrix(&B);
    r4x4_mul(&A, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sum = 0.0;
            for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[k][j]; }
            rn_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
              "r4x4_mul error"
            );
          }
      }

    if (verbose) { fprintf(stderr, "--- r4x4_mul_tr ---\n"); }
    throw_matrix(&A);
    throw_matrix(&B);
    r4x4_mul_tr(&A, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sum = 0.0;
            for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[j][k]; }
            rn_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
              "r4x4_mul error"
            );
          }
      }

    if (verbose) { fprintf(stderr, "--- r4x4_transp ---\n"); }
    throw_matrix(&A);
    r4x4_transp(&A, &B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(B.c[i][j],A.c[j][i], NO, NO, "r4x4_transp error (1)"); }
      }
    /* In-place transpose: */
    B = A;
    r4x4_transp(&B, &B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(B.c[i][j],A.c[j][i], NO, NO, "r4x4_transp error (2)"); }
      }

    if (verbose) { fprintf(stderr, "--- r4x4_det ---\n"); }
    throw_matrix(&A);
    for (i = 0; i < N; i++)
      { int k = (i + 1) % N;
        for (j = 0; j < N; j++)
          { /* Check for linearity */
            r = drandom();
            A.c[i][j] = r;
            rr = r4x4_det(&A);

            s = drandom();
            A.c[i][j] = s;
            ss = r4x4_det(&A);

            t = drandom();
            A.c[i][j] = r*(1-t) + s*t;
            tt = r4x4_det(&A);
            mag = fabs(rr) + fabs(ss) + fabs(tt);
            rn_check_eps(tt,(rr*(1 - t) + ss*t),000000001 * mag, NO, NO,
              "r4x4_det error(1)"
            );
          }

        /* Row swap test: */
        r = r4x4_det(&A);
        for (j = 0; j < N; j++)
          { double t = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = t; }
        rr = r4x4_det(&A);
        mag = fabs(r) + fabs(rr);
        rn_check_eps(r,-rr,000000001 * mag, NO, NO, "r4x4_det error(2)");

        /* Col swap test: */
        r = r4x4_det(&A);
        for (j = 0; j < N; j++)
          { double t = A.c[j][i]; A.c[j][i] = A.c[j][k]; A.c[j][k] = t; }
        rr = r4x4_det(&A);
        mag = fabs(r) + fabs(rr);
        rn_check_eps(r,-rr,000000001 * mag, NO, NO, "r4x4_det error(3)");
      }

    if (verbose) { fprintf(stderr, "--- r4x4_inv ---\n"); }
    throw_matrix(&A);
    r4x4_inv(&A, &B);
    r4x4_mul(&A, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double val = (i == j ? 1.0 : 0.0);
            affirm((C.c[i][j] - val) < 000000001, "r4x4_inv error");
          }
      }

    if (verbose) { fprintf(stderr, "--- r4x4_norm,r4x4_norm_sqr,r4x4_mod_norm ---\n"); }
    throw_matrix(&A);
    s = r4x4_norm_sqr(&A);
    r = r4x4_norm(&A);
    t = r4x4_mod_norm_sqr(&A);
    ss = 0; tt = 0;
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double Aij = A.c[i][j];
            double Dij = (i == j ? Aij - 1 : Aij);
            ss += Aij*Aij;
            tt += Dij*Dij;
          }
      }
    affirm(ss >= 0, "r4x4_norm_sqr error");
    affirm(fabs(ss - s) < 000000001, "r4x4_norm_sqr error");
    rr = sqrt(ss);
    affirm(fabs(rr - r) < 000000001, "r4x4_norm error");
    affirm(tt >= 0, "r4x4_mod_norm_sqr error");
    affirm(fabs(tt - t) < 000000001, "r4x4_mod_norm_sqr error");
 
    if (verbose) { fprintf(stderr, "--- r4x4_print ---\n"); }
    if (verbose)
      { throw_matrix (&A);
        fprintf(stderr, "A = ");
        r4x4_print(stderr, &A);
        fputc('\n', stderr);
      }
  }  

void throw_matrix(r4x4_t *m)
  { int i, j;
    r4_t a;
    for (i = 0; i < N; i++)
      { r4_throw_cube(&a);
        for (j = 0; j < N; j++) { m->c[i][j] = a.c[j]; }
      }
  }
    
