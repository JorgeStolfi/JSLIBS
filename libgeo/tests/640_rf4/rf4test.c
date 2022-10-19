/* rf4test --- test program for rf4.h, rf4x4.h  */
/* Last edited on 2022-01-04 08:43:21 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>
#include <jsrandom.h>
#include <flt.h>

#include <rf4.h>
#include <rf4x4.h>
#include <rn_test_tools.h>
#include <rfn_check.h>

#define N 4
#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_rf4(int32_t verbose);
void test_rf4x4(int32_t verbose);
void throw_matrix(rf4x4_t *m);

int32_t main (int32_t argc, char **argv)
  { int32_t i;
    srand(1993);
    srandom(1993);

    for (i = 0; i < 100; i++) test_rf4(i < 3);
    for (i = 0; i < 100; i++) test_rf4x4(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_rf4(int32_t verbose)
  { rf4_t a, b, c, d, e, para, perp;
    double r, s, t;
    float sf, tf;
    double rr, ss, tt;
    double mag;
    int32_t i, j, k;

    if (verbose)
      { fprintf(stderr,
          "sizeof(rf4_t) = %lu  %d*sizeof(double) = %lu\n",
          sizeof(rf4_t), N, N*sizeof(double)
        );
      }


    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf4_zero ---\n"); }
    a = rf4_zero();
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i], 0.0, NO, NO, "rf4_zero error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf4_all ---\n"); }
    a = rf4_all(3.14f);
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i], 3.14f, NO, NO, "rf4_all error"); }

    if (verbose) { fprintf(stderr, "--- rf4_axis ---\n"); }
    for (k = 0; k < N; k++)
      { a = rf4_axis(k);
        for (i = 0; i < N; i++)
          { rn_check_eq(a.c[i],(i == k ? 1.0 : 0.0), NO, NO, "rf4_axis error"); }
      }

    if (verbose) { fprintf(stderr, "--- rf4_throw_cube ---\n"); }
    a = rf4_throw_cube();
    for (i = 0; i < N; i++)
      { affirm(a.c[i] != a.c[(i+1)%N], "rf4_throw probable error(1)"); 
        /* Check whether there are more than 8 nonzero bits: */
        double vv = a.c[i]*256.0;
        affirm(vv != floor(vv), "rf4_throw error(3)"); 
        affirm((a.c[i] > -1.0) && (a.c[i] < 1.0), "rf4_throw error(2)"); 
      }

    if (verbose) { fprintf(stderr, "--- rf4_throw_dir ---\n"); }
    /* Should check uniformity... */
    a = rf4_throw_dir();
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "rf4_throw_dir error(1)"); }
    /* Check whether the norm is 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    rn_check_eps(1,rr,0.0000001 * rr, NO, NO, "rf4_throw_dir error (2)");

    if (verbose) { fprintf(stderr, "--- rf4_throw_ball ---\n"); }
    /* Should check uniformity... */
    a = rf4_throw_ball();
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "rf4_throw_ball error(1)"); }
    /* Check whether the norm is at most 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    demand(rr <= 1 + 0.0000001*rr, "rf4_throw_ball error (2)");

    if (verbose) { fprintf(stderr, "--- rf4_add ---\n"); }
    a = rf4_throw_cube();
    b = rf4_throw_cube();
    d = rf4_add(&a, &b);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] + b.c[i], NO, NO, "rf4_add error"); }

    if (verbose) { fprintf(stderr, "--- rf4_sub ---\n"); }
    a = rf4_throw_cube();
    b = rf4_throw_cube();
    d = rf4_sub(&a, &b);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] - b.c[i], NO, NO, "rf4_sub error"); }

    if (verbose) { fprintf(stderr, "--- rf4_neg ---\n"); }
    a = rf4_throw_cube();
    d = rf4_neg(&a);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],- a.c[i], NO, NO, "rf4_neg error"); }

    if (verbose) { fprintf(stderr, "--- rf4_scale ---\n"); }
    s = drandom();
    a = rf4_throw_cube();
    d = rf4_scale(s, &a);
    for (i = 0; i < N; i++)
      { float zi = (float)(s*a.c[i]);
        rn_check_eq(d.c[i], zi, NO, NO, "rf4_scale error(1)");
      }

    if (verbose) { fprintf(stderr, "--- rf4_mix ---\n"); }
    s = drandom();
    t = drandom();
    a = rf4_throw_cube();
    b = rf4_throw_cube();
    d = rf4_mix(s, &a, t, &b);
    for (i = 0; i < N; i++)
      { float ddi = (float)(s * a.c[i] + t * b.c[i]);
        rn_check_eq(d.c[i],ddi, NO, NO, "rf4_mix error");
      }

    if (verbose) { fprintf(stderr, "--- rf4_mix_in ---\n"); }
    s = drandom();
    a = rf4_throw_cube();
    b = rf4_throw_cube();
    d = b;
    rf4_mix_in(s, &a, &d);
    for (i = 0; i < N; i++)
      { float ddi = (float)(b.c[i] + s * a.c[i]);
        rn_check_eq(d.c[i],ddi, NO, NO, "rf4_mix_in error");
      }

    if (verbose) { fprintf(stderr, "--- rf4_weigh ---\n"); }
    a = rf4_throw_cube();
    b = rf4_throw_cube();
    d = rf4_weigh(&a, &b);
    for (i = 0; i < N; i++)
      { float ddi = (float)(a.c[i] * b.c[i]);
        rn_check_eq(d.c[i],ddi, NO, NO, "rf4_weigh error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf4_unweigh ---\n"); }
    a = rf4_throw_cube();
    b = rf4_throw_cube();
    d = rf4_unweigh(&a, &b);
    for (i = 0; i < N; i++)
      { float ddi = (float)(((double)a.c[i]) / b.c[i]);
        rn_check_eq(d.c[i],ddi, NO, NO, "rf4_unweigh error");
      }

    if (verbose) { fprintf(stderr, "--- rf4_rot_axis ---\n"); }
    { a = rf4_throw_cube();
      int32_t i = int32_abrandom(0, N-1);
      int32_t j = int32_abrandom(0, N-2); if (j >= i) { j++; }
      double ang = 2.1*M_PI*drandom();
      d = rf4_rot_axis(&a, i, j, ang);
      rfn_check_rot_axis(N, a.c, i, j, ang, d.c, "rf4_rot_axis error");
    }

    if (verbose) { fprintf(stderr, "--- rf4_norm, rf4_norm_sqr, rf4_L_inf_norm ---\n"); }
    a = rf4_throw_cube();
    r = rf4_norm(&a);
    s = rf4_norm_sqr(&a);
    t = rf4_L_inf_norm(&a);
    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < N; i++)
      { float ai = fabsf(a.c[i]);
        ss += ((double)ai)*ai; 
        if (ai > tt) { tt = ai; }
      }
    rr = sqrt(ss);
    rn_check_eps(r, rr, 0.0000001 * rr, NO, NO, "rf4_norm error");
    rn_check_eps(s, ss, 0.0000001 * ss, NO, NO, "rf4_norm_sqr error");
    rn_check_eq(t, tt, NO, NO, "rf4_L_inf_norm error");

    if (verbose) { fprintf(stderr, "--- rf4_dist, rf4_dist_sqr, rf4_L_inf_dist ---\n"); }
    a = rf4_throw_cube();
    b = rf4_throw_cube();
    r = rf4_dist(&a, &b);
    s = rf4_dist_sqr(&a, &b);
    t = rf4_L_inf_dist(&a, &b);
    ss = 0.0;
    tt = 0.0;

    for (i = 0; i < N; i++)
      { double di = fabs((double)(a.c[i]) - b.c[i]);
        ss += di*di; 
        float dfi = (float)di;
        if (dfi > tt) { tt = dfi; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.0000001 * rr, NO, NO, "rf4_dist error");
    rn_check_eps(s,ss,0.0000001 * ss, NO, NO, "rf4_dist_sqr error");
    rn_check_eq(t,tt, NO, NO, "rf4_L_inf_dist error");

    if (verbose) { fprintf(stderr, "--- rf4_dir, rf4_L_inf_dir ---\n"); }
    a = rf4_throw_cube();
    b = rf4_dir(&a, &r);
    d = rf4_L_inf_dir(&a, &sf);
    ss = rf4_norm(&a);
    tf = rf4_L_inf_norm(&a);
    for (i = 0; i < N; i++)
      { float bi = (float)(a.c[i]/ss);
        rn_check_eps(b.c[i], bi, 0.0000001 * ss, NO, NO, "rf4_dir error");
        float di = (float)((double)(a.c[i])/tf);
        rn_check_eps(d.c[i], di, 0.0000001 * tf, NO, NO, "rf4_L_inf_dir error");
      }

    if (verbose) { fprintf(stderr, "--- rf4_dot, rf4_cos, rf4_sin, rf4_angle ---\n"); }
    a = rf4_throw_cube();
    b = rf4_throw_cube();
    r = rf4_dot(&a, &b);
    double S = rf4_sin(&a, &b);
    double C = rf4_cos(&a, &b);
    double A = rf4_angle(&a, &b);
    mag = sqrt(rf4_dot(&a,&a)*rf4_dot(&b,&b));
    rr = 0.0;
    for (i = 0; i < N; i++) { rr += a.c[i]*b.c[i]; }
    double CC = rr/(rf4_norm(&a)*rf4_norm(&b));
    rn_check_eps(r,rr,0.0000001 * mag, NO, NO, "rf4_dot error(1)");
    rn_check_eps(C,CC,0.0000001, NO, NO, "rf4_cos error(1)");
    d = a;
    rf4_mix_in(-rr/rf4_norm_sqr(&b), &b, &d);
    double SS = rf4_norm(&d)/rf4_norm(&a);
    rn_check_eps(S,SS,0.0000001, NO, NO, "rf4_sin error(1)");
    double AA = atan2(SS, CC);
    rn_check_eps(A,AA,0.0000001, NO, NO, "rf4_angle error(1)");
    for (i = 0; i < N; i++)
      { a = rf4_axis(i);
        for (j = 0; j < N; j++)
          { b = rf4_axis(j);
            r = rf4_dot(&a, &b);
            s = rf4_sin(&a, &b);
            t = rf4_cos(&a, &b);
            rr = (i == j ? 1.0 : 0.0);
            rn_check_eq(r,rr, NO, NO, "rf4_dot error(2)");
            rn_check_eq(t,rr, NO, NO, "rf4_dot error(2)");
            rn_check_eq(s,1.0 - rr, NO, NO, "rf4_dot error(2)");
          }
      }

    if (verbose) { fprintf(stderr, "--- rf4_cross ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        int32_t i2 = (i + 2) % N;
        int32_t i3 = (i + 3) % N;
        float sign = ((i % 2) == 0 ? 1.0 : -1.0);
        int32_t p;
        a = rf4_axis(i0);
        b = rf4_axis(i1);
        c = rf4_axis(i2);
        d = rf4_cross(&a, &b, &c);
        e = rf4_axis(i3);
        for (p = 0; p < N; p++)
          { float ep = sign*e.c[p];
            rn_check_eq(d.c[p],ep, NO, NO, "rf4_cross error(x)");
          }
      }
    /* Test on random vectors: */
    a = rf4_throw_cube();
    b = rf4_throw_cube();
    c = rf4_throw_cube();
    d = rf4_cross(&a, &b, &c);
    mag = rf4_norm(&a)*rf4_norm(&b)*rf4_norm(&c);
    r = rf4_dot(&a, &d);
    rn_check_eps(r,0.0,0.0000001 * mag*rf4_norm(&a), NO, NO, "rf4_cross error(1)");
    r = rf4_dot(&b, &d);
    rn_check_eps(r,0.0,0.0000001 * mag*rf4_norm(&b), NO, NO, "rf4_cross error(2)");
    r = rf4_dot(&c, &d);
    rn_check_eps(r,0.0,0.0000001 * mag*rf4_norm(&c), NO, NO, "rf4_cross error(3)");
 
    if (verbose) { fprintf(stderr, "--- rf4_det ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        int32_t i2 = (i + 2) % N;
        int32_t i3 = (i + 3) % N;
        float sign = ((i % 2) == 0 ? 1.0 : -1.0);
        a = rf4_axis(i0);
        b = rf4_axis(i1);
        c = rf4_axis(i2);
        d = rf4_axis(i3);
        r = rf4_det(&a, &b, &c, &d);
        rn_check_eq(r, sign, NO, NO, "rf4_det error(2)");
      }
    /* Test on random vectors: */
    a = rf4_throw_cube();
    b = rf4_throw_cube();
    c = rf4_throw_cube();
    d = rf4_throw_cube();
    r = rf4_det(&a, &b, &c, &d);
    e = rf4_cross(&a, &b, &c);
    rr = rf4_dot(&d, &e);
    mag = rf4_norm(&a)*rf4_norm(&b)*rf4_norm(&c)*rf4_norm(&d);
    rn_check_eps(r,rr,0.0000001 * mag, NO, NO, "rf4_det error(1)");

    if (verbose) { fprintf(stderr, "--- rf4_decomp ---\n"); }
    a = rf4_throw_cube();
    b = rf4_throw_cube();
    r = rf4_decomp(&a, &b, &para, &perp);
    rr = rf4_dot(&a, &b)/rf4_norm_sqr(&b);  
    rn_check_eps(r,rr,0.0000001 * (fabs(r) + fabs(rr)), NO, NO, "rf4_decomp error(1)");
    c = rf4_add(&para, &perp);
    s = rf4_dist(&a, &c);
    affirm (s <= 0.0000001 * rf4_norm(&a), "rf4_decomp error(2)");
    s = rf4_dot(&perp, &b);
    rn_check_eps(s,0.0,0.0000001 * rf4_norm(&b), NO, NO, "rf4_decomp error(3)");
    t = rf4_dot(&para, &perp);
    rn_check_eps(t,0.0,0.0000001 * rf4_norm(&a), NO, NO, "rf4_decomp error(4)");

    if (verbose) { fprintf(stderr, "--- rf4_print ---\n"); }
    if (verbose)
      { a = rf4_throw_cube ();
        fprintf(stderr, "a = ");
        rf4_print(stderr, &a);
        fputc('\n', stderr);
      }
  }

void test_rf4x4(int32_t verbose)
  { rf4x4_t A, B, C;
    rf4_t a, b, c, bb, cc;
    double r, s, t;
    float rf, sf;
    double rr, ss, tt;
    double mag;
    int32_t i, j, k;

    if (verbose) { fprintf(stderr, "--- Size and allocation ---\n"); }
    if (verbose)
      { fprintf(stderr,
          "sizeof(rf4x4_t) = %lu  %d*%d*sizeof(double) = %lu\n",
          sizeof(rf4x4_t), N, N, N*N*sizeof(double)
        );
        fprintf(stderr, "&B = %016lx\n", (long unsigned)&B);
        fprintf(stderr, "&A-&B = %lu\n", ((long unsigned)(&A))-((long unsigned)(&B)));
        fprintf(stderr, "&B-&C = %lu\n", ((long unsigned)(&B))-((long unsigned)(&C)));
        fprintf(stderr, "&(B.c) = %016lx\n", (long unsigned)&(B.c));
        fprintf(stderr, "B.c = %016lx\n", (long unsigned)(B.c));
        fprintf(stderr, "&(B.c[0]) = %016lx\n", (long unsigned)&(B.c[0]));
        fprintf(stderr, "B.c[0] = %016lx\n", (long unsigned)(B.c[0]));
        fprintf(stderr, "&(B.c[0][0]) = %016lx\n", (long unsigned)&(B.c[0][0]));
      }

    if (verbose) { fprintf(stderr, "--- Indexing and addressing ---\n"); }
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { float *Aij = &(A.c[i][j]); 
          affirm(Aij == ((float *)&(A.c))+(N*i)+j, "rf4x4_t indexing error");
        }

    if (verbose) { fprintf(stderr, "--- rf4x4_zero, rf4x4_ident ---\n"); }
    A = rf4x4_zero();
    B = rf4x4_ident();
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(A.c[i][j],0.0, NO, NO, "rf4x4_zero error"); 
            rn_check_eq(B.c[i][j],(i == j ? 1.0 : 0.0), NO, NO, "rf4x4_ident error");
          }
      }

    if (verbose) { fprintf(stderr, "--- rf4x4_get_row, rf4x4_set_row, rf4x4_get_col, rf4x4_set_col ---\n"); }
    throw_matrix(&A);
    int32_t dir; /* 0 for row, 1 for col. */
    for (dir = 0; dir < 2; dir++)
      { for (i = 0; i < N; i++)
          { /* Check {rf4x4_get_row,rf4x4_get_col}: */
            a = rf4_throw_cube();
            if (dir == 0) { a = rf4x4_get_row(&A, i); } else { a = rf4x4_get_col(&A, i); }
            for (j = 0; j < N; j++)
              { float vj = (dir == 0 ? A.c[i][j] : A.c[j][i]);
                affirm(vj = a.c[j], "rf4x4_get_row/rf4x4_get_col error");
              }
            /* Check {rf4x4_set_row,rf4x4_set_col}: */
            a = rf4_throw_cube();
            if (dir == 0) { rf4x4_set_row(&A, i, &a); } else { rf4x4_set_col(&A, i, &a); }
            for (j = 0; j < N; j++)
              { float vj = (dir == 0 ? A.c[i][j] : A.c[j][i]);
                affirm(vj = a.c[j], "rf4x4_set_row/rf4x4_set_col error");
              }
          }
      }

    if (verbose) { fprintf(stderr, "--- rf4x4_map_row, rf4x4_map_col ---\n"); }
    throw_matrix(&A);
    a = rf4_throw_cube();
    b = rf4x4_map_row(&a, &A);
    c = rf4x4_map_col(&A, &a);
    for (i = 0; i < N; i++)
      { double sb = 0;
        double sc = 0;
        for (j = 0; j < N; j++)
          { sb += (double)(a.c[j]) * A.c[j][i];
            sc += (double)(A.c[i][j]) * a.c[j];
          }
        bb.c[i] = (float)sb;
        cc.c[i] = (float)sc;
      }
    r = rf4_dist(&b, &bb);
    affirm(r < 0.0000001 * rf4_norm(&bb), "rf4_map_row error");
    s = rf4_dist(&c, &cc);
    affirm(s < 0.0000001 * rf4_norm(&cc), "rf4_map_col error");

    if (verbose) { fprintf(stderr, "--- rf4x4_scale ---\n"); }
    throw_matrix(&A);
    r = drandom();
    C = rf4x4_scale(r, &A);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { float sel = (float)(r * A.c[i][j]);
            rn_check_eps(C.c[i][j], sel, 0.0000001 * fabs(sel), NO, NO,
              "rf4x4_scale error"
            );
          }
      }

    if (verbose) { fprintf(stderr, "--- rf4x4_mul ---\n"); }
    throw_matrix(&A);
    throw_matrix(&B);
    C = rf4x4_mul(&A, &B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sum = 0.0;
            for (k = 0; k < N; k++) { sum += (double)(A.c[i][k])*B.c[k][j]; }
            rn_check_eps(C.c[i][j],sum,0.0000001 * fabs(sum), NO, NO,
              "rf4x4_mul error"
            );
          }
      }

    if (verbose) { fprintf(stderr, "--- rf4x4_mul_tr ---\n"); }
    throw_matrix(&A);
    throw_matrix(&B);
    C = rf4x4_mul_tr(&A, &B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sum = 0.0;
            for (k = 0; k < N; k++) { sum += (double)(A.c[i][k])*B.c[j][k]; }
            rn_check_eps(C.c[i][j],sum,0.0000001 * fabs(sum), NO, NO,
              "rf4x4_mul error"
            );
          }
      }

    if (verbose) { fprintf(stderr, "--- rf4x4_transp ---\n"); }
    throw_matrix(&A);
    B = rf4x4_transp(&A);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(B.c[i][j], A.c[j][i], NO, NO, "rf4x4_transp error (1)"); }
      }
    /* In-place transpose: */
    B = A;
    B = rf4x4_transp(&B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(B.c[i][j], A.c[j][i], NO, NO, "rf4x4_transp error (2)"); }
      }

    if (verbose) { fprintf(stderr, "--- rf4x4_det ---\n"); }
    throw_matrix(&A);
    for (i = 0; i < N; i++)
      { int32_t k = (i + 1) % N;
        for (j = 0; j < N; j++)
          { /* Check for linearity */
            rf = frandom();
            A.c[i][j] = rf;
            rr = rf4x4_det(&A);

            sf = frandom();
            A.c[i][j] = sf;
            ss = rf4x4_det(&A);

            t = drandom();
            A.c[i][j] = (float)(rf*(1-t) + sf*t);
            tt = rf4x4_det(&A);
            mag = fabs(rr) + fabs(ss) + fabs(tt);
            rn_check_eps(tt,(rr*(1 - t) + ss*t),000000001 * mag, NO, NO,
              "rf4x4_det error(1)"
            );
          }

        /* Row swap test: */
        r = rf4x4_det(&A);
        for (j = 0; j < N; j++)
          { float tf = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = tf; }
        rr = rf4x4_det(&A);
        mag = fabs(r) + fabs(rr);
        rn_check_eps(r,-rr,000000001 * mag, NO, NO, "rf4x4_det error(2)");

        /* Col swap test: */
        r = rf4x4_det(&A);
        for (j = 0; j < N; j++)
          { float tf = A.c[j][i]; A.c[j][i] = A.c[j][k]; A.c[j][k] = tf; }
        rr = rf4x4_det(&A);
        mag = fabs(r) + fabs(rr);
        rn_check_eps(r,-rr,000000001 * mag, NO, NO, "rf4x4_det error(3)");
      }

    if (verbose) { fprintf(stderr, "--- rf4x4_inv ---\n"); }
    throw_matrix(&A);
    B = rf4x4_inv(&A);
    C = rf4x4_mul(&A, &B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double val = (i == j ? 1.0 : 0.0);
            affirm((C.c[i][j] - val) < 000000001, "rf4x4_inv error");
          }
      }

    if (verbose) { fprintf(stderr, "--- rf4x4_norm,rf4x4_norm_sqr,rf4x4_mod_norm ---\n"); }
    throw_matrix(&A);
    s = rf4x4_norm_sqr(&A);
    r = rf4x4_norm(&A);
    t = rf4x4_mod_norm_sqr(&A);
    ss = 0; tt = 0;
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double Aij = A.c[i][j];
            double Dij = (i == j ? Aij - 1 : Aij);
            ss += Aij*Aij;
            tt += Dij*Dij;
          }
      }
    affirm(ss >= 0, "rf4x4_norm_sqr error");
    affirm(fabs(ss - s) < 000000001, "rf4x4_norm_sqr error");
    rr = sqrt(ss);
    affirm(fabs(rr - r) < 000000001, "rf4x4_norm error");
    affirm(tt >= 0, "rf4x4_mod_norm_sqr error");
    affirm(fabs(tt - t) < 000000001, "rf4x4_mod_norm_sqr error");
 
    if (verbose) { fprintf(stderr, "--- rf4x4_print ---\n"); }
    if (verbose)
      { throw_matrix (&A);
        fprintf(stderr, "A = ");
        rf4x4_print(stderr, &A);
        fputc('\n', stderr);
      }
  }  

void throw_matrix(rf4x4_t *m)
  { int32_t i, j;
    rf4_t a;
    for (i = 0; i < N; i++)
      { a = rf4_throw_cube();
        for (j = 0; j < N; j++) { m->c[i][j] = a.c[j]; }
      }
  }
    
