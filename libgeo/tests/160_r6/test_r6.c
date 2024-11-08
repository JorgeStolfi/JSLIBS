/* test_r6 --- test program for r6.h, r6x6.h  */
/* Last edited on 2024-11-08 12:44:44 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>
#include <jsrandom.h>

#include <r6.h>
#include <r6_extra.h>
#include <r6x6.h>
#include <rn_test_tools.h>

#define N 6
#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_r6(int32_t verbose);
void test_r6x6(int32_t verbose);
void throw_matrix(r6x6_t *m);

int32_t main (int32_t argc, char **argv)
  { int32_t i;
    srand(1993);
    srandom(1993);

    for (i = 0; i < 100; i++) test_r6(i < 3);
    for (i = 0; i < 100; i++) test_r6x6(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_r6(int32_t verbose)
  { r6_t a, b, c, d, e, f, g, para, perp;
    double r, s, t;
    double rr, ss, tt;
    double mag;
    int32_t i, j, k;

    /* ---------------------------------------------------------------------- */
    if (verbose)
      { fprintf(stderr,
          "sizeof(r6_t) = %lu  %d*sizeof(double) = %lu\n",
          sizeof(r6_t), N, N*sizeof(double)
        );
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_zero ---\n"); }
    r6_zero(&a);
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i],0.0, NO, NO, "r6_zero error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_all ---\n"); }
    r6_all(3.14, &a);
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i],3.14, NO, NO, "r6_all error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_axis ---\n"); }
    for (k = 0; k < N; k++)
      { r6_axis(k, &a);
        for (i = 0; i < N; i++)
          { rn_check_eq(a.c[i],(i == k ? 1.0 : 0.0), NO, NO, "r6_axis error"); }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_throw_cube ---\n"); }
    r6_throw_cube(&a);
    for (i = 0; i < N; i++)
      { affirm(a.c[i] != a.c[(i+1)%N], "r6_throw probable error(1)"); 
        /* Check whether there are more than 8 nonzero bits: */
        double vv = a.c[i]*256.0;
        affirm(vv != floor(vv), "r6_throw error(3)"); 
        affirm((a.c[i] > -1.0) && (a.c[i] < 1.0), "r6_throw error(2)"); 
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_throw_dir ---\n"); }
    /* Should check uniformity... */
    r6_throw_dir(&a);
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "r6_throw_dir error(1)"); }
    /* Check whether the norm is 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    rn_check_eps(1,rr,0.000000001 * rr, NO, NO, "r6_throw_dir error (2)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_throw_ball ---\n"); }
    /* Should check uniformity... */
    r6_throw_ball(&a);
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "r6_throw_ball error(1)"); }
    /* Check whether the norm is at most 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    demand(rr <= 1 + 0.000000001*rr, "r6_throw_ball error (2)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_add ---\n"); }
    r6_throw_cube(&a);
    r6_throw_cube(&b);
    r6_add(&a, &b, &d);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] + b.c[i], NO, NO, "r6_add error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_sub ---\n"); }
    r6_throw_cube(&a);
    r6_throw_cube(&b);
    r6_sub(&a, &b, &d);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] - b.c[i], NO, NO, "r6_sub error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_neg ---\n"); }
    r6_throw_cube(&a);
    r6_neg(&a, &d);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],- a.c[i], NO, NO, "r6_neg error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_scale ---\n"); }
    s = drandom();
    r6_throw_cube(&a);
    r6_scale(s, &a, &d);
    for (i = 0; i < N; i++)
      { double zi = s*a.c[i];
        rn_check_eq(d.c[i],zi, NO, NO, "r6_scale error(1)");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_mix ---\n"); }
    s = drandom();
    t = drandom();
    r6_throw_cube(&a);
    r6_throw_cube(&b);
    r6_mix(s, &a, t, &b, &d);
    for (i = 0; i < N; i++)
      { double ddi = s * a.c[i] + t * b.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r6_mix error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_mix_in ---\n"); }
    s = drandom();
    r6_throw_cube(&a);
    r6_throw_cube(&b);
    d = b;
    r6_mix_in(s, &a, &d);
    for (i = 0; i < N; i++)
      { double ddi = b.c[i] + s * a.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r6_mix_in error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_weigh ---\n"); }
    r6_throw_cube(&a);
    r6_throw_cube(&b);
    r6_weigh(&a, &b, &d);
    for (i = 0; i < N; i++)
      { double ddi = a.c[i] * b.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r6_weigh error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_unweigh ---\n"); }
    r6_throw_cube(&a);
    r6_throw_cube(&b);
    r6_unweigh(&a, &b, &d);
    for (i = 0; i < N; i++)
      { double ddi = a.c[i] / b.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r6_unweigh error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_rot_axis ---\n"); }
    { r6_throw_cube(&a);
      int32_t i = int32_abrandom(0, N-1);
      int32_t j = int32_abrandom(0, N-2); if (j >= i) { j++; }
      double ang = 2.1*M_PI*drandom();
      r6_rot_axis(&a, i, j, ang, &d);
      rn_test_rot_axis(N, a.c, i, j, ang, d.c, "r6_rot_axis error");
    }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_norm, r6_norm_sqr, r6_L_inf_norm ---\n"); }
    r6_throw_cube(&a);
    r = r6_norm(&a);
    s = r6_norm_sqr(&a);
    t = r6_L_inf_norm(&a);
    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < N; i++)
      { double ai = fabs(a.c[i]);
        ss += ai*ai; 
        if (ai > tt) { tt = ai; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.000000001 * rr, NO, NO, "r6_norm error");
    rn_check_eps(s,ss,0.000000001 * ss, NO, NO, "r6_norm_sqr error");
    rn_check_eq(t,tt, NO, NO, "r6_L_inf_norm error");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_dist, r6_dist_sqr, r6_L_inf_dist ---\n"); }
    r6_throw_cube(&a);
    r6_throw_cube(&b);
    r = r6_dist(&a, &b);
    s = r6_dist_sqr(&a, &b);
    t = r6_L_inf_dist(&a, &b);
    ss = 0.0;
    tt = 0.0;

    for (i = 0; i < N; i++)
      { double di = fabs(a.c[i] - b.c[i]);
        ss += di*di; 
        if (di > tt) { tt = di; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.000000001 * rr, NO, NO, "r6_dist error");
    rn_check_eps(s,ss,0.000000001 * ss, NO, NO, "r6_dist_sqr error");
    rn_check_eq(t,tt, NO, NO, "r6_L_inf_dist error");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_dir, r6_L_inf_dir ---\n"); }
    r6_throw_cube(&a);
    r = r6_dir(&a, &b);
    s = r6_L_inf_dir(&a, &d);
    ss = r6_norm(&a);
    tt = r6_L_inf_norm(&a);
    for (i = 0; i < N; i++)
      { rn_check_eps(b.c[i],a.c[i]/ss,0.000000001 * ss, NO, NO, "r6_dir error");
        rn_check_eps(d.c[i],a.c[i]/tt,0.000000001 * tt, NO, NO, "r6_L_inf_dir error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_dot, r6_cos, r6_sin, r6_angle ---\n"); }
    r6_throw_cube(&a);
    r6_throw_cube(&b);
    r = r6_dot(&a, &b);
    double S = r6_sin(&a, &b);
    double C = r6_cos(&a, &b);
    double A = r6_angle(&a, &b);
    mag = sqrt(r6_dot(&a,&a)*r6_dot(&b,&b));
    rr = 0.0;
    for (i = 0; i < N; i++) { rr += a.c[i]*b.c[i]; }
    double CC = rr/(r6_norm(&a)*r6_norm(&b));
    rn_check_eps(r,rr,0.000000001 * mag, NO, NO, "r6_dot error(1)");
    rn_check_eps(C,CC,0.000000001, NO, NO, "r6_cos error(1)");
    d = a;
    r6_mix_in(-rr/r6_norm_sqr(&b), &b, &d);
    double SS = r6_norm(&d)/r6_norm(&a);
    rn_check_eps(S,SS,0.000000001, NO, NO, "r6_sin error(1)");
    double AA = atan2(SS, CC);
    rn_check_eps(A,AA,0.000000001, NO, NO, "r6_angle error(1)");
    for (i = 0; i < N; i++)
      { r6_axis(i, &a);
        for (j = 0; j < N; j++)
          { r6_axis(j, &b);
            r = r6_dot(&a, &b);
            s = r6_sin(&a, &b);
            t = r6_cos(&a, &b);
            rr = (i == j ? 1.0 : 0.0);
            rn_check_eq(r,rr, NO, NO, "r6_dot error(2)");
            rn_check_eq(t,rr, NO, NO, "r6_dot error(2)");
            rn_check_eq(s,1.0 - rr, NO, NO, "r6_dot error(2)");
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_cross ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        int32_t i2 = (i + 2) % N;
        int32_t i3 = (i + 3) % N;
        int32_t i4 = (i + 4) % N;
        int32_t i5 = (i + 5) % N;
        double sign = ((i % 2) == 0 ? 1.0 : -1.0);
        int32_t p;
        r6_axis(i0, &a);
        r6_axis(i1, &b);
        r6_axis(i2, &c);
        r6_axis(i3, &d);
        r6_axis(i4, &e);
        r6_cross(&a, &b, &c, &d, &e, &f);
        r6_axis(i5, &g);
        for (p = 0; p < N; p++)
          { double gp = sign*g.c[p];
            rn_check_eq(f.c[p],gp, NO, NO, "r6_cross error(x)");
          }
      }
    /* Test on random vectors: */
    r6_throw_cube(&a);
    r6_throw_cube(&b);
    r6_throw_cube(&c);
    r6_throw_cube(&d);
    r6_throw_cube(&e);
    r6_cross(&a, &b, &c, &d, &e, &f);
    mag = r6_norm(&a)*r6_norm(&b)*r6_norm(&c)*r6_norm(&d)*r6_norm(&e);
    r = r6_dot(&a, &f);
    rn_check_eps(r,0.0,0.00000001 * mag*r6_norm(&a), NO, NO, "r6_cross error(1)");
    r = r6_dot(&b, &f);
    rn_check_eps(r,0.0,0.00000001 * mag*r6_norm(&b), NO, NO, "r6_cross error(2)");
    r = r6_dot(&c, &f);
    rn_check_eps(r,0.0,0.00000001 * mag*r6_norm(&c), NO, NO, "r6_cross error(3)");
    r = r6_dot(&d, &f);
    rn_check_eps(r,0.0,0.00000001 * mag*r6_norm(&d), NO, NO, "r6_cross error(3)");
    r = r6_dot(&e, &f);
    rn_check_eps(r,0.0,0.00000001 * mag*r6_norm(&e), NO, NO, "r6_cross error(3)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_pick_ortho ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        r6_axis(i0, &a);
        double ma = r6_L_inf_norm(&a);
        r6_pick_ortho(&a, &d);
        double mo = r6_L_inf_norm(&d);
        rn_check_eq(ma, mo, NO, NO, "r6_pick_ortho error(0)");
        r = r6_dot(&a, &d);
        rn_check_eps(r, 0.0, 0.00000001*ma, NO, NO, "r6_pick_ortho error(1)");
      }
    /* Test on random vectors: */
    { r6_throw_cube(&a);
      double ma = r6_L_inf_norm(&a);
      r6_pick_ortho(&a, &d);
      double mo = r6_L_inf_norm(&d); 
      rn_check_eq(ma, mo, NO, NO, "r6_pick_ortho error(2)");
      r = r6_dot(&a, &d);
      rn_check_eps(r, 0.0, 0.00000001*ma, NO, NO, "r6_pick_ortho error(3)");
    }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_det ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        int32_t i2 = (i + 2) % N;
        int32_t i3 = (i + 3) % N;
        int32_t i4 = (i + 4) % N;
        int32_t i5 = (i + 5) % N;
        double sign = ((i % 2) == 0 ? 1.0 : -1.0);
        r6_axis(i0, &a);
        r6_axis(i1, &b);
        r6_axis(i2, &c);
        r6_axis(i3, &d);
        r6_axis(i4, &e);
        r6_axis(i5, &f);
        r = r6_det(&a, &b, &c, &d, &e, &f);
        rn_check_eq(r,sign, NO, NO, "r6_det error(2)");
      }
    /* Test on random vectors: */
    r6_throw_cube(&a);
    r6_throw_cube(&b);
    r6_throw_cube(&c);
    r6_throw_cube(&d);
    r6_throw_cube(&e);
    r6_throw_cube(&f);
    r = r6_det(&a, &b, &c, &d, &e, &f);
    r6_cross(&a, &b, &c, &d, &e, &g);
    rr = r6_dot(&f, &g);
    mag = r6_norm(&a)*r6_norm(&b)*r6_norm(&c)*r6_norm(&d)*r6_norm(&e)*r6_norm(&f);
    rn_check_eps(r,rr,0.00000001 * mag, NO, NO, "r6_det error(1)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_decomp ---\n"); }
    r6_throw_cube(&a);
    r6_throw_cube(&b);
    r = r6_decomp(&a, &b, &para, &perp);
    rr = r6_dot(&a, &b)/r6_norm_sqr(&b);  
    rn_check_eps(r,rr,0.000000001 * (fabs(r) + fabs(rr)), NO, NO, "r6_decomp error(1)");
    r6_add(&para, &perp, &c);
    s = r6_dist(&a, &c);
    affirm (s <= 0.000000001 * r6_norm(&a), "r6_decomp error(2)");
    s = r6_dot(&perp, &b);
    rn_check_eps(s,0.0,0.000000001 * r6_norm(&b), NO, NO, "r6_decomp error(3)");
    t = r6_dot(&para, &perp);
    rn_check_eps(t,0.0,0.000000001 * r6_norm(&a), NO, NO, "r6_decomp error(4)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6_print ---\n"); }
    if (verbose)
      { r6_throw_cube (&a);
        fprintf(stderr, "a = ");
        r6_print(stderr, &a);
        fputc('\n', stderr);
      }
  }

void test_r6x6(int32_t verbose)
  { r6x6_t A, B, C;
    r6_t a, b, c, bb, cc;
    double r, s, t;
    double rr, ss, tt;
    double mag;
    int32_t i, j, k;

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- Size and allocation ---\n"); }
    if (verbose)
      { fprintf(stderr,
          "sizeof(r6x6_t) = %lu  %d*%d*sizeof(double) = %lu\n",
          sizeof(r6x6_t), N, N, N*N*sizeof(double)
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

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- Indexing and addressing ---\n"); }
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { double *Bij = &(B.c[i][j]); 
          affirm(Bij = ((double *)&B)+(N*i)+j, "r6x6_t indexing error");
        }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6x6_zero, r6x6_ident ---\n"); }
    r6x6_zero(&A);
    r6x6_ident(&B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(A.c[i][j],0.0, NO, NO, "r6x6_zero error"); 
            rn_check_eq(B.c[i][j],(i == j ? 1.0 : 0.0), NO, NO, "r6x6_ident error");
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6x6_get_row, r6x6_set_row, r6x6_get_col, r6x6_set_col ---\n"); }
    throw_matrix(&A);
    int32_t dir; /* 0 for row, 1 for col. */
    for (dir = 0; dir < 2; dir++)
      { for (i = 0; i < N; i++)
          { /* Check {r6x6_get_row,r6x6_get_col}: */
            r6_throw_cube(&a);
            if (dir == 0) { r6x6_get_row(&A, i, &a); } else { r6x6_get_col(&A, i, &a); }
            for (j = 0; j < N; j++)
              { double vj = (dir == 0 ? A.c[i][j] : A.c[j][i]);
                affirm(vj = a.c[j], "r6x6_get_row/r6x6_get_col error");
              }
            /* Check {r6x6_set_row,r6x6_set_col}: */
            r6_throw_cube(&a);
            if (dir == 0) { r6x6_set_row(&A, i, &a); } else { r6x6_set_col(&A, i, &a); }
            for (j = 0; j < N; j++)
              { double vj = (dir == 0 ? A.c[i][j] : A.c[j][i]);
                affirm(vj = a.c[j], "r6x6_set_row/r6x6_set_col error");
              }
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6x6_map_row, r6x6_map_col ---\n"); }
    throw_matrix(&A);
    r6_throw_cube(&a);
    r6x6_map_row(&a, &A, &b);
    r6x6_map_col(&A, &a, &c);
    r6_zero(&bb);
    r6_zero(&cc);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { bb.c[j] += a.c[i] * A.c[i][j];
            cc.c[i] += A.c[i][j] * a.c[j];
          }
      }
    r = r6_dist(&b, &bb);
    affirm(r < 0.000000001 * r6_norm(&bb), "r6_map_row error");
    s = r6_dist(&c, &cc);
    affirm(s < 0.000000001 * r6_norm(&cc), "r6_map_col error");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6x6_scale ---\n"); }
    throw_matrix(&A);
    r = drandom();
    r6x6_scale(r, &A, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sel = r * A.c[i][j];
            rn_check_eps(C.c[i][j],sel,0.000000001 * fabs(sel), NO, NO,
              "r6x6_scale error"
            );
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6x6_mul ---\n"); }
    throw_matrix(&A);
    throw_matrix(&B);
    r6x6_mul(&A, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sum = 0.0;
            for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[k][j]; }
            rn_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
              "r6x6_mul error"
            );
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6x6_mul_tr ---\n"); }
    throw_matrix(&A);
    throw_matrix(&B);
    r6x6_mul_tr(&A, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sum = 0.0;
            for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[j][k]; }
            rn_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
              "r6x6_mul error"
            );
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6x6_transp ---\n"); }
    throw_matrix(&A);
    r6x6_transp(&A, &B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(B.c[i][j],A.c[j][i], NO, NO, "r6x6_transp error (1)"); }
      }
    /* In-place transpose: */
    B = A;
    r6x6_transp(&B, &B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(B.c[i][j],A.c[j][i], NO, NO, "r6x6_transp error (2)"); }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6x6_det ---\n"); }
    throw_matrix(&A);
    for (i = 0; i < N; i++)
      { int32_t k = (i + 1) % N;
        for (j = 0; j < N; j++)
          { /* Check for linearity */
            r = drandom();
            A.c[i][j] = r;
            rr = r6x6_det(&A);

            s = drandom();
            A.c[i][j] = s;
            ss = r6x6_det(&A);

            t = drandom();
            A.c[i][j] = r*(1-t) + s*t;
            tt = r6x6_det(&A);
            mag = fabs(rr) + fabs(ss) + fabs(tt);
            rn_check_eps(tt,(rr*(1 - t) + ss*t),000000001 * mag, NO, NO,
              "r6x6_det error(1)"
            );
          }

        /* Row swap test: */
        r = r6x6_det(&A);
        for (j = 0; j < N; j++)
          { double t = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = t; }
        rr = r6x6_det(&A);
        mag = fabs(r) + fabs(rr);
        rn_check_eps(r,-rr,000000001 * mag, NO, NO, "r6x6_det error(2)");

        /* Col swap test: */
        r = r6x6_det(&A);
        for (j = 0; j < N; j++)
          { double t = A.c[j][i]; A.c[j][i] = A.c[j][k]; A.c[j][k] = t; }
        rr = r6x6_det(&A);
        mag = fabs(r) + fabs(rr);
        rn_check_eps(r,-rr,000000001 * mag, NO, NO, "r6x6_det error(3)");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6x6_inv ---\n"); }
    throw_matrix(&A);
    r6x6_inv(&A, &B);
    r6x6_mul(&A, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double val = (i == j ? 1.0 : 0.0);
            affirm((C.c[i][j] - val) < 000000001, "r6x6_inv error");
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6x6_norm,r6x6_norm_sqr,r6x6_mod_norm ---\n"); }
    throw_matrix(&A);
    s = r6x6_norm_sqr(&A);
    r = r6x6_norm(&A);
    t = r6x6_mod_norm_sqr(&A);
    ss = 0; tt = 0;
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double Aij = A.c[i][j];
            double Dij = (i == j ? Aij - 1 : Aij);
            ss += Aij*Aij;
            tt += Dij*Dij;
          }
      }
    affirm(ss >= 0, "r6x6_norm_sqr error");
    affirm(fabs(ss - s) < 000000001, "r6x6_norm_sqr error");
    rr = sqrt(ss);
    affirm(fabs(rr - r) < 000000001, "r6x6_norm error");
    affirm(tt >= 0, "r6x6_mod_norm_sqr error");
    affirm(fabs(tt - t) < 000000001, "r6x6_mod_norm_sqr error");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6x6_normalize ---\n"); }
    throw_matrix(&A);
    B = A;
    s = r6x6_norm(&B);
    ss = r6x6_normalize(&B);
    affirm(fabs(ss - s) < 000000001, "r6x6_normalize result error");
    t = r6x6_norm(&B);
    tt = 1.0;
    affirm(fabs(tt - t) < 000000001, "r6x6_normalize norm error");
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double Aij = A.c[i][j];
            double Bij = B.c[i][j];
            affirm(fabs(Bij*ss - Aij) < 000000001, "r6x6_normalize elem error");
          }
      }
    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r6x6_print ---\n"); }
    if (verbose)
      { throw_matrix (&A);
        fprintf(stderr, "A = ");
        r6x6_print(stderr, &A);
        fputc('\n', stderr);
      }
  }  

void throw_matrix(r6x6_t *m)
  { int32_t i, j;
    r6_t a;
    for (i = 0; i < N; i++)
      { r6_throw_cube(&a);
        for (j = 0; j < N; j++) { m->c[i][j] = a.c[j]; }
      }
  }

