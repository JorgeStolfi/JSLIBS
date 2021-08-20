/* rf3test --- test program for rf3.h, rf3x3.h  */
/* Last edited on 2021-08-20 17:24:34 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>
#include <jsrandom.h>
#include <flt.h>

#include <rf3.h>
#include <rn_test_tools.h>
#include <rfn_check.h>

#define N 3
#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_rf3(int32_t verbose);

int32_t main (int32_t argc, char **argv)
  { int32_t i;
    srand(1993);
    srandom(1993);

    for (i = 0; i < 100; i++) test_rf3(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_rf3(int32_t verbose)
  { rf3_t a, b, c, d, e, para, perp;
    double r, s, t;
    float sf, tf;
    double rr, ss, tt;
    double mag;
    int32_t i, j, k;

    if (verbose)
      { fprintf(stderr,
          "sizeof(rf3_t) = %lud  %d*sizeof(double) = %lud\n",
          sizeof(rf3_t), N, N*sizeof(double)
        );
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf3_zero ---\n"); }
    a = rf3_zero();
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i], 0.0, NO, NO, "rf3_zero error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf3_all ---\n"); }
    a = rf3_all(3.14f);
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i], 3.14f, NO, NO, "rf3_all error"); }

    if (verbose) { fprintf(stderr, "--- rf3_axis ---\n"); }
    for (k = 0; k < N; k++)
      { a = rf3_axis(k);
        for (i = 0; i < N; i++)
          { rn_check_eq(a.c[i],(i == k ? 1.0 : 0.0), NO, NO, "rf3_axis error"); }
      }

    if (verbose) { fprintf(stderr, "--- rf3_throw_cube ---\n"); }
    a = rf3_throw_cube();
    for (i = 0; i < N; i++)
      { affirm(a.c[i] != a.c[(i+1)%N], "rf3_throw probable error(1)"); 
        /* Check whether there are more than 8 nonzero bits: */
        double vv = a.c[i]*256.0;
        affirm(vv != floor(vv), "rf3_throw error(3)"); 
        affirm((a.c[i] > -1.0) && (a.c[i] < 1.0), "rf3_throw error(2)"); 
      }

    if (verbose) { fprintf(stderr, "--- rf3_throw_dir ---\n"); }
    /* Should check uniformity... */
    a = rf3_throw_dir();
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "rf3_throw_dir error(1)"); }
    /* Check whether the norm is 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    rn_check_eps(1,rr,0.0000001 * rr, NO, NO, "rf3_throw_dir error (2)");

    if (verbose) { fprintf(stderr, "--- rf3_throw_ball ---\n"); }
    /* Should check uniformity... */
    a = rf3_throw_ball();
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "rf3_throw_ball error(1)"); }
    /* Check whether the norm is at most 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    demand(rr <= 1 + 0.0000001*rr, "rf3_throw_ball error (2)");

    if (verbose) { fprintf(stderr, "--- rf3_add ---\n"); }
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    d = rf3_add(&a, &b);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] + b.c[i], NO, NO, "rf3_add error"); }

    if (verbose) { fprintf(stderr, "--- rf3_sub ---\n"); }
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    d = rf3_sub(&a, &b);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] - b.c[i], NO, NO, "rf3_sub error"); }

    if (verbose) { fprintf(stderr, "--- rf3_neg ---\n"); }
    a = rf3_throw_cube();
    d = rf3_neg(&a);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],- a.c[i], NO, NO, "rf3_neg error"); }

    if (verbose) { fprintf(stderr, "--- rf3_scale ---\n"); }
    s = drandom();
    a = rf3_throw_cube();
    d = rf3_scale(s, &a);
    for (i = 0; i < N; i++)
      { float zi = (float)(s*a.c[i]);
        rn_check_eq(d.c[i],zi, NO, NO, "rf3_scale error(1)");
      }

    if (verbose) { fprintf(stderr, "--- rf3_mix ---\n"); }
    s = drandom();
    t = drandom();
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    d = rf3_mix(s, &a, t, &b);
    for (i = 0; i < N; i++)
      { float ddi = (float)(s * a.c[i] + t * b.c[i]);
        rn_check_eq(d.c[i],ddi, NO, NO, "rf3_mix error");
      }

    if (verbose) { fprintf(stderr, "--- rf3_mix_in ---\n"); }
    s = drandom();
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    d = b;
    rf3_mix_in(s, &a, &d);
    for (i = 0; i < N; i++)
      { float ddi = (float)(b.c[i] + s * a.c[i]);
        rn_check_eq(d.c[i],ddi, NO, NO, "rf3_mix_in error");
      }

    if (verbose) { fprintf(stderr, "--- rf3_weigh ---\n"); }
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    d = rf3_weigh(&a, &b);
    for (i = 0; i < N; i++)
      { float ddi = (float)((double)(a.c[i]) * b.c[i]);
        rn_check_eq(d.c[i],ddi, NO, NO, "rf3_weigh error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf3_unweigh ---\n"); }
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    d = rf3_unweigh(&a, &b);
    for (i = 0; i < N; i++)
      { float ddi = (float)(((double)a.c[i]) / b.c[i]);
        rn_check_eq(d.c[i],ddi, NO, NO, "rf3_unweigh error");
      }

    if (verbose) { fprintf(stderr, "--- rf3_rot_axis ---\n"); }
    { a = rf3_throw_cube();
      int32_t i = int32_abrandom(0, N-1);
      int32_t j = int32_abrandom(0, N-2); if (j >= i) { j++; }
      double ang = 2.1*M_PI*drandom();
      d = rf3_rot_axis(&a, i, j, ang);
      rfn_check_rot_axis(N, a.c, i, j, ang, d.c, "rf3_rot_axis error");
    }

    if (verbose) { fprintf(stderr, "--- rf3_norm, rf3_norm_sqr, rf3_L_inf_norm ---\n"); }
    a = rf3_throw_cube();
    r = rf3_norm(&a);
    s = rf3_norm_sqr(&a);
    t = rf3_L_inf_norm(&a);
    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < N; i++)
      { double ai = fabs(a.c[i]);
        ss += ai*ai; 
        if (ai > tt) { tt = ai; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.0000001 * rr, NO, NO, "rf3_norm error");
    rn_check_eps(s,ss,0.0000001 * ss, NO, NO, "rf3_norm_sqr error");
    rn_check_eq(t,tt, NO, NO, "rf3_L_inf_norm error");

    if (verbose) { fprintf(stderr, "--- rf3_dist, rf3_dist_sqr, rf3_L_inf_dist ---\n"); }
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    r = rf3_dist(&a, &b);
    s = rf3_dist_sqr(&a, &b);
    t = rf3_L_inf_dist(&a, &b);
    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < N; i++)
      { double di = fabs(a.c[i] - b.c[i]);
        ss += di*di; 
        if (di > tt) { tt = di; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.0000001 * rr, NO, NO, "rf3_dist error");
    rn_check_eps(s,ss,0.0000001 * ss, NO, NO, "rf3_dist_sqr error");
    rn_check_eq(t,tt, NO, NO, "rf3_L_inf_dist error");

    if (verbose) { fprintf(stderr, "--- rf3_dir, rf3_L_inf_dir ---\n"); }
    a = rf3_throw_cube();
    b = rf3_dir(&a, &r);
    d = rf3_L_inf_dir(&a, &sf);
    ss = rf3_norm(&a);
    tf = rf3_L_inf_norm(&a);
    for (i = 0; i < N; i++)
      { float bi =  (float)(a.c[i]/ss);
        rn_check_eps(b.c[i], bi, 0.0000001 * ss, NO, NO, "rf3_dir error");
        float di = (float)((double)(a.c[i])/tf);
        rn_check_eps(d.c[i], di, 0.0000001 * tf, NO, NO, "rf3_L_inf_dir error");
      }

    if (verbose) { fprintf(stderr, "--- rf3_dot, rf3_cos, rf3_sin, rf3_angle ---\n"); }
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    r = rf3_dot(&a, &b);
    double S = rf3_sin(&a, &b);
    double C = rf3_cos(&a, &b);
    double A = rf3_angle(&a, &b);
    mag = sqrt(rf3_dot(&a,&a)*rf3_dot(&b,&b));
    rr = 0.0;
    for (i = 0; i < N; i++) { rr += a.c[i]*b.c[i]; }
    double CC = rr/(rf3_norm(&a)*rf3_norm(&b));
    rn_check_eps(r,rr,0.0000001 * mag, NO, NO, "rf3_dot error(1)");
    rn_check_eps(C,CC,0.0000001, NO, NO, "rf3_cos error(1)");
    d = a;
    rf3_mix_in(-rr/rf3_norm_sqr(&b), &b, &d);
    double SS = rf3_norm(&d)/rf3_norm(&a);
    rn_check_eps(S,SS,0.0000001, NO, NO, "rf3_sin error(1)");
    double AA = atan2(SS, CC);
    rn_check_eps(A,AA,0.0000001, NO, NO, "rf3_angle error(1)");
    for (i = 0; i < N; i++)
      { a = rf3_axis(i);
        for (j = 0; j < N; j++)
          { b = rf3_axis(j);
            r = rf3_dot(&a, &b);
            s = rf3_sin(&a, &b);
            t = rf3_cos(&a, &b);
            rr = (i == j ? 1.0 : 0.0);
            rn_check_eq(r,rr, NO, NO, "rf3_dot error(2)");
            rn_check_eq(t,rr, NO, NO, "rf3_dot error(2)");
            rn_check_eq(s,1.0 - rr, NO, NO, "rf3_dot error(2)");
          }
      }

    if (verbose) { fprintf(stderr, "--- rf3_cross ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        int32_t i2 = (i + 2) % N;
        int32_t p;
        a = rf3_axis(i0);
        b = rf3_axis(i1);
        d = rf3_cross(&a, &b);
        e = rf3_axis(i2);
        for (p = 0; p < N; p++)
          { float ep = e.c[p];
            rn_check_eq(d.c[p],ep, NO, NO, "rf3_cross error(x)");
          }
      }
    /* Test on random vectors: */
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    d = rf3_cross(&a, &b);
    mag = rf3_norm(&a)*rf3_norm(&b);
    r = rf3_dot(&a, &d);
    rn_check_eps(r,0.0,0.0000001 * mag*rf3_norm(&a), NO, NO, "rf3_cross error(1)");
    r = rf3_dot(&b, &d);
    rn_check_eps(r,0.0,0.0000001 * mag*rf3_norm(&b), NO, NO, "rf3_cross error(2)");

    if (verbose) { fprintf(stderr, "--- rf3_det ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        int32_t i2 = (i + 2) % N;
        a = rf3_axis(i0);
        b = rf3_axis(i1);
        c = rf3_axis(i2);
        r = rf3_det(&a, &b, &c);
        rn_check_eq(r,1.0, NO, NO, "rf3_det error(2)");
      }
    /* Test on random vectors: */
    a = rf3_throw_cube();
    b = rf3_throw_cube();
    c = rf3_throw_cube();
    r = rf3_det(&a, &b, &c);
    e = rf3_cross(&a, &b);
    rr = rf3_dot(&e, &c);
    mag = rf3_norm(&a)*rf3_norm(&b)*rf3_norm(&c);
    rn_check_eps(r,rr,0.0000001 * mag, NO, NO, "rf3_det error(1)");

    if (verbose) { fprintf(stderr, "--- rf3_decomp ---\n"); }

    a = rf3_throw_cube();
    b = rf3_throw_cube();
    r = rf3_decomp(&a, &b, &para, &perp);
    rr = rf3_dot(&a, &b)/rf3_norm_sqr(&b);  
    rn_check_eps(r,rr,0.0000001 * (fabs(r) + fabs(rr)), NO, NO, "rf3_decomp error(1)");
    c = rf3_add(&para, &perp);
    s = rf3_dist(&a, &c);
    affirm (s <= 0.0000001 * rf3_norm(&a), "rf3_decomp error(2)");
    s = rf3_dot(&perp, &b);
    rn_check_eps(s,0.0,0.0000001 * rf3_norm(&b), NO, NO, "rf3_decomp error(3)");
    t = rf3_dot(&para, &perp);
    rn_check_eps(t,0.0,0.0000001 * rf3_norm(&a), NO, NO, "rf3_decomp error(4)");

    if (verbose) { fprintf(stderr, "--- rf3_print ---\n"); }
    if (verbose)
      { a = rf3_throw_cube();
        fprintf(stderr, "a = ");
        rf3_print(stderr, &a);
        fputc('\n', stderr);
      }
  }
