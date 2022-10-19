/* rf2test --- test program for rf2.h, rf2x2.h  */
/* Last edited on 2022-01-04 08:43:13 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>
#include <jsrandom.h>
#include <flt.h>
#include <rf2.h>
#include <rn_test_tools.h>
#include <rfn_check.h>


#define N 2
#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_rf2(int32_t verbose);
void test_rf2maps(int32_t verbose);

int32_t main (int32_t argc, char **argv)
  {
    srand(1993);
    srandom(1993);

    for (int32_t i = 0; i < 100; i++) test_rf2(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_rf2(int32_t verbose)
  {
    rf2_t a, b, c, d, e, para, perp;
    double r, s, t;
    float sf, tf;
    double rr, ss, tt;
    double mag;
    int32_t i, j, k;

    if (verbose)
      { fprintf(stderr,
          "sizeof(rf2_t) = %lu  %d*sizeof(double) = %lu\n",
          sizeof(rf2_t), N, N*sizeof(double)
        );
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_zero ---\n"); }
    a = rf2_zero();
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i], 0.0, NO, NO, "rf2_zero error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_all ---\n"); }
    a = rf2_all(3.14f);
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i], 3.14f, NO, NO, "rf2_all error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_axis ---\n"); }
    for (k = 0; k < N; k++)
      { a = rf2_axis(k);
        for (i = 0; i < N; i++)
          { rn_check_eq(a.c[i],(i == k ? 1.0 : 0.0), NO, NO, "rf2_axis error"); }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_throw_cube ---\n"); }
    a = rf2_throw_cube();
    for (i = 0; i < N; i++)
      { affirm(a.c[i] != a.c[(i+1)%N], "rf2_throw probable error(1)"); 
        /* Check whether there are more than 8 nonzero bits: */
        double vv = a.c[i]*256.0;
        affirm(vv != floor(vv), "rf2_throw error(3)"); 
        affirm((a.c[i] > -1.0) && (a.c[i] < 1.0), "rf2_throw error(2)"); 
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_throw_dir ---\n"); }
    /* Should check uniformity... */
    a = rf2_throw_dir();
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "rf2_throw_dir error(1)"); }
    /* Check whether the norm is 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    rn_check_eps(1,rr,0.0000001 * rr, NO, NO, "rf2_throw_dir error (2)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_throw_ball ---\n"); }
    /* Should check uniformity... */
    a = rf2_throw_ball();
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "rf2_throw_ball error(1)"); }
    /* Check whether the norm is at most 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    demand(rr <= 1 + 0.0000001*rr, "rf2_throw_ball error (2)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_add ---\n"); }
    a = rf2_throw_cube();
    b = rf2_throw_cube();
    d = rf2_add(&a, &b);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] + b.c[i], NO, NO, "rf2_add error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_sub ---\n"); }
    a = rf2_throw_cube();
    b = rf2_throw_cube();
    d = rf2_sub(&a, &b);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] - b.c[i], NO, NO, "rf2_sub error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_neg ---\n"); }
    a = rf2_throw_cube();
    d = rf2_neg(&a);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],- a.c[i], NO, NO, "rf2_neg error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_scale ---\n"); }
    s = drandom();
    a = rf2_throw_cube();
    d = rf2_scale(s, &a);
    for (i = 0; i < N; i++)
      { float zi = (float)(s*a.c[i]);
        rn_check_eq(d.c[i],zi, NO, NO, "rf2_scale error(1)");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_mix ---\n"); }
    s = drandom();
    t = drandom();
    a = rf2_throw_cube();
    b = rf2_throw_cube();
    d = rf2_mix(s, &a, t, &b);
    for (i = 0; i < N; i++)
      { float ddi = (float)(s * a.c[i] + t * b.c[i]);
        rn_check_eq(d.c[i],ddi, NO, NO, "rf2_mix error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_mix_in ---\n"); }
    s = drandom();
    a = rf2_throw_cube();
    b = rf2_throw_cube();
    d = b;
    rf2_mix_in(s, &a, &d);
    for (i = 0; i < N; i++)
      { float ddi = (float)(b.c[i] + s * a.c[i]);
        rn_check_eq(d.c[i],ddi, NO, NO, "rf2_mix_in error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_weigh ---\n"); }
    a = rf2_throw_cube();
    b = rf2_throw_cube();
    d = rf2_weigh(&a, &b);
    for (i = 0; i < N; i++)
      { float ddi = (float)(((double)a.c[i]) * b.c[i]);
        rn_check_eq(d.c[i],ddi, NO, NO, "rf2_weigh error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_unweigh ---\n"); }
    a = rf2_throw_cube();
    b = rf2_throw_cube();
    d = rf2_unweigh(&a, &b);
    for (i = 0; i < N; i++)
      { float ddi = (float)(((double)a.c[i]) / b.c[i]);
        rn_check_eq(d.c[i],ddi, NO, NO, "rf2_unweigh error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_rot_axis ---\n"); }
    { a = rf2_throw_cube();
      double ang = 2.1*M_PI*drandom();
      d = rf2_rot(&a, ang);
      rfn_check_rot_axis(N, a.c, 0, 1, ang, d.c, "rf2_rot_axis error");
    }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_norm, rf2_norm_sqr, rf2_L_inf_norm ---\n"); }
    a = rf2_throw_cube();
    r = rf2_norm(&a);
    s = rf2_norm_sqr(&a);
    t = rf2_L_inf_norm(&a);
    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < N; i++)
      { float ai = fabsf(a.c[i]);
        ss += ((double)ai)*ai; 
        if (ai > tt) { tt = ai; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.0000001 * rr, NO, NO, "rf2_norm error");
    rn_check_eps(s,ss,0.0000001 * ss, NO, NO, "rf2_norm_sqr error");
    rn_check_eq(t,tt, NO, NO, "rf2_L_inf_norm error");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_dist, rf2_dist_sqr, rf2_L_inf_dist ---\n"); }
    a = rf2_throw_cube();
    b = rf2_throw_cube();
    r = rf2_dist(&a, &b);
    s = rf2_dist_sqr(&a, &b);
    t = rf2_L_inf_dist(&a, &b);
    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < N; i++)
      { double di = fabs(((double)a.c[i]) - b.c[i]);
        ss += di*di; 
        if (di > tt) { tt = di; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.0000001 * rr, NO, NO, "rf2_dist error");
    rn_check_eps(s,ss,0.0000001 * ss, NO, NO, "rf2_dist_sqr error");
    rn_check_eq(t,tt, NO, NO, "rf2_L_inf_dist error");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_dir, rf2_L_inf_dir ---\n"); }
    a = rf2_throw_cube();
    b = rf2_dir(&a, &r);
    d = rf2_L_inf_dir(&a, &sf);
    ss = rf2_norm(&a);
    tf = rf2_L_inf_norm(&a);
    for (i = 0; i < N; i++)
      { float bi =  (float)(a.c[i]/ss);
        rn_check_eps(b.c[i], bi, 0.0000001 * ss, NO, NO, "rf2_dir error");
        float di = (float)(((double)a.c[i])/tf);
        rn_check_eps(d.c[i], di, 0.0000001 * tf, NO, NO, "rf2_L_inf_dir error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_dot, rf2_cos, rf2_sin, rf2_angle ---\n"); }
    a = rf2_throw_cube();
    b = rf2_throw_cube();
    r = rf2_dot(&a, &b);
    double S = rf2_sin(&a, &b);
    double C = rf2_cos(&a, &b);
    double A = rf2_angle(&a, &b);
    mag = sqrt(rf2_dot(&a,&a)*rf2_dot(&b,&b));
    rr = 0.0;
    for (i = 0; i < N; i++) { rr += ((double)a.c[i])*b.c[i]; }
    double CC = rr/(rf2_norm(&a)*rf2_norm(&b));
    rn_check_eps(r,rr,0.0000001 * mag, NO, NO, "rf2_dot error(1)");
    rn_check_eps(C,CC,0.0000001, NO, NO, "rf2_cos error(1)");
    d = a;
    rf2_mix_in(-rr/rf2_norm_sqr(&b), &b, &d);
    double SS = rf2_norm(&d)/rf2_norm(&a);
    rn_check_eps(S,SS,0.0000001, NO, NO, "rf2_sin error(1)");
    double AA = atan2(SS, CC);
    rn_check_eps(A,AA,0.0000001, NO, NO, "rf2_angle error(1)");
    for (i = 0; i < N; i++)
      { a = rf2_axis(i);
        for (j = 0; j < N; j++)
          { b = rf2_axis(j);
            r = rf2_dot(&a, &b);
            s = rf2_sin(&a, &b);
            t = rf2_cos(&a, &b);
            rr = (i == j ? 1.0 : 0.0);
            rn_check_eq(r,rr, NO, NO, "rf2_dot error(2)");
            rn_check_eq(t,rr, NO, NO, "rf2_dot error(2)");
            rn_check_eq(s,1.0 - rr, NO, NO, "rf2_dot error(2)");
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_cross ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        double sign = ((i % 2) == 0 ? 1.0 : -1.0);
        int32_t p;
        a = rf2_axis(i0);
        d = rf2_cross(&a);
        e = rf2_axis(i1);
        for (p = 0; p < N; p++)
          { double ep = sign*e.c[p];
            rn_check_eq(d.c[p],ep, NO, NO, "rf2_cross error(x)");
          }
      }
    /* Test on random vectors: */
    a = rf2_throw_cube();
    d = rf2_cross(&a);
    mag = rf2_norm(&a);
    r = rf2_dot(&a, &d);
    rn_check_eps(r,0.0,0.00000001 * mag, NO, NO, "rf2_cross error(1)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_det ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        double sign = ((i % 2) == 0 ? 1.0 : -1.0);
        a = rf2_axis(i0);
        b = rf2_axis(i1);
        r = rf2_det(&a, &b);
        rn_check_eq(r,sign, NO, NO, "rf2_det error(2)");
      }
    /* Test on random vectors: */
    a = rf2_throw_cube();
    b = rf2_throw_cube();
    r = rf2_det(&a, &b);
    e = rf2_cross(&a);
    rr = rf2_dot(&e, &b);
    mag = rf2_norm(&a)*rf2_norm(&b);
    rn_check_eps(r,rr,0.00000001 * mag, NO, NO, "rf2_det error(1)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_decomp ---\n"); }
    a = rf2_throw_cube();
    b = rf2_throw_cube();
    r = rf2_decomp(&a, &b, &para, &perp);
    rr = rf2_dot(&a, &b)/rf2_norm_sqr(&b);  
    rn_check_eps(r,rr,0.0000001 * (fabs(r) + fabs(rr)), NO, NO, "rf2_decomp error(1)");
    c = rf2_add(&para, &perp);
    s = rf2_dist(&a, &c);
    affirm (s <= 0.0000001 * rf2_norm(&a), "rf2_decomp error(2)");
    s = rf2_dot(&perp, &b);
    rn_check_eps(s,0.0,0.0000001 * rf2_norm(&b), NO, NO, "rf2_decomp error(3)");
    t = rf2_dot(&para, &perp);
    rn_check_eps(t,0.0,0.0000001 * rf2_norm(&a), NO, NO, "rf2_decomp error(4)");
    
    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- rf2_print ---\n"); }
    if (verbose)
      { a = rf2_throw_cube ();
        fprintf(stderr, "a = ");
        rf2_print(stderr, &a);
        fputc('\n', stderr);
      }

    if (verbose)
      { 
        fprintf(stderr, "!! rf2_is_finite NOT TESTED\n");
        fprintf(stderr, "!! rf2_eq NOT TESTED\n");
        fprintf(stderr, "!! rf2_throw_normal NOT TESTED\n");
      }
  }

