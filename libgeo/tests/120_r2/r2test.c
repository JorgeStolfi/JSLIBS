/* r2test --- test program for r2.h, r2x2.h  */
/* Last edited on 2024-08-30 21:00:08 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>
#include <jsrandom.h>
#include <bool.h>
#include <flt.h>

#include <rn_test_tools.h>

#include <r2.h>
#include <r2x2.h>
#include <r2_extra.h>
#include <r2_bezier.h>

#define N 2
#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);

void test_r2(bool_t verbose);
void test_r2x2(bool_t verbose);
void test_r2_extra(bool_t verbose);
void test_r2_bezier(bool_t verbose);

void throw_matrix(r2x2_t *M);
void throw_diag_matrix(r2x2_t *M);
void throw_symmetric_matrix(r2x2_t *M);

void test_r2_zero(bool_t verbose);
void test_r2_all(bool_t verbose);
void test_r2_axis(bool_t verbose);
void test_r2_add(bool_t verbose);
void test_r2_sub(bool_t verbose);
void test_r2_neg(bool_t verbose);
void test_r2_scale(bool_t verbose);
void test_r2_mix(bool_t verbose);
void test_r2_mix_in(bool_t verbose);
void test_r2_weigh(bool_t verbose);
void test_r2_unweigh(bool_t verbose);
void test_r2_rot(bool_t verbose);
void test_r2_norm_sqr(bool_t verbose);
void test_r2_norm(bool_t verbose);
void test_r2_L_inf_norm(bool_t verbose);
void test_r2_dist_sqr(bool_t verbose);
void test_r2_dist(bool_t verbose);
void test_r2_L_inf_dist(bool_t verbose);
void test_r2_dir(bool_t verbose);
void test_r2_L_inf_dir(bool_t verbose);
void test_r2_dot(bool_t verbose);
void test_r2_cos(bool_t verbose);
void test_r2_sin(bool_t verbose);
void test_r2_angle(bool_t verbose);
void test_r2_cross(bool_t verbose);
void test_r2_det(bool_t verbose);
void test_r2_decomp(bool_t verbose);
void test_r2_throw_cube(bool_t verbose);
void test_r2_throw_dir(bool_t verbose);
void test_r2_throw_ball(bool_t verbose);
void test_r2_print(bool_t verbose);
void test_r2_gen_print(bool_t verbose);

void test_r2x2_size_allocation(bool_t verbose);
void test_r2x2_indexing_addressing(bool_t verbose);

void test_r2x2_zero(bool_t verbose);
void test_r2x2_ident(bool_t verbose);
void test_r2x2_transp(bool_t verbose);
void test_r2x2_get_row(bool_t verbose);
void test_r2x2_set_row(bool_t verbose);
void test_r2x2_get_col(bool_t verbose);
void test_r2x2_set_col(bool_t verbose);
void test_r2x2_map_row(bool_t verbose);
void test_r2x2_map_col(bool_t verbose);
void test_r2x2_scale(bool_t verbose);
void test_r2x2_mul(bool_t verbose);
void test_r2x2_mul_tr(bool_t verbose);
void test_r2x2_tr_mul(bool_t verbose);
void test_r2x2_det(bool_t verbose);
void test_r2x2_inv(bool_t verbose);
void test_r2x2_sqrt(bool_t verbose);
void test_r2x2_norm_sqr(bool_t verbose);
void test_r2x2_norm(bool_t verbose);
void test_r2x2_normalize(bool_t verbose);
void test_r2x2_mod_norm_sqr(bool_t verbose);
void test_r2x2_sym_eigen(bool_t verbose);
void test_r2x2_print(bool_t verbose);
void test_r2x2_gen_print(bool_t verbose);

void test_r2_map_projective(bool_t verbose);
void test_r2_map_twirl(bool_t verbose);
void test_r2_map_expand__r2_map_contract(bool_t verbose);

int32_t main (int32_t argc, char **argv)
  {
    srand(1993);
    srandom(1993);

    for (int32_t i = 0; i < 100; i++) test_r2(i < 3);
    for (int32_t i = 0; i < 100; i++) test_r2x2(i < 3);
    for (int32_t i = 0; i < 100; i++) test_r2_extra(i < 3);
    for (int32_t i = 0; i < 100; i++) test_r2_bezier(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

/* TESTS OF {r2.h} FUNCTIONS */

void test_r2(bool_t verbose)
  {
    if (verbose)
      { fprintf(stderr,
          "sizeof(r2_t) = %lu  %d*sizeof(double) = %lu\n",
          sizeof(r2_t), N, N*sizeof(double)
        );
      }

    test_r2_zero(verbose);
    test_r2_all(verbose);
    test_r2_axis(verbose);
    test_r2_add(verbose);
    test_r2_sub(verbose);
    test_r2_neg(verbose);
    test_r2_scale(verbose);
    test_r2_mix(verbose);
    test_r2_mix_in(verbose);
    test_r2_weigh(verbose);
    test_r2_unweigh(verbose);
    test_r2_rot(verbose);
    test_r2_norm_sqr(verbose);
    test_r2_norm(verbose);
    test_r2_L_inf_norm(verbose);
    test_r2_dist_sqr(verbose);
    test_r2_dist(verbose);
    test_r2_L_inf_dist(verbose);
    test_r2_dir(verbose);
    test_r2_L_inf_dir(verbose);
    test_r2_dot(verbose);
    test_r2_cos(verbose);
    test_r2_sin(verbose);
    test_r2_angle(verbose);
    test_r2_cross(verbose);
    test_r2_det(verbose);
    test_r2_decomp(verbose);
    test_r2_throw_cube(verbose);
    test_r2_throw_dir(verbose);
    test_r2_throw_ball(verbose);
    test_r2_gen_print(verbose);
    test_r2_print(verbose);

    if (verbose)
      { 
        fprintf(stderr, "!! r2_is_finite NOT TESTED\n");
        fprintf(stderr, "!! r2_eq NOT TESTED\n");
        fprintf(stderr, "!! r2_barycenter NOT TESTED\n");
        fprintf(stderr, "!! r2_bbox NOT TESTED\n");
        fprintf(stderr, "!! r2_orient NOT TESTED\n");
        fprintf(stderr, "!! r2_circumcenter NOT TESTED\n");
        fprintf(stderr, "!! r2_incircle NOT TESTED\n");
        fprintf(stderr, "!! r2_throw_normal NOT TESTED\n");
     }
  }

/* TESTS OF {r2.h} FUNCTIONS */

void test_r2_zero(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_zero ---\n"); }
    r2_t a;
    r2_zero(&a);
    for (int32_t i = 0; i < N; i++)
      { rn_check_eq(a.c[i],0.0, NO, NO, "r2_zero error"); }
  }

void test_r2_all(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_all ---\n"); }
    r2_t a;
    double val = M_PI/3;
    r2_all(val, &a);
    for (int32_t i = 0; i < N; i++)
      { rn_check_eq(a.c[i], val, NO, NO, "r2_all error"); }
  }

void test_r2_axis(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_axis ---\n"); }
    r2_t a;
    for (int32_t k = 0; k < N; k++)
      { r2_axis(k, &a);
        for (int32_t i = 0; i < N; i++)
          { double vi = (i == k ? 1.0 : 0.0);
            rn_check_eq(a.c[i], vi, NO, NO, "r2_axis error");
          }
      }
  }

void test_r2_add(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_add ---\n"); }
    r2_t a, b, d;
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r2_add(&a, &b, &d);
    for (int32_t i = 0; i < N; i++)
      { double vi = a.c[i] + b.c[i];
        rn_check_eq(d.c[i], vi, NO, NO, "r2_add error");
      }
  }
  
void test_r2_sub(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_sub ---\n"); }
    r2_t a, b, d;
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r2_sub(&a, &b, &d);
    for (int32_t i = 0; i < N; i++)
      { double vi = a.c[i] - b.c[i];
        rn_check_eq(d.c[i], vi, NO, NO, "r2_sub error");
      }
  }

void test_r2_neg(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_neg ---\n"); }
    r2_t a, d;
    r2_throw_cube(&a);
    r2_neg(&a, &d);
    for (int32_t i = 0; i < N; i++)
      { double vi = - a.c[i];
        rn_check_eq(d.c[i], vi, NO, NO, "r2_neg error");
      }
  }

void test_r2_scale(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_scale ---\n"); }
    r2_t a, d;
    double s = drandom();
    r2_throw_cube(&a);
    r2_scale(s, &a, &d);
    for (int32_t i = 0; i < N; i++)
      { double vi = s*a.c[i];
        rn_check_eq(d.c[i], vi, NO, NO, "r2_scale error(1)");
      }
  }

void test_r2_mix(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_mix ---\n"); }
    r2_t a, b, d;
    double s = drandom();
    double t = drandom();
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r2_mix(s, &a, t, &b, &d);
    for (int32_t i = 0; i < N; i++)
      { double vi = s * a.c[i] + t * b.c[i];
        rn_check_eq(d.c[i], vi, NO, NO, "r2_mix error");
      }
  }

void test_r2_mix_in(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_mix_in ---\n"); }
    r2_t a, b, d;
    double s = drandom();
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    d = b;
    r2_mix_in(s, &a, &d);
    for (int32_t i = 0; i < N; i++)
      { double vi = b.c[i] + s * a.c[i];
        rn_check_eq(d.c[i], vi, NO, NO, "r2_mix_in error");
      }
  }

void test_r2_weigh(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_weigh ---\n"); }
    r2_t a, b, d;
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r2_weigh(&a, &b, &d);
    for (int32_t i = 0; i < N; i++)
      { double vi = a.c[i] * b.c[i];
        rn_check_eq(d.c[i], vi, NO, NO, "r2_weigh error");
      }
  }

void test_r2_unweigh(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_unweigh ---\n"); }
    r2_t a, b, d;
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r2_unweigh(&a, &b, &d);
    for (int32_t i = 0; i < N; i++)
      { double vi = a.c[i] / b.c[i];
        rn_check_eq(d.c[i], vi, NO, NO, "r2_unweigh error");
      }
  }

void test_r2_rot(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_rot ---\n"); }
    r2_t a, d;
    { r2_throw_cube(&a);
      double ang = 2.1*M_PI*drandom();
      r2_rot(&a, ang, &d);
      rn_test_rot_axis(N, a.c, 0, 1, ang, d.c, "r2_rot error");
    }
  }

void test_r2_norm_sqr(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_norm_sqr ---\n"); }
    r2_t a;
    r2_throw_cube(&a);
    double s = r2_norm_sqr(&a);
    double ss = 0.0;
    for (int32_t i = 0; i < N; i++)
      { double ai = fabs(a.c[i]);
        ss += ai*ai; 
      }
    rn_check_eps(s,ss, 0.000000001*ss, NO, NO, "r2_norm_sqr error");
  }

void test_r2_norm(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_norm ---\n"); }
    r2_t a;
    r2_throw_cube(&a);
    double r = r2_norm(&a);
    double ss = 0.0;
    for (int32_t i = 0; i < N; i++)
      { double ai = fabs(a.c[i]);
        ss += ai*ai; 
      }
    double rr = sqrt(ss);
    rn_check_eps(r,rr, 0.000000001*rr, NO, NO, "r2_norm error");
  }

void test_r2_L_inf_norm(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_L_inf_norm ---\n"); }
    r2_t a;
    r2_throw_cube(&a);
    double t = r2_L_inf_norm(&a);
    double tt = 0.0;
    for (int32_t i = 0; i < N; i++)
      { double ai = fabs(a.c[i]);
        if (ai > tt) { tt = ai; }
      }
    rn_check_eq(t,tt, NO, NO, "r2_L_inf_norm error");
  }

void test_r2_dist_sqr(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_dist_sqr ---\n"); }
    r2_t a, b;
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    double s = r2_dist_sqr(&a, &b);
    double ss = 0.0;
    for (int32_t i = 0; i < N; i++)
      { double di = fabs(a.c[i] - b.c[i]);
        ss += di*di; 
      }
    rn_check_eps(s,ss, 0.000000001*ss, NO, NO, "r2_dist_sqr error");
  }

/* ---------------------------------------------------------------------- */
void test_r2_dist(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_dist ---\n"); }
    r2_t a, b;
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    double r = r2_dist(&a, &b);
    double ss = 0.0;
    for (int32_t i = 0; i < N; i++)
      { double di = fabs(a.c[i] - b.c[i]);
        ss += di*di; 
      }
    double rr = sqrt(ss);
    rn_check_eps(r,rr, 0.000000001*r, NO, NO, "r2_dist error");
  }

void test_r2_L_inf_dist(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_L_inf_dist ---\n"); }
    r2_t a, b;
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    double t = r2_L_inf_dist(&a, &b);
    double tt = 0.0;
    for (int32_t i = 0; i < N; i++)
      { double di = fabs(a.c[i] - b.c[i]);
        if (di > tt) { tt = di; }
      }
    rn_check_eq(t,tt, NO, NO, "r2_L_inf_dist error");
  }

void test_r2_dir(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_dir ---\n"); }
    r2_t a, d;
    r2_throw_cube(&a);
    double r = r2_dir(&a, &d);
    double rr = r2_norm(&a);
    rn_check_eps(r, rr, 0.000000001*rr, NO, NO, "r2_dir error (1)");
    for (int32_t i = 0; i < N; i++)
      { double vi = a.c[i]/rr;
        rn_check_eps(d.c[i], vi, 0.000000001*rr, NO, NO, "r2_dir error (2)");
      }
  }

void test_r2_L_inf_dir(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_L_inf_dir ---\n"); }
    r2_t a, d;
    r2_throw_cube(&a);
    double t = r2_L_inf_dir(&a, &d);
    double tt = r2_L_inf_norm(&a);
    rn_check_eq(t, tt, NO, NO, "r2_L_inf_dir error (1)");
    for (int32_t i = 0; i < N; i++)
      { double vi = a.c[i]/tt;
        rn_check_eps(d.c[i],vi, 0.000000001*tt, NO, NO, "r2_L_inf_dir error (2)");
      }
  }

void test_r2_dot(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_dot ---\n"); }
    r2_t a, b;
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    double r = r2_dot(&a, &b);
    double mag = sqrt(r2_dot(&a,&a)*r2_dot(&b,&b));
    double rr = 0.0;
    for (int32_t i = 0; i < N; i++) { rr += a.c[i]*b.c[i]; }
    rn_check_eps(r,rr, 0.000000001*mag, NO, NO, "r2_dot error(1)");
    for (int32_t i = 0; i < N; i++)
      { r2_axis(i, &a);
        for (int32_t j = 0; j < N; j++)
          { r2_axis(j, &b);
            double rij = r2_dot(&a, &b);
            double rrij = (i == j ? 1.0 : 0.0);
            rn_check_eq(rij,rrij, NO, NO, "r2_dot error(2)");
          }
      }
  }

void test_r2_cos(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_cos ---\n"); }
    r2_t a, b;
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    double C = r2_cos(&a, &b);
    double r = r2_dot(&a, &b);
    double CC = r/(r2_norm(&a)*r2_norm(&b));
    rn_check_eps(C,CC, 0.000000001, NO, NO, "r2_cos error(1)");
    for (int32_t i = 0; i < N; i++)
      { r2_axis(i, &a);
        for (int32_t j = 0; j < N; j++)
          { r2_axis(j, &b);
            double Cij = r2_cos(&a, &b);
            double CCij = (i == j ? 1.0 : 0.0);
            rn_check_eq(Cij,CCij,  NO, NO, "r2_cos error(2)");
          }
      }
  }

void test_r2_sin(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_sin ---\n"); }
    r2_t a, b, d;
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    double S = r2_sin(&a, &b);
    d = a;
    r2_mix_in(-r2_dot(&a, &b)/r2_norm_sqr(&b), &b, &d);
    double SS = r2_norm(&d)/r2_norm(&a);
    rn_check_eps(S,SS, 0.000000001, NO, NO, "r2_sin error(1)");
    for (int32_t i = 0; i < N; i++)
      { r2_axis(i, &a);
        for (int32_t j = 0; j < N; j++)
          { r2_axis(j, &b);
            double Sij = r2_sin(&a, &b);
            double SSij = (i == j ? 0.0 : 1.0);
            rn_check_eq(Sij,SSij, NO, NO, "r2_sin error(2)");
          }
      }
  }

void test_r2_angle(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_angle ---\n"); }
    r2_t a, b;
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    double S = r2_sin(&a, &b);
    double C = r2_cos(&a, &b);
    double A = r2_angle(&a, &b);
    double AA = atan2(S, C);
    rn_check_eps(A,AA, 0.000000001, NO, NO, "r2_angle error(1)");
    for (int32_t i = 0; i < N; i++)
      { r2_axis(i, &a);
        for (int32_t j = 0; j < N; j++)
          { r2_axis(j, &b);
            double Aij = r2_angle(&a, &b);
            double AAij = (i == j ? 0.0 : M_PI/2);
            rn_check_eps(Aij, AAij, 0.000000001, NO, NO, "r2_angle error(2)");
          }
      }
  }

void test_r2_cross(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_cross ---\n"); }
    r2_t a, d, e;
    /* Test on basis vectors: */
    for (int32_t i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        double sign = ((i % 2) == 0 ? 1.0 : -1.0);
        r2_axis(i0, &a);
        r2_cross(&a, &d);
        r2_axis(i1, &e);
        for (int32_t p = 0; p < N; p++)
          { double ep = sign*e.c[p];
            rn_check_eq(d.c[p],ep, NO, NO, "r2_cross error(x)");
          }
      }
    /* Test on random vectors: */
    r2_throw_cube(&a);
    r2_cross(&a, &d);
    double mag = r2_norm(&a);
    double r = r2_dot(&a, &d);
    rn_check_eps(r, 0.0, 0.00000001*mag, NO, NO, "r2_cross error(1)");
  }

void test_r2_det(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_det ---\n"); }
    r2_t a, b, e;
    /* Test on basis vectors: */
    for (int32_t i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        double sign = ((i % 2) == 0 ? 1.0 : -1.0);
        r2_axis(i0, &a);
        r2_axis(i1, &b);
        double r = r2_det(&a, &b);
        rn_check_eq(r,sign, NO, NO, "r2_det error(2)");
      }
    /* Test on random vectors: */
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    double r = r2_det(&a, &b);
    r2_cross(&a, &e);
    double rr = r2_dot(&e, &b);
    double mag = r2_norm(&a)*r2_norm(&b);
    rn_check_eps(r, rr, 0.00000001*mag, NO, NO, "r2_det error(1)");
  }

void test_r2_decomp(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_decomp ---\n"); }
    r2_t a, b, c, para, perp;
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    double r = r2_decomp(&a, &b, &para, &perp);
    double rr = r2_dot(&a, &b)/r2_norm_sqr(&b);  
    rn_check_eps(r, rr, 0.000000001*(fabs(r) + fabs(rr)), NO, NO, "r2_decomp error(1)");
    r2_add(&para, &perp, &c);
    double u = r2_dist(&a, &c);
    affirm (u <= 0.000000001*r2_norm(&a), "r2_decomp error(2)");
    double s = r2_dot(&perp, &b);
    rn_check_eps(s, 0.0, 0.000000001*r2_norm(&b), NO, NO, "r2_decomp error(3)");
    double t = r2_dot(&para, &perp);
    rn_check_eps(t, 0.0, 0.000000001*r2_norm(&a), NO, NO, "r2_decomp error(4)");
  }

void test_r2_throw_cube(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_throw_cube ---\n"); }
    r2_t a;
    r2_throw_cube(&a);
    for (int32_t i = 0; i < N; i++)
      { affirm(a.c[i] != a.c[(i+1)%N], "r2_throw probable error(1)"); 
        /* Check whether there are more than 8 nonzero bits: */
        double vv = a.c[i]*256.0;
        affirm(vv != floor(vv), "r2_throw error(3)"); 
        affirm((a.c[i] > -1.0) && (a.c[i] < 1.0), "r2_throw error(2)"); 
      }
  }

void test_r2_throw_dir(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_throw_dir ---\n"); }
    r2_t a;
    /* Should check uniformity... */
    r2_throw_dir(&a);
    /* Check variation: */
    for (int32_t i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "r2_throw_dir error(1)"); }
    /* Check whether the norm is 1: */
    double rr = 0;
    for (int32_t i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    rn_check_eps(1, rr, 0.000000001*rr, NO, NO, "r2_throw_dir error (2)");
  }

void test_r2_throw_ball(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_throw_ball ---\n"); }
    r2_t a;
    /* Should check uniformity... */
    r2_throw_ball(&a);
    /* Check variation: */
    for (int32_t i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "r2_throw_ball error(1)"); }
    /* Check whether the norm is at most 1: */
    double rr = 0;
    for (int32_t i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    demand(rr <= 1 + 0.000000001*rr, "r2_throw_ball error (2)");
  }

void test_r2_gen_print(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_gen_print ---\n"); }
    r2_t a;
    if (verbose)
      { r2_throw_cube (&a);
        fprintf(stderr, "a = ");
        r2_gen_print(stderr, &a, "%+8.3f", "« ", " : ", " »");
        fputc('\n', stderr);
      }
  }

void test_r2_print(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_print ---\n"); }
    r2_t a;
    if (verbose)
      { r2_throw_cube (&a);
        fprintf(stderr, "a = ");
        r2_print(stderr, &a);
        fputc('\n', stderr);
      }
  }

/* TESTS OF {r2x2.h} FUNCTIONS */

void test_r2x2(bool_t verbose)
  {
    test_r2x2_size_allocation(verbose);
    test_r2x2_indexing_addressing(verbose);

    test_r2x2_zero(verbose);
    test_r2x2_ident(verbose);
    test_r2x2_transp(verbose);
    test_r2x2_get_row(verbose);
    test_r2x2_get_col(verbose);
    test_r2x2_set_row(verbose);
    test_r2x2_set_col(verbose);
    test_r2x2_map_row(verbose);
    test_r2x2_map_col(verbose);
    test_r2x2_scale(verbose);
    test_r2x2_mul(verbose);
    test_r2x2_mul_tr(verbose);
    test_r2x2_tr_mul(verbose);
    test_r2x2_det(verbose);
    test_r2x2_inv(verbose);
    test_r2x2_sqrt(verbose);
    test_r2x2_norm_sqr(verbose);
    test_r2x2_norm(verbose);
    test_r2x2_normalize(verbose);
    test_r2x2_mod_norm_sqr(verbose);
    test_r2x2_gen_print(verbose);
    test_r2x2_print(verbose);
    test_r2x2_sym_eigen(verbose);

    if (verbose)
      { 
        fprintf(stderr, "!! r2x2_add NOT TESTED\n");
        fprintf(stderr, "!! r2x2_sub NOT TESTED\n");
        fprintf(stderr, "!! r2x2_neg NOT TESTED\n");
        fprintf(stderr, "!! r2x2_mix NOT TESTED\n");
        fprintf(stderr, "!! r2x2_adj NOT TESTED\n");
        fprintf(stderr, "!! r2x2_is_unif_scaling NOT TESTED\n");
        fprintf(stderr, "!! r2x2_rot90 NOT TESTED\n");
        fprintf(stderr, "!! r2x2_rot_and_scale NOT TESTED\n");
        fprintf(stderr, "!! r2x2_from_point_pairs NOT TESTED\n");
        fprintf(stderr, "!! r2x2_L_inf_norm NOT TESTED\n");
        fprintf(stderr, "!! r2x2_L_inf_normalize NOT TESTED\n");
        fprintf(stderr, "!! r2x2_diff_sqr NOT TESTED\n");
        fprintf(stderr, "!! r2x2_from_cols NOT TESTED\n");
        fprintf(stderr, "!! r2x2_from_rows NOT TESTED\n");
      }
  }  
  
void test_r2x2_size_allocation(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- Size and allocation ---\n"); }
    r2x2_t A, B, C;
    if (verbose)
      { fprintf(stderr,
          "sizeof(r2x2_t) = %lu  %d*%d*sizeof(double) = %lu\n",
          sizeof(r2x2_t), N, N, N*N*sizeof(double)
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
  }

void test_r2x2_indexing_addressing(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- Indexing and addressing ---\n"); }
    r2x2_t A;
    for (int32_t i = 0; i < N; i++)
      for (int32_t j = 0; j < N; j++)
        { double *Aij = &(A.c[i][j]); 
          affirm(Aij == ((double *)&A)+(N*i)+j, "r2x2_t indexing error");
        }
  }

void test_r2x2_zero(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_zero ---\n"); }
    r2x2_t A;
    r2x2_zero(&A);
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double vij = 0.0;
            rn_check_eq(A.c[i][j], vij, NO, NO, "r2x2_zero error"); 
          }
      }
  }

void test_r2x2_ident(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_ident ---\n"); }
    r2x2_t A;
    r2x2_ident(&A);
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double vij = (i == j ? 1.0 : 0.0);
            rn_check_eq(A.c[i][j], vij, NO, NO, "r2x2_ident error");
          }
      }
  }
  
void test_r2x2_transp(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_transp ---\n"); }
    r2x2_t A, B;
    throw_matrix(&A);
    r2x2_transp(&A, &B);
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double vij = A.c[j][i];
            rn_check_eq(B.c[i][j], vij, NO, NO, "r2x2_transp error (1)");
          }
      }
      
    /* In-place transpose: */
    B = A;
    r2x2_transp(&B, &B);
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double vij = A.c[j][i];
            rn_check_eq(B.c[i][j], vij, NO, NO, "r2x2_transp error (2)");
          }
      }
  }

void test_r2x2_get_row(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_get_row ---\n"); }
    r2x2_t A;
    r2_t a;
    throw_matrix(&A);
    for (int32_t i = 0; i < N; i++)
      { r2_throw_cube(&a);
        r2x2_get_row(&A, i, &a);
        for (int32_t j = 0; j < N; j++)
          { double vj = A.c[i][j];
            affirm(vj = a.c[j], "r2x2_get_row error");
          }
      }
  }

void test_r2x2_set_row(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_set_row ---\n"); }
    r2x2_t A;
    r2_t a;
    throw_matrix(&A);
    for (int32_t i = 0; i < N; i++)
      { r2_throw_cube(&a);
        r2x2_set_row(&A, i, &a);
        for (int32_t j = 0; j < N; j++)
          { double vj = A.c[i][j];
            affirm(vj = a.c[j], "r2x2_set_row error");
          }
      }
  }

void test_r2x2_get_col(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_get_col ---\n"); }
    r2x2_t A;
    r2_t a;
    throw_matrix(&A);
    for (int32_t i = 0; i < N; i++)
      { r2_throw_cube(&a);
        r2x2_get_col(&A, i, &a);
        for (int32_t j = 0; j < N; j++)
          { double vj = A.c[j][i];
            affirm(vj = a.c[j], "r2x2_get_col error");
          }
      }
  }

void test_r2x2_set_col(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_set_col ---\n"); }
    r2x2_t A;
    r2_t a;
    throw_matrix(&A);
    for (int32_t i = 0; i < N; i++)
      { r2_throw_cube(&a);
        r2x2_set_col(&A, i, &a);
        for (int32_t j = 0; j < N; j++)
          { double vj = A.c[j][i];
            affirm(vj = a.c[j], "r2x2_set_col error");
          }
      }
  }

void test_r2x2_map_row(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_map_row ---\n"); }
    r2x2_t A;
    r2_t a, b, bb;
    throw_matrix(&A);
    r2_throw_cube(&a);
    r2x2_map_row(&a, &A, &b);
    r2_zero(&bb);
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { bb.c[j] += a.c[i] * A.c[i][j]; }
      }
    double r = r2_dist(&b, &bb);
    affirm(r < 0.000000001*r2_norm(&bb), "r2_map_row error");
  }

void test_r2x2_map_col(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_map_col ---\n"); }
    r2x2_t A;
    r2_t a, b, bb;
    throw_matrix(&A);
    r2_throw_cube(&a);
    r2x2_map_col(&A, &a, &b);
    r2_zero(&bb);
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { bb.c[i] += A.c[i][j] * a.c[j];
          }
      }
    double s = r2_dist(&b, &bb);
    affirm(s < 0.000000001*r2_norm(&bb), "r2_map_col error");
  }

void test_r2x2_scale(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_scale ---\n"); }
    r2x2_t A, C;
    throw_matrix(&A);
    double r = drandom();
    r2x2_scale(r, &A, &C);
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double vij = r * A.c[i][j];
            rn_check_eps(C.c[i][j], vij, 0.000000001*fabs(vij), NO, NO,
              "r2x2_scale error"
            );
          }
      }
  }

void test_r2x2_mul(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_mul ---\n"); }
    r2x2_t A, B, C;
    throw_matrix(&A);
    throw_matrix(&B);
    r2x2_mul(&A, &B, &C);
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double sum = 0.0;
            for (int32_t k = 0; k < N; k++) { sum += A.c[i][k]*B.c[k][j]; }
            rn_check_eps(C.c[i][j], sum, 0.000000001*fabs(sum), NO, NO,
              "r2x2_mul error"
            );
          }
      }
  }

void test_r2x2_mul_tr(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_mul_tr ---\n"); }
    r2x2_t A, B, C;
    throw_matrix(&A);
    throw_matrix(&B);
    r2x2_mul_tr(&A, &B, &C);
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double sum = 0.0;
            for (int32_t k = 0; k < N; k++) { sum += A.c[i][k]*B.c[j][k]; }
            rn_check_eps(C.c[i][j], sum, 0.000000001*fabs(sum), NO, NO,
              "r2x2_mul_tr error"
            );
          }
      }
  }

void test_r2x2_tr_mul(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_tr_mul ---\n"); }
    r2x2_t A, B, C;
    throw_matrix(&A);
    throw_matrix(&B);
    r2x2_tr_mul(&A, &B, &C);
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double sum = 0.0;
            for (int32_t k = 0; k < N; k++) { sum += A.c[k][i]*B.c[k][j]; }
            rn_check_eps(C.c[i][j], sum, 0.000000001*fabs(sum), NO, NO,
              "r2x2_tr_mul error"
            );
          }
      }
  }

void test_r2x2_det(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_det ---\n"); }
    r2x2_t A;
    throw_matrix(&A);
    for (int32_t i = 0; i < N; i++)
      { int32_t k = (i + 1) % N;
        for (int32_t j = 0; j < N; j++)
          { /* Check for linearity */
            double r = drandom();
            A.c[i][j] = r;
            double rr = r2x2_det(&A);

            double s = drandom();
            A.c[i][j] = s;
            double ss = r2x2_det(&A);

            double t = drandom();
            A.c[i][j] = r*(1-t) + s*t;
            double tt = r2x2_det(&A);
            double mag = fabs(rr) + fabs(ss) + fabs(tt);
            rn_check_eps(tt, rr*(1.0-t)+ss*t, 0.000000001*mag, NO, NO,
              "r2x2_det error(1)"
            );
          }

        /* Row swap test: */
        double r = r2x2_det(&A);
        for (int32_t j = 0; j < N; j++)
          { double t = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = t; }
        double rr = r2x2_det(&A);
        double rmag = fabs(r) + fabs(rr);
        rn_check_eps(r, -rr, 0.000000001*rmag, NO, NO, "r2x2_det error(2)");

        /* Col swap test: */
        double s = r2x2_det(&A);
        for (int32_t j = 0; j < N; j++)
          { double t = A.c[j][i]; A.c[j][i] = A.c[j][k]; A.c[j][k] = t; }
        double ss = r2x2_det(&A);
        double smag = fabs(s) + fabs(ss);
        rn_check_eps(s, -ss, 0.000000001*smag, NO, NO, "r2x2_det error(3)");
      }
  }

void test_r2x2_inv(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_inv ---\n"); }
    r2x2_t A, B, C;
    throw_matrix(&A);
    r2x2_inv(&A, &B);
    r2x2_mul(&A, &B, &C);
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double vij = (i == j ? 1.0 : 0.0);
            affirm(fabs(C.c[i][j] - vij) < 0.000000001, "r2x2_inv error");
          }
      }
  }

void test_r2x2_sqrt(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_sqrt ---\n"); }
    r2x2_t A, B, C;
    do 
      { throw_matrix(&A); }
    while 
      ((r2x2_det(&A) < 0) || (A.c[0][0] + A.c[1][1] + 2*sqrt(r2x2_det(&A)) <= 1.0e-10));
    r2x2_sqrt(&A, &B);
    r2x2_mul(&B, &B, &C);
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { affirm(fabs(C.c[i][j] - A.c[i][j]) < 0.000000001, "r2x2_sqrt error"); }
      }
  }

void test_r2x2_norm_sqr(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_norm_sqr ---\n"); }
    r2x2_t A;
    throw_matrix(&A);
    double s = r2x2_norm_sqr(&A);
    double ss = 0;
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double Aij = A.c[i][j];
            ss += Aij*Aij;
          }
      }
    affirm(ss >= 0, "r2x2_norm_sqr error");
    affirm(fabs(ss - s) < 0.000000001, "r2x2_norm_sqr error");
  }

void test_r2x2_norm(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_norm ---\n"); }
    r2x2_t A;
    throw_matrix(&A);
    double r = r2x2_norm(&A);
    double ss = 0;
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double Aij = A.c[i][j];
            ss += Aij*Aij;
          }
      }
    double rr = sqrt(ss);
    affirm(fabs(rr - r) < 0.000000001, "r2x2_norm error");
  }

void test_r2x2_normalize(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_normalize ---\n"); }
    r2x2_t A, B;
    throw_matrix(&A);
    B = A;
    double s = r2x2_norm(&B);
    double ss = r2x2_normalize(&B);
    affirm(fabs(ss - s) < 0.000000001, "r2x2_normalize result error");
    double t = r2x2_norm(&B);
    double tt = 1.0;
    affirm(fabs(tt - t) < 0.000000001, "r2x2_normalize norm error");
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double Aij = A.c[i][j];
            double Bij = B.c[i][j];
            affirm(fabs(Bij*ss - Aij) < 0.000000001, "r2x2_normalize elem error");
          }
      }
  }

void test_r2x2_mod_norm_sqr(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_mod_norm_sqr ---\n"); }
    r2x2_t A;
    throw_matrix(&A);
    double t = r2x2_mod_norm_sqr(&A);
    double tt = 0;
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double Aij = A.c[i][j];
            double Dij = (i == j ? Aij - 1 : Aij);
            tt += Dij*Dij;
          }
      }
    affirm(tt >= 0, "r2x2_mod_norm_sqr error");
    affirm(fabs(tt - t) < 0.000000001, "r2x2_mod_norm_sqr error");
  }


void test_r2x2_sym_eigen(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_sym_eigen ---\n"); }
    r2x2_t A, B;
    r2_t a;
    if (drandom() < 0.2) 
      { /* Test with a diagonal matrix: */
        throw_diag_matrix(&A);
      }
    else
      { /* Test with a general symmetric matrix: */
        throw_symmetric_matrix(&A);
      }
    r2x2_sym_eigen(&A, &a, &B);
   /* Check order of eigenvalues: */
    for (int32_t i = 1; i < N; i++)
      { affirm(a.c[i-1] >= a.c[i], "r2x2_sym_eigen error: order"); }
    /* Check whether {B} is orthonormal: */
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double sum = 0.0;
            for (int32_t k = 0; k < N; k++) 
              { sum += B.c[k][i]*B.c[k][j]; }
            double val = (i == j ? 1.0 : 0.0);
            rn_check_eps(val, sum, 0.000000001 * fabs(sum), NO, NO, 
              "r2x2_sym_eigen error: not orthormal"
            );
          }
      }
    /* Check whether {B} is right-handed: */
    double rr = r2x2_det(&B);
    affirm(fabs(rr - 1.0) < 0.000000001, "r2x2_sym_eigen error: not right-handed");
    /* Check whether {A = B'*DIAG(e)*B}: */
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { double sum = 0.0;
            for (int32_t k = 0; k < N; k++) 
              { sum += B.c[k][i]*a.c[k]*B.c[k][j]; }
            rn_check_eps(A.c[i][j], sum, 0.000000001 * fabs(sum), NO, NO, 
              "r2x2_sym_eigen error: decomp"
            );
          }
      }
  }
void test_r2x2_gen_print(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_gen_print ---\n"); }
    r2x2_t A;
    if (verbose)
      { throw_matrix (&A);
        fprintf(stderr, "A = ");
        r2x2_gen_print(stderr, &A, "%+6.3f", " [|\n", "\n", "|]", "  « ", " | ", " »");
        fputc('\n', stderr);
      }
  }

void test_r2x2_print(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2x2_print ---\n"); }
    r2x2_t A;
    if (verbose)
      { throw_matrix (&A);
        fprintf(stderr, "A = ");
        r2x2_print(stderr, &A);
        fputc('\n', stderr);
      }
  }

/* TESTS OF {r2_extra.h} FUNCTIONS */

void test_r2_extra(bool_t verbose)
  {
    test_r2_map_projective(verbose);
    test_r2_map_twirl(verbose);
    test_r2_map_expand__r2_map_contract(verbose);
      
    if (verbose)
      { 
        fprintf(stderr, "!! r2_map_radial NOT TESTED\n");
        fprintf(stderr, "!! r2_map_compute_numeric_jacobian NOT TESTED\n");
        fprintf(stderr, "!! r2_map_check_jacobian NOT TESTED\n");
        fprintf(stderr, "!! r2_get_persp_rectangle_bbox NOT TESTED\n");
        fprintf(stderr, "!! r2_get_persp_disk_bbox NOT TESTED\n");
        fprintf(stderr, "!! r2_pixel_is_inside_persp_rectangle NOT TESTED\n");
        fprintf(stderr, "!! r2_pixel_is_inside_persp_disk NOT TESTED\n");
        fprintf(stderr, "!! r2_clip_seg_to_unit_disk NOT TESTED\n");
        fprintf(stderr, "!! r2_debug_point_jacobian NOT TESTED\n");
      }
  }

void test_r2_map_projective(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_map_projective ---\n"); }

    bool_t debug = FALSE;

    r3x3_t M = (r3x3_t){{{ 1.0, 2.0, 3.0 }, { 1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0 }}};
    
    auto void do_project(r2_t *p, r2x2_t *J);
    void do_project(r2_t *p, r2x2_t *J)
      { r2_map_projective(p, &M, p, J); }
    
    for (int32_t ii = 0; ii <= 2; ii++)
      { for (int32_t jj = 0; jj <= 2; jj++)
          { double X = ii/2.0;
            double Y = jj/2.0;
            r2_t a = (r2_t){{ X, Y }};
            r2x2_t J;    
            r2_t b = a;
            r2x2_ident(&J);
            r2_map_projective(&b, &M, &b, &J);
            r2_map_check_jacobian(&a, &do_project, "r2_map_projective", 1.0e-6, debug);
            r2_t c = (r2_t){{ (X + 2.0)/(X + 1.0), (Y + 3.0)/(X + 1.0) }};
            for (int32_t k = 0; k < N; k++)
              { rn_check_eps(b.c[k],c.c[k], 1.0e-8, NO, NO, "r2_map_projective error"); }
          }
      }
  }

void test_r2_map_twirl(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_map_twirl ---\n"); }
    
    bool_t debug = FALSE;

    int32_t NX = 640;
    int32_t NY = 480;
    r2_t ctr = (r2_t){{ 0.5*NX, 0.5*NY }};
    double rad = 0.25*fmin(NX,NY);
    double ang = 0.5*M_PI*(2*drandom()-1);
    
    auto void do_twirl(r2_t *p, r2x2_t *J);
    void do_twirl(r2_t *p, r2x2_t *J)
      { r2_map_twirl(p, &ctr, rad, ang, J); }
    
    for (int32_t ii = 0; ii <= 2; ii++)
      { for (int32_t jj = 0; jj <= 2; jj++)
          { double X = (0.10 + 0.80*ii/2.0)*NX;
            double Y = (0.10 + 0.80*jj/2.0)*NY;
            r2_t a = (r2_t){{ X, Y }};
            r2x2_t J;    
            r2_t b = a;
            r2x2_ident(&J);
            r2_map_twirl(&b, &ctr, rad, ang, &J);
            r2_map_check_jacobian(&a, &do_twirl, "r2_map_twirl", 1.0e-6, debug);
          }
      }
  }

void test_r2_map_expand__r2_map_contract(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_map_expand, r2_map_contract ---\n"); }
    
    bool_t debug = FALSE;
    
    double xlo = 2.0, xhi = 5.0;
    double ylo = 1.0, yhi = 3.0;
    
    auto void do_expand(r2_t *p, r2x2_t *J);
    void do_expand(r2_t *p, r2x2_t *J)
      { r2_map_expand(p, xlo, xhi, ylo, yhi, J); }
    
    auto void do_contract(r2_t *p, r2x2_t *J);
    void do_contract(r2_t *p, r2x2_t *J)
      { r2_map_contract(p, xlo, xhi, ylo, yhi, J); }
    
    for (int32_t ii = 0; ii <= 2; ii++)
      { for (int32_t jj = 0; jj <= 2; jj++)
          { double X = xlo + (0.05 + 0.90*jj/2)*(xhi - xlo);
            double Y = ylo + (0.05 + 0.90*ii/2)*(yhi - ylo);
            
            r2_t a = (r2_t){{ X, Y }};
            r2_map_check_jacobian(&a, &do_expand, "r2_map_expand", 1.0e-6, debug);
            
            r2x2_t JE;    
            r2_t b = (r2_t){{ X, Y }};
            
            r2_map_check_jacobian(&b, &do_contract, "r2_map_contract", 1.0e-6, debug);
            
            r2x2_ident(&JE);
            r2_map_expand(&b, xlo, xhi, ylo, yhi, &JE);
            r2_t c = b;
            r2x2_t JC = JE;    
            r2_map_contract(&c, xlo, xhi, ylo, yhi, &JC);
            double r = r2_dist(&a, &c);  
            rn_check_eps(r, 0.0, 0.0000001, NO, NO, "r2_map_expand/r2_map_contract error(1)");

            r2x2_t K;
            r2x2_ident(&K);
            double s = 0;
            for (int32_t i = 0; i < 2; i++)
              { for (int32_t j = 0; j < 2; j++)
                  { s += fabs(K.c[i][j] - JC.c[i][j]); }
              }
            rn_check_eps(s, 0.0, 0.000000001 * r2x2_norm(&JE), NO, NO, "r2_map_expand/r2_map_contract error(2)");
          }
      }
  }

/* TESTS OF {r2_bezier.h} FUNCTIONS */

void test_r2_bezier(bool_t verbose)
  {
    if (verbose)
      { 
        fprintf(stderr, "!! r2_bezier_length_estimate NOT TESTED\n");
        fprintf(stderr, "!! r2_bezier_eval NOT TESTED\n");
        fprintf(stderr, "!! r2_bezier_split NOT TESTED\n");
        fprintf(stderr, "!! r2_bezier_from_bend NOT TESTED\n");
      }
  }

/* TESTING TOOLS */

void throw_matrix(r2x2_t *M)
  {
    r2_t a;
    for (int32_t i = 0; i < N; i++)
      { r2_throw_cube(&a);
        for (int32_t j = 0; j < N; j++) { M->c[i][j] = a.c[j]; }
      }
  }

void throw_diag_matrix(r2x2_t *M)
  {
    for (int32_t i = 0; i < N; i++)
      { for (int32_t j = 0; j < N; j++)
          { M->c[i][j] = (i == j ? 2*drandom() - 1.0 : 0.0); }
      }
  }

void throw_symmetric_matrix(r2x2_t *M)
  {
    for (int32_t i = 0; i < N; i++)
      { /* Note: {j} runs to {i} not {N-1}! */
        for (int32_t j = 0; j <= i; j++)
          { M->c[i][j] = 2*drandom() - 1.0;
            if (j != i) { M->c[j][i] = M->c[i][j]; }
          }
      }
  }
