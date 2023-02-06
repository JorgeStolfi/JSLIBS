/* rntest --- test program for rn.h, rmxn.h  */
/* Last edited on 2023-02-03 05:37:09 by stolfi */

/* We need to set these in order to get {isnan}. What a crock... */
#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#include <affirm.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <flt.h>
#include <bool.h>

#include <rn.h>
#include <rmxn.h>
#include <rmxn_extra.h>
#include <rn_test_tools.h>

#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_rn(bool_t verbose);
void test_rmxn(bool_t verbose);
void test_rmxn_square(int32_t m, bool_t verbose);
void test_rmxn_rectangular(int32_t m, int32_t n, bool_t verbose);
void test_rmxn_rectangular_matched(int32_t m, int32_t p, int32_t n, bool_t verbose);

void test_rn_zero_all_axis(int32_t n, double a[], bool_t verbose);
void test_rn_throw(int32_t n, double a[], bool_t verbose);
void test_rn_add_sub_neg(int32_t n, double a[], double b[], double c[], bool_t verbose);
void test_rn_scale_shift_mix(int32_t n, double a[], double b[], double c[], bool_t verbose);
void test_rn_weigh_unweigh(int32_t n, double a[], double b[], double c[], bool_t verbose);
void test_rn_rot_axis(int32_t n, double a[], double c[], bool_t verbose);
void test_rn_sum(int32_t n, double a[], bool_t verbose);
void test_rn_norm_dist(int32_t n, double a[], double b[], bool_t verbose);
void test_rn_dir(int32_t n, double a[], double b[], double c[], bool_t verbose);
void test_rn_dot_trig(int32_t n, double a[], double b[], double c[], bool_t verbose);
void test_rn_cross(int32_t n, double b[], double c[], bool_t verbose);
void test_rn_det(int32_t n, double a[], bool_t verbose);
void test_rn_decomp(int32_t n, double a[], double b[], double c[], double d[], bool_t verbose);
void test_rn_print(int32_t n, double a[], bool_t verbose);
void test_rn_ball_vol(int32_t n, bool_t verbose);

void test_rmxn_zero_ident(int32_t m, int32_t n, double Amn[], double Bmn[], bool_t verbose);
void test_rmxn_copy(int32_t m, int32_t n, double Amn[], double Bmn[], bool_t verbose);
void test_rmxn_get_set_col(int32_t m, int32_t n, double Amn[], double Bmn[], double am[], double bm[], double cm[], bool_t verbose);
void test_rmxn_get_set_row(int32_t m, int32_t n, double Amn[], double Bmn[], double an[], double bn[], double cn[], bool_t verbose);
void test_rmxn_scale(int32_t m, int32_t n, double Amn[], double Bmn[], bool_t verbose);
void test_rmxn_mix(int32_t m, int32_t n, double Amn[], double Bmn[], double Cmn[], bool_t verbose);
void test_rmxn_rel_diff(int32_t m, int32_t n, double Amn[], double Bmn[], double Cmn[], bool_t verbose);
void test_rmxn_map_row(int32_t m, int32_t n, double Amn[], double am[], double bn[], double cn[], bool_t verbose);
void test_rmxn_map_col(int32_t m, int32_t n, double Amn[], double an[], double bm[], double cm[], bool_t verbose);
void test_rmxn_mul(int32_t m, int32_t p, int32_t n, double Amp[], double Bpn[], double Cmn[], bool_t verbose);
void test_rmxn_mul_tr(int32_t m, int32_t p, int32_t n, double Amn[], double Bpn[], double Amp[], bool_t verbose);
void test_rmxn_tr_mul(int32_t m, int32_t p, int32_t n, double Amp[], double Amn[], double Bpn[], bool_t verbose);
void test_rmxn_det(int32_t m, double Qmm[], bool_t verbose);
void test_rmxn_inv_cof_adj(int32_t m, double Qmm[], double Rmm[], double Smm[], bool_t verbose);
void test_rmxn_cholesky(int32_t m, double Qmm[], double Rmm[], double Smm[], bool_t verbose);
void test_rmxn_LT_inv_map(int32_t m, double Qmm[], double Rmm[], double am[], double bm[], double cm[], bool_t verbose);
void test_rmxn_LT_div(int32_t m, int32_t n, double Qmm[], double Rmm[], double Smm[], double Amn[], double Bmn[], double Cmn[], bool_t verbose);
void test_rmxn_mod_norm(int32_t m, double Qmm[], bool_t verbose);
void test_rmxn_norms(int32_t m, int32_t n, double Amn[], bool_t verbose);
void test_rmxn_print(int32_t m, int32_t n, bool_t verbose);
void test_rmxn_canonical_simplex_and_measures(int32_t n, bool_t verbose);
void test_rmxn_throw_canonical_simplex(int32_t n, bool_t verbose);
void test_rmxn_regular_simplex_and_measures(int32_t n, bool_t verbose);
void test_rmxn_throw_ortho_spin(int32_t n, bool_t verbose);
void test_rmxn_max_abs_elem(int32_t m, int32_t n, double Amn[], bool_t verbose);
void test_rmxn_perturb_unif(int32_t m, int32_t n, double Amn[], bool_t verbose);
void test_rmxn_transform_quadratic(int32_t m, int32_t n, double Amn[], bool_t verbose);


void throw_matrix(int32_t m, int32_t n, double *Amn);
void throw_LT_matrix(int32_t m, double *Lmm);
void print_matrix(FILE *wr, int32_t m, int32_t n, double *Amn);

void check_simplex(int32_t d, int32_t n, double V[], double rExp, double iExp, double sExp, double hExp, double mExp);
  /* Checks whether the {d}-dimensional simplex {V} of {R^n} is
    regular (modulo some roundoff). Assumes that {V} has {d+1} rows
    (vertices) and {n} columns (coordinates). Also checks the
    circum-radius {rExp}, the in-radius {iExp}, the edge length
    {sExp}, the height {hExp}, and the measure {mExp}. */

void check_ortho_matrix(int32_t n, double M[]);
  /* Checks whether the {n x n} matrix {M} orthonormal; that is
    whether the rows are pairwise orthogonal and have length 1. */

int32_t main (int32_t argc, char **argv)
  { int32_t i;
    srand(1993);
    srandom(1993);
    for (i = 0; i < 100; i++) test_rn(i <= 3);
    for (i = 0; i < 100; i++) test_rmxn(i <= 3);
    fclose(stderr);
    fclose(stdout);
    return 0;
  }


void test_rn (bool_t verbose)
  { double *a, *b, *c, *d;
    int32_t maxsize = (verbose ? 5 : 10);
    int32_t n = rand()/(RAND_MAX/maxsize) + 1;  /* Vector size. */

    fprintf(stderr, "test_rn:  n = %d\n", n);
    
    /* TEST: double *rn_alloc(int32_t n); */

    a = rn_alloc(n);
    b = rn_alloc(n);
    c = rn_alloc(n);
    d = rn_alloc(n);
    
    test_rn_zero_all_axis(n, a, verbose);
    test_rn_throw(n, a, verbose);
    test_rn_add_sub_neg(n, a, b, c, verbose);
    test_rn_scale_shift_mix(n, a, b, c, verbose);
    test_rn_weigh_unweigh(n, a, b, c, verbose);
    test_rn_rot_axis(n, a, c, verbose);
    test_rn_sum(n, a, verbose);
    test_rn_norm_dist(n, a, b, verbose);
    test_rn_dir(n, a, b, c, verbose);
    test_rn_dot_trig(n, a, b, c, verbose);
    test_rn_cross(n, b, c, verbose);
    test_rn_det(n, a, verbose);
    test_rn_decomp(n, a, b, c, d, verbose);
    test_rn_print(n, a, verbose);
    test_rn_ball_vol(n, verbose);
      
    /* NOT TESTED: */
    /* TEST: rn_copy (int32_t n, double *a, double *r); */
    /* TEST: double rn_mirror (int32_t n, double *a, double *u, double *r); */
    /* TEST: void rn_throw_normal (int32_t n, double *r); */
    /* TEST: double rn_abs_rel_diff(int32_t n, double *a, double *b, double abs_tol, double rel_tol); */
    /* TEST: void rn_do_check_eq(double x, double y, int32_t *i, int32_t *j, char *msg, rn_LOCPARMS); */
    /* TEST: void rn_do_check_eps(double x, double y, double eps, int32_t *i, int32_t *j, char *msg, rn_LOCPARMS); */

    free(a);
    free(b);
    free(c);
    free(d);
  }
  
void test_rn_zero_all_axis(int32_t n, double a[], bool_t verbose)
  {
    /* TEST: void rn_zero (int32_t n, double *r); */
    /* TEST: void rn_all (int32_t n, double x, double *r); */
    /* TEST: void rn_axis (int32_t n, int32_t i, double *r); */

    if (verbose) { fprintf(stderr, "--- rn_zero ---\n"); }
    rn_zero(n, a);
    for (int32_t i = 0; i < n; i++)
      { rn_check_eq(a[i],0.0, &i, NO, "rn_zero error"); }

    if (verbose) { fprintf(stderr, "--- rn_all ---\n"); }
    rn_all(n, 3.14, a);
    for (int32_t i = 0; i < n; i++)
      { rn_check_eq(a[i],3.14, &i, NO, "rn_all error"); }

    if (verbose) { fprintf(stderr, "--- rn_axis ---\n"); }
    for (int32_t k = 0; k < n; k++)
      { rn_axis(n, k, a);
        for (int32_t i = 0; i < n; i++)
          { rn_check_eq(a[i],(i == k ? 1.0 : 0.0), &i, NO, "rn_axis error"); }
      }
 }
    
void test_rn_throw(int32_t n, double a[], bool_t verbose)
  {

    /* TEST: void rn_throw_cube (int32_t n, double *r); */
    /* TEST: void rn_throw_dir (int32_t n, double *r); */
    /* TEST: void rn_throw_ball (int32_t n, double *r); */

    if (verbose) { fprintf(stderr, "--- rn_throw_cube ---\n"); }
    /* Should check uniformity... */
    rn_throw_cube(n, a);
    for (int32_t i = 0; i < n; i++)
      { affirm((n == 1) || (a[i] != a[(i+1)%n]), "rn_throw_cube probable error(1)"); 
        /* Check whether there are more than 8 nonzero bits: */
        double vv = a[i]*256.0;
        affirm(vv != floor(vv), "rn_throw_cube error(3)"); 
        affirm((a[i] > -1.0) && (a[i] < 1.0), "rn_throw_cube error(2)"); 
      }

    if (verbose) { fprintf(stderr, "--- rn_throw_dir ---\n"); }
    /* Should check uniformity... */
    rn_throw_dir(n, a);
    /* Check variation: */
    for (int32_t i = 0; i < n; i++) { affirm((n == 1) || (a[i] != a[(i+1)%n]), "rn_throw_dir error(1)"); }
    /* Check whether the norm is 1: */
    double rr = 0;
    for (int32_t i = 0; i < n; i++) { double ai = a[i]; rr += ai*ai; }
    rn_check_eps(1,rr,0.000000001 * rr, NO, NO, "rn_throw_dir error (2)");

    if (verbose) { fprintf(stderr, "--- rn_throw_ball ---\n"); }
    /* Should check uniformity... */
    rn_throw_ball(n, a);
    /* Check variation: */
    for (int32_t i = 0; i < n; i++) { affirm((n == 1) || (a[i] != a[(i+1)%n]), "rn_throw_ball error(1)"); }
    /* Check whether the norm is at most 1: */
    rr = 0;
    for (int32_t i = 0; i < n; i++) { double ai = a[i]; rr += ai*ai; }
    demand(rr <= 1 + 0.000000001*rr, "rn_throw_ball error (2)");
 }

void test_rn_add_sub_neg(int32_t n, double a[], double b[], double c[], bool_t verbose)
  {

    /* TEST: void rn_add (int32_t n, double *a, double *b, double *r); */
    /* TEST: void rn_sub (int32_t n, double *a, double *b, double *r); */
    /* TEST: void rn_neg (int32_t n, double *a, double *r); */

    if (verbose) { fprintf(stderr, "--- rn_add ---\n"); }
    rn_throw_cube(n, a);
    rn_throw_cube(n, b);
    rn_add(n, a, b, c);
    for (int32_t i = 0; i < n; i++)
      { rn_check_eq(c[i],a[i] + b[i], &i, NO, "rn_add error"); }

    if (verbose) { fprintf(stderr, "--- rn_sub ---\n"); }
    rn_throw_cube(n, a);
    rn_throw_cube(n, b);
    rn_sub(n, a, b, c);
    for (int32_t i = 0; i < n; i++)
      { rn_check_eq(c[i],a[i] - b[i], &i, NO, "rn_sub error"); }

    if (verbose) { fprintf(stderr, "--- rn_neg ---\n"); }
    rn_throw_cube(n, a);
    rn_neg(n, a, c);
    for (int32_t i = 0; i < n; i++)
      { double dxi = -a[i];
        rn_check_eq(c[i],dxi, &i, NO, "rn_neg error"); }
 }
    
void test_rn_scale_shift_mix(int32_t n, double a[], double b[], double c[], bool_t verbose)
  {
    /* TEST: void rn_scale (int32_t n, double s, double *a, double *r); */
    /* TEST: void rn_shift (int32_t n, double s, double *a, double *r); */
    /* TEST: void rn_mix (int32_t n, double s, double *a, double t, double *b, double *r); */
    /* TEST: void rn_mix_in (int32_t n, double s, double *a, double *r); */

    if (verbose) { fprintf(stderr, "--- rn_scale ---\n"); }
    double s = drandom();
    rn_throw_cube(n, a);
    rn_scale(n, s, a, c);
    for (int32_t i = 0; i < n; i++)
      { double dxi = s*a[i];
        fprintf(stderr, "%5d %24.16e %24.16e %24.16e %24.16e\n", i, s, a[i], c[i], dxi);
        rn_check_eq(c[i],dxi, &i, NO, "rn_scale error(1)");
      }

    if (verbose) { fprintf(stderr, "--- rn_shift ---\n"); }
    s = drandom();
    rn_throw_cube(n, a);
    rn_shift(n, s, a, c);
    for (int32_t i = 0; i < n; i++)
      { double dxi = s + a[i];
        fprintf(stderr, "%5d %24.16e %24.16e %24.16e %24.16e\n", i, s, a[i], c[i], dxi);
        rn_check_eq(c[i],dxi, &i, NO, "rn_scale error(1)");
      }

    if (verbose) { fprintf(stderr, "--- rn_mix ---\n"); }
    s = drandom();
    double t = drandom();
    rn_throw_cube(n, a);
    rn_throw_cube(n, b);
    rn_mix(n, s, a, t, b, c);
    for (int32_t i = 0; i < n; i++)
      { double dxi = s * a[i] + t * b[i];
        rn_check_eq(c[i],dxi, &i, NO, "rn_mix error");
      }

    if (verbose) { fprintf(stderr, "--- rn_mix_in ---\n"); }
    s = drandom();
    rn_throw_cube(n, a);
    rn_throw_cube(n, b);
    for (int32_t i = 0; i < n; i++) { c[i] = b[i]; }
    rn_mix_in(n, s, a, c);
    for (int32_t i = 0; i < n; i++)
      { double dxi = b[i] + s * a[i];
        rn_check_eq(c[i],dxi, &i, NO, "rn_mix_in error");
      }
 }
    
void test_rn_weigh_unweigh(int32_t n, double a[], double b[], double c[], bool_t verbose)
  {

    /* TEST: void rn_weigh (int32_t n, double *a, double *w, double *r); */
    /* TEST: void rn_unweigh (int32_t n, double *a, double *w, double *r); */

    if (verbose) { fprintf(stderr, "--- rn_weigh ---\n"); }
    rn_throw_cube(n, a);
    rn_throw_cube(n, b);
    rn_weigh(n, a, b, c);
    for (int32_t i = 0; i < n; i++)
      { double dxi = a[i] * b[i];
        rn_check_eq(c[i],dxi, &i, NO, "rn_weigh error");
      }

    if (verbose) { fprintf(stderr, "--- rn_unweigh ---\n"); }
    rn_throw_cube(n, a);
    rn_throw_cube(n, b);
    rn_unweigh(n, a, b, c);
    for (int32_t i = 0; i < n; i++)
      { double dxi = a[i] / b[i];
        rn_check_eq(c[i],dxi, &i, NO, "rn_unweigh error");
      }
 }
    
void test_rn_rot_axis(int32_t n, double a[], double c[], bool_t verbose)
  {
    /* TEST: void rn_rot_axis (int32_t n, double *a, int32_t i, int32_t j, double ang, double *r); */
    /* TEST: void rn_test_rot_axis(int32_t n, double *a, int32_t i, int32_t j, double ang, double *r, char *msg); */

    if (n >= 2)
      { if (verbose) { fprintf(stderr, "--- rn_rot_axis ---\n"); }
        { rn_throw_cube(n, a);
          int32_t i = int32_abrandom(0, n-1);
          int32_t j = int32_abrandom(0, n-2); if (j >= i) { j++; }
          double ang = 2.1*M_PI*drandom();
          rn_rot_axis(n, a, i, j, ang, c);
          rn_test_rot_axis(n, a, i, j, ang, c, "rn_rot_axis error");
        }
      }
 }
    
void test_rn_sum(int32_t n, double a[], bool_t verbose)
  {

    /* TEST: double rn_sum (int32_t n, double *a); */

    if (verbose) { fprintf(stderr, "--- rn_sum ---\n"); }
    rn_throw_cube(n, a);
    double s = rn_sum(n, a);
    double ss = 0.0;
    for (int32_t i = 0; i < n; i++) { ss += a[i]; }
    double rr = rn_L_inf_norm(n, a);
    rn_check_eps(s,ss, 0.000000001*rr, NO, NO, "rn_sum error");
 }
    
void test_rn_norm_dist(int32_t n, double a[], double b[], bool_t verbose)
  {

    /* TEST: double rn_norm (int32_t n, double *a); */
    /* TEST: double rn_norm_sqr (int32_t n, double *a); */
    /* TEST: double rn_L_inf_norm (int32_t n, double *a); */
    /* TEST: double rn_dist (int32_t n, double *a, double *b); */
    /* TEST: double rn_dist_sqr (int32_t n, double *a, double *b); */
    /* TEST: double rn_L_inf_dist (int32_t n, double *a, double *b); */

    if (verbose) { fprintf(stderr, "--- rn_norm, rn_norm_sqr, rn_L_inf_norm ---\n"); }
    rn_throw_cube(n, a);
    double r = rn_norm(n, a);
    double s = rn_norm_sqr(n, a);
    double t = rn_L_inf_norm(n, a);
    double ss = 0.0;
    double tt = 0.0;
    for (int32_t i = 0; i < n; i++)
      { double ai = fabs(a[i]);
        ss += ai*ai; 
        if (ai > tt) { tt = ai; }
      }
    double rr = sqrt(ss);
    rn_check_eps(r,rr,0.000000001 * rr, NO, NO, "rn_norm error");
    rn_check_eps(s,ss,0.000000001 * ss, NO, NO, "rn_norm_sqr error");
    rn_check_eq(t,tt, NO, NO, "rn_L_inf_norm error");

    if (verbose) { fprintf(stderr, "--- rn_dist, rn_dist_sqr, rn_L_inf_dist ---\n"); }
    rn_throw_cube(n, a);
    rn_throw_cube(n, b);
    r = rn_dist(n, a, b);
    s = rn_dist_sqr(n, a, b);
    t = rn_L_inf_dist(n, a, b);

    ss = 0.0;
    tt = 0.0;
    for (int32_t i = 0; i < n; i++)
      { double di = fabs(a[i] - b[i]);
        ss += di*di; 
        if (di > tt) { tt = di; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.000000001 * rr, NO, NO, "rn_dist error");
    rn_check_eps(s,ss,0.000000001 * ss, NO, NO, "rn_dist_sqr error");
    rn_check_eq(t,tt, NO, NO, "rn_L_inf_dist error");
 }
    
void test_rn_dir(int32_t n, double a[], double b[], double c[], bool_t verbose)
  {

    /* TEST: double rn_dir (int32_t n, double *a, double *r); */
    /* TEST: double rn_L_inf_dir (int32_t n, double *a, double *r); */

    if (verbose) { fprintf(stderr, "--- rn_dir, rn_L_inf_dir ---\n"); }
    rn_throw_cube(n, a);
    double r = rn_dir(n, a, b);
    assert(fabs(r) < sqrt(n) + 1.0e-12);
    double s = rn_L_inf_dir(n, a, c);
    assert(fabs(s) < 1.0 + 1.0e-12);
    double ss = rn_norm(n, a);
    double tt = rn_L_inf_norm(n, a);
    for (int32_t i = 0; i < n; i++)
      { rn_check_eps(b[i],a[i]/ss,0.000000001 * ss, NO, NO, "rn_dir error");
        rn_check_eps(c[i],a[i]/tt,0.000000001 * tt, NO, NO, "rn_L_inf_dir error");
      }
 }

void test_rn_dot_trig(int32_t n, double a[], double b[], double c[], bool_t verbose)
  {

    /* TEST: double rn_dot (int32_t n, double *a, double *b); */
    /* TEST: double rn_cos (int32_t n, double *a, double *b); */
    /* TEST: double rn_sin (int32_t n, double *a, double *b); */
    /* TEST: double rn_angle (int32_t n, double *a, double *b); */

    if (verbose) { fprintf(stderr, "--- rn_dot, rn_cos, rn_sin, rn_angle ---\n"); }
    rn_throw_cube(n, a);
    rn_throw_cube(n, b);
    double r = rn_dot(n, a, b);
    double S = rn_sin(n, a, b);
    double C = rn_cos(n, a, b);
    double A = rn_angle(n, a, b);

    double mag = sqrt(rn_dot(n, a,a)*rn_dot(n, b,b));
    double rr = 0.0;
    for (int32_t i = 0; i < n; i++) { rr += a[i]*b[i]; }
    double CC = rr/(rn_norm(n, a)*rn_norm(n, b));
    rn_check_eps(r,rr,0.000000001 * mag, NO, NO, "rn_dot error(1)");
    rn_check_eps(C,CC,0.000000001, NO, NO, "rn_cos error(1)");
    for (int32_t i = 0; i < n; i++) { c[i] = a[i]; }
    rn_mix_in(n, -rr/rn_norm_sqr(n, b), b, c);
    double SS = rn_norm(n, c)/rn_norm(n, a);
    rn_check_eps(S,SS,0.000000001, NO, NO, "rn_sin error(1)");
    double AA = atan2(SS, CC);
    rn_check_eps(A,AA,0.000000001, NO, NO, "rn_angle error(1)");
    for (int32_t i = 0; i < n; i++)
      { rn_axis(n, i, a);
        for (int32_t j = 0; j < n; j++)
          { rn_axis(n, j, b);
            double r = rn_dot(n, a, b);
            double s = rn_sin(n, a, b);
            double t = rn_cos(n, a, b);
            double rr = (i == j ? 1.0 : 0.0);
            rn_check_eq(r,rr, &i, &j, "rn_dot error(2)");
            rn_check_eq(t,rr, &i, &j, "rn_dot error(3)");
            rn_check_eq(s,1.0 - rr, &i, &j, "rn_dot error(4)");
          }
      }
 }

void test_rn_cross(int32_t n, double b[], double c[], bool_t verbose)
  {

    /* TEST: void rn_cross (int32_t n, double **a, double *r); */

    if (verbose) { fprintf(stderr, "--- rn_cross ---\n"); }
    { double **z = (double **)malloc((n-1)*sizeof(double *));
      double *magz = rn_alloc(n);
      /* Test on basis vectors: */
      for (int32_t i = 0; i < n; i++)
        { double sign = ((n-1)*i % 2 == 0 ? +1.0 : -1.0);
          for (int32_t k = 0; k < n-1; k++)
            { int32_t ik = (i + k) % n; 
              z[k] = rn_alloc(n);
              rn_axis(n, ik, z[k]);
            }
          { int32_t in1 = (i + n-1) % n; rn_axis(n, in1, b); }
          rn_cross(n, z, c);
          for (int32_t j = 0; j < n; j++)
            { double cxj = sign*b[j];
              rn_check_eq(c[j],cxj, &i, &j, "rn_cross error(x)");
            }
        }
      /* Test on random vectors: */
      double mag = 1.0;
      for (int32_t k = 0; k < n-1; k++)
        { z[k] = rn_alloc(n);
          rn_throw_cube(n, z[k]);
          { magz[k] = rn_norm(n, z[k]); mag *= magz[k]; }
        }
      rn_cross(n, z, c);
      for (int32_t k = 0; k < n-1; k++)
        { double r = rn_dot(n, z[k], c);
          rn_check_eps(r,0.0,0.00000001 * mag*magz[k], NO, NO, "rn_cross error(1)");
        }
      for (int32_t k = 0; k < n-1; k++) { free(z[k]); }
      free(z); free(magz);
    }
 }
    
void test_rn_det(int32_t n, double c[], bool_t verbose)
  {
    
    /* TEST: double rn_det (int32_t n, double **a); */

    if (verbose) { fprintf(stderr, "--- rn_det ---\n"); }
    { double **z = (double **)malloc(n*sizeof(double *));
      double *magz = rn_alloc(n);
      /* Test on basis vectors: */
      for (int32_t i = 0; i < n; i++)
        { double sign = ((n-1)*i % 2 == 0 ? +1.0 : -1.0);
          for (int32_t k = 0; k < n; k++)
            { int32_t ik = (i + k) % n; 
              z[k] = rn_alloc(n);
              rn_axis(n, ik, z[k]);
            }
          double r = rn_det(n, z);
          rn_check_eq(r,sign, &i, NO, "rn_det error(2)");
        }
      
      /* Test on random vectors (consistency with {rn_cross}): */
      double mag = 1.0;
      for (int32_t k = 0; k < n; k++)
        { z[k] = rn_alloc(n);
          rn_throw_cube(n, z[k]);
          { magz[k] = rn_norm(n, z[k]); mag *= magz[k]; }
        }
      double r = rn_det(n, z);
      rn_cross(n, z, c);
      double rr = rn_dot(n, c, z[n-1]);
      rn_check_eps(r,rr,0.00000001 * mag, NO, NO, "rn_det error(1)");

      for (int32_t k = 0; k < n; k++) { free(z[k]); }
      free(z); free(magz);
    }
 }

void test_rn_decomp(int32_t n, double a[], double b[], double c[], double d[], bool_t verbose)
  {
    double *para, *perp;
    para = rn_alloc(n);
    perp = rn_alloc(n);

    /* TEST: double rn_decomp (int32_t n, double *a, double *u, double *para, double *perp); */

    if (verbose) { fprintf(stderr, "--- rn_decomp ---\n"); }
    rn_throw_cube(n, a);
    rn_throw_cube(n, b);
    double r = rn_decomp(n, a, b, para, perp);
    double rr = rn_dot(n, a, b)/rn_norm_sqr(n, b);  
    rn_check_eps(r,rr,0.000000001 * (fabs(r) + fabs(rr)), NO, NO, "rn_decomp error(1)");
    rn_add(n, para, perp, c);
    double s = rn_dist(n, a, c);
    rn_check_eps(s,0.0,0.000000001 * rn_norm(n, a), NO, NO, "rn_decomp error(2)");
    s = rn_dot(n, perp, b);
    rn_check_eps(s,0.0,0.000000001 * rn_norm(n, b), NO, NO, "rn_decomp error(3)");
    double t = rn_dot(n, para, perp);
    rn_check_eps(t,0.0,0.000000001 * rn_norm(n, a), NO, NO, "rn_decomp error(4)");

    if (verbose) { fprintf(stderr, "--- rn_mirror ---\n"); }
    rn_throw_cube(n, a);
    rn_throw_dir(n, c);
    r = rn_mirror(n, a, c, b);
    /* The dot products must be equal and opposite: */
    double tol = 1.0e-6*rn_norm(n,a);
    r = rn_dot(n, a, c);  
    s = rn_dot(n, b, c);
    rn_check_eps(r,-s,tol, NO, NO, "rn_mirror error(1)");
    /* Compute the average {d} of {a} and {b}: */
    rn_mix(n, 0.5, a, 0.5, b, d);
    /* The average of {a} and {b} must be orthogonal to {c}: */
    t = rn_dot(n, c, d);
    rn_check_eps(t,0,tol, NO, NO, "rn_mirror error(2)");

    free(para);
    free(perp);
 }

void test_rn_print(int32_t n, double a[], bool_t verbose)
  {

    /* TEST: void rn_print (FILE *f, int32_t n, double *a); */
    /* TEST: void rn_gen_print */

    if (verbose) { fprintf(stderr, "--- rn_print, rn_gen_print ---\n"); }
    if (verbose)
      { rn_throw_cube (n, a);
        fprintf(stderr, "a = ");
        rn_print(stderr, n, a);
        fputc('\n', stderr);
        rn_gen_print(stderr, n, a, "%6.3f", "point = < ", " | ", " >");
        fputc('\n', stderr);
      }
  }
     
void test_rn_ball_vol(int32_t n, bool_t verbose)
  {
    /* TEST: double rn_ball_vol(double r, int32_t d); */
    /* TEST: double rn_ball_cap_vol_frac_pos(int32_t d, double u); */
    /* TEST: double rn_ball_cap_vol_frac_ang(int32_t d, double w); */

    fprintf(stderr, "test_rn_ball_vol:  n = %d\n", n);
    
    if (verbose) { fprintf(stderr, "--- testing rn_ball_vol ---\n"); }
    { 
      for (int32_t ir = 1; ir <= 3; ir++)
        { double r = (double)ir;
          double v = rn_ball_vol(r, n);
          if (verbose) { fprintf(stderr, "  measure of %d-ball with radius %.3f = %8.5f\n", n, r, v); } 
          /* Checking: */
          if (n <= 4)
            { double vv;
              if (n == 1)
                { vv = 2*r; }
              else if (n == 2)
                { vv = M_PI*r*r; }
              else if (n == 3)
                { vv = M_PI*4/3*r*r*r; }
              else if (n == 4)
                { vv = M_PI*M_PI/2*r*r*r*r; }
              else { assert(FALSE); vv = 0.0; }
              affirm(fabs(vv - v) <= 1.0e-6, "rn_ball_vol error");
            }
        }
    }
    
    if (verbose) { fprintf(stderr, "--- testing rn_ball_cap_vol_frac_ang ---\n"); }
    { int32_t NW = 10; /* Number of latitude steps. */
      int32_t imin = 0;
      int32_t imax = NW;
      double wmax = M_PI;
      if (verbose) { fprintf(stderr, "  volume fraction between lat = 0 and lat = z:\n"); }
      for (int32_t i = imin; i <= imax; i++)
        { double w = wmax*((double)i)/((double)NW);
          double f = rn_ball_cap_vol_frac_ang(n, w);
          if (verbose) { fprintf(stderr, "    %8.5f  %8.5f\n", w, f); } 
          /* Checking: */
          if (n <= 4)
            { double ff;
              if (n == 1)
                { ff = sin(w)/2; }
              else if (n == 2)
                { ff = (w + sin(w)*cos(w))/M_PI; }
              else if (n == 3)
                { ff = sin(w)*(3 - sin(w)*sin(w))/4; }
              else if (n == 4)
                { ff = (w + 2*sin(2*w)/3 + sin(4*w)/12)/M_PI; }
              else { assert(FALSE); ff = 0.0; }
              affirm(fabs(ff - f) <= 1.0e-6, "rn_ball_cap_vol_frac_ang error");
            }
        }
      fprintf(stderr, "\n");
    }

    if (verbose) { fprintf(stderr, "--- testing rn_ball_cap_vol_frac_pos ---\n"); }
    { int32_t NX = 10; /* Number of position steps in each hemisphere. */
      int32_t imin = -NX-1;
      int32_t imax = +NX+1;
      double xmax = 1.0;
      if (verbose) { fprintf(stderr, "  volume fraction between x = -1 and x = z:\n"); }
      for (int32_t i = imin; i <= imax; i++)
        { double u = xmax*((double)i)/((double)NX);
          double f = rn_ball_cap_vol_frac_pos(n, u);
          if (verbose) { fprintf(stderr, "  %8.5f  %8.5f\n", u, f); }  
          /* !!! Should check the value of {f} for {n <= 4} !!! */
        }
      fprintf(stderr, "\n");
    }
      
  }

void test_rmxn(bool_t verbose)
  {
    int32_t maxsize = (verbose ? 5 : 10);
    int32_t m = rand()/(RAND_MAX/maxsize) + 1;  /* Number of rows. */
    int32_t p = rand()/(RAND_MAX/maxsize) + 1;  /* Middle dimension for rmxn_mul. */
    int32_t n = rand()/(RAND_MAX/maxsize) + 1;  /* Number of columns. */

    test_rmxn_square(m, verbose);
    test_rmxn_rectangular(m, n, verbose);
    test_rmxn_rectangular_matched(m, p, n, verbose);
  }

void test_rmxn_square(int32_t m, bool_t verbose)
  { 
    /* Tests functions of {rmxn.h} that use only square matrices. */

    fprintf(stderr, "test_rmxn_square:  m = %d\n", m);

    double *Qmm = rmxn_alloc(m, m);
    double *Rmm = rmxn_alloc(m, m);
    double *Smm = rmxn_alloc(m, m);

    double *am = rn_alloc(m);
    double *bm = rn_alloc(m);
    double *cm = rn_alloc(m);

    test_rmxn_mod_norm(m, Qmm, verbose);
    test_rmxn_det(m, Qmm, verbose);
    test_rmxn_inv_cof_adj(m, Qmm, Rmm, Smm, verbose);
    test_rmxn_cholesky(m, Qmm, Rmm, Smm, verbose);
    test_rmxn_LT_inv_map(m, Qmm, Rmm, am, bm, cm, verbose);
    
    free(Qmm); free(Rmm); free(Smm);
    free(am); free(bm); free(cm);
  } 

void test_rmxn_rectangular(int32_t m, int32_t n, bool_t verbose)
  { 
    /* Tests functions of {rmxn.h} that deal with rectangular {m} by {n} matrices. */

    /* TEST: double *rmxn_alloc(int32_t m, int32_t n); */

    fprintf(stderr, "test_rmxn_rectangular:  m = %d  n = %d\n", m, n);
    
    double *Amn = rmxn_alloc(m, n);
    double *Bmn = rmxn_alloc(m, n);
    double *Cmn = rmxn_alloc(m, n);

    double *Qmm = rmxn_alloc(m, m);
    double *Rmm = rmxn_alloc(m, m);
    double *Smm = rmxn_alloc(m, m);

    double *am = rn_alloc(m);
    double *bm = rn_alloc(m);
    double *cm = rn_alloc(m);
    
    double *an = rn_alloc(n);
    double *bn = rn_alloc(n);
    double *cn = rn_alloc(n);

    test_rmxn_zero_ident(m, n, Amn, Bmn, verbose);
    test_rmxn_copy(m, n, Amn, Bmn, verbose);

    test_rmxn_get_set_col(m, n, Amn, Bmn, am, bm, cm, verbose);
    test_rmxn_get_set_row(m, n, Amn, Bmn, an, bn, cn, verbose);

    test_rmxn_scale(m, n, Amn, Bmn, verbose);
    test_rmxn_mix(m, n, Amn, Bmn, Cmn, verbose);
    
    test_rmxn_rel_diff(m, n, Amn, Bmn, Cmn, verbose);
    test_rmxn_map_row(m, n, Amn, am, bn, cn, verbose);
    test_rmxn_map_col(m, n, Amn, an, bm, cm, verbose);

    test_rmxn_LT_div(m, n, Qmm, Rmm, Smm, Amn, Bmn, Cmn, verbose);
    test_rmxn_norms(m, n, Amn, verbose);
    test_rmxn_canonical_simplex_and_measures(n, verbose);
    test_rmxn_throw_canonical_simplex(n, verbose);
    test_rmxn_regular_simplex_and_measures(n, verbose);
    test_rmxn_throw_ortho_spin(n, verbose);
    
    test_rmxn_max_abs_elem(m, n, Amn, verbose);
    test_rmxn_perturb_unif(m, n, Amn, verbose);
    test_rmxn_transform_quadratic(m, n, Amn, verbose);
    
    test_rmxn_print(m, n, verbose);

    free(Amn); free(Bmn); free(Cmn);
    free(Qmm); free(Rmm); free(Smm);
    free(am); free(bm); free(cm);
    free(an); free(bn); free(cn);
  } 

void test_rmxn_rectangular_matched(int32_t m, int32_t p, int32_t n, bool_t verbose)
  { 
    /* Tests functions of {rmxn.h} that deal with rectangular matrices with matched dimensions. */

    fprintf(stderr, "test_rmxn_rectangular_matched:  m = %d  p = %d  n = %d\n", m, p, n);
    
    double *Amn = rmxn_alloc(m, n);
    double *Bmn = rmxn_alloc(m, n);
    double *Cmn = rmxn_alloc(m, n);

    double *Amp = rmxn_alloc(m, p);
    double *Bpn = rmxn_alloc(p, n);

    test_rmxn_mul(m, p, n, Amp, Bpn, Cmn, verbose);
    test_rmxn_mul_tr(m, p, n, Amn, Bpn, Amp, verbose);
    test_rmxn_tr_mul(m, p, n, Amp, Amn, Bpn, verbose);
    
    free(Amn); free(Bmn); free(Cmn);
    free(Amp); free(Bpn);
  }  
     
void test_rmxn_zero_ident(int32_t m, int32_t n, double Amn[], double Bmn[], bool_t verbose)
  {
    /* TEST: void rmxn_zero(int32_t m, int32_t n, double *M); */
    /* TEST: void rmxn_ident(int32_t m, int32_t n, double *M); */

    if (verbose) { fprintf(stderr, "--- rmxn_zero, rmxn_ident ---\n"); }
    rmxn_zero(m, n, Amn);
    rmxn_ident(m, n, Bmn);
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { rn_check_eq(Amn[n*i + j],0.0, &i, &j, "rmxn_zero error"); 
            rn_check_eq(Bmn[n*i + j],(i == j ? 1.0 : 0.0), &i, &j, "rmxn_ident error");
          }
      }
  }
       
void test_rmxn_copy(int32_t m, int32_t n, double Amn[], double Bmn[], bool_t verbose)
  {
    /* TEST: void rmxn_copy(int32_t m, int32_t n, double *A, double *M); */

    if (verbose) { fprintf(stderr, "--- rmxn_copy ---\n"); }
    throw_matrix(m, n, Amn);
    rmxn_copy(m, n, Amn, Bmn);
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { rn_check_eq(Amn[n*i + j],Bmn[n*i + j], &i, &j, "rmxn_copy error"); }
      }
  }
    
void test_rmxn_get_set_col(int32_t m, int32_t n, double Amn[], double Bmn[], double am[], double bm[], double cm[], bool_t verbose)
  {
    /* TEST: void rmxn_get_col(int32_t m, int32_t n, double *A, int32_t j, double *r); */
    /* TEST: void rmxn_set_col(int32_t m, int32_t n, double *A, int32_t j, double *r); */

    if (verbose) { fprintf(stderr, "--- rmxn_get_col, rmxn_set_col ---\n"); }
    throw_matrix(m, n, Amn);
    rn_throw_cube(m, am);
    rmxn_copy(m, n, Amn, Bmn);
    rn_copy(m, am, bm);
    int32_t k = int32_abrandom(0, n-1);
    rmxn_get_col(m, n, Bmn, k, cm);
    rmxn_set_col(m, n, Bmn, k, bm);
    for (int32_t i = 0; i < m; i++)
      { /* Check whether {rmxn_get_col} copied the column correctly: */
        double ciObs = cm[i];
        double ciExp = Amn[i*n + k];
        rn_check_eq(ciObs,ciExp, &i, &k, "rmxn_get_col error (1)"); 
        /* Check whether {rmxn_set_col} modified the vector arg: */
        double biObs = bm[i];
        double biExp = am[i];
        rn_check_eq(biObs,biExp, &i, &k, "rmxn_set_col error (1)");
        for (int32_t j = 0; j < n; j++)
          { double BijObs = Bmn[i*n + j];
            double BijExp = (j == k ? am[i] : Amn[i*n + j]);
            rn_check_eq(BijObs,BijExp, &i, &j, "rmxn_set_col error (2)");
          }
      }
  }

void test_rmxn_get_set_row(int32_t m, int32_t n, double Amn[], double Bmn[], double an[], double bn[], double cn[], bool_t verbose)
  {
    /* TEST: void rmxn_get_row(int32_t m, int32_t n, double *A, int32_t i, double *r); */
    /* TEST: void rmxn_set_row(int32_t m, int32_t n, double *A, int32_t i, double *r); */
 
    if (verbose) { fprintf(stderr, "--- rmxn_get_row, rmxn_set_row ---\n"); }
    throw_matrix(m, n, Amn);
    rn_throw_cube(n, an);
    rmxn_copy(m, n, Amn, Bmn);
    rn_copy(n, an, bn);
    int32_t k = int32_abrandom(0, m-1);
    rmxn_get_row(m, n, Bmn, k, cn);
    rmxn_set_row(m, n, Bmn, k, bn);
    for (int32_t j = 0; j < n; j++)
      { /* Check whether {rmxn_get_row} copied the row correctly: */
        double cjObs = cn[j];
        double cjExp = Amn[k*n + j];
        rn_check_eq(cjObs,cjExp, &k, &j, "rmxn_get_row error (1)"); 
        /* Check whether {rmxn_set_row} modified the vector arg: */
        double bjObs = bn[j];
        double bjExp = an[j];
        rn_check_eq(bjObs,bjExp, &k, &j, "rmxn_set_row error (1)");
        for (int32_t i = 0; i < m; i++)
          { double BijObs = Bmn[i*n + j];
            double BijExp = (i == k ? an[j] : Amn[i*n + j]);
            rn_check_eq(BijObs,BijExp, &i, &j, "rmxn_set_row error (2)");
          }
      }
  }
    
void test_rmxn_scale(int32_t m, int32_t n, double Amn[], double Bmn[], bool_t verbose)
  {
    /* TEST: void rmxn_scale(int32_t m, int32_t n, double s, double *A, double *M); */

    if (verbose) { fprintf(stderr, "--- rmxn_scale ---\n"); }
    double s = drandom();
    throw_matrix(m, n, Amn);
    rmxn_scale(m, n, s, Amn, Bmn);
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { double zij = s*Amn[n*i + j];
            rn_check_eq(Bmn[n*i + j],zij, &i, &j, "rmxn_scale error(1)");
          }
      }
  }
    
void test_rmxn_mix(int32_t m, int32_t n, double Amn[], double Bmn[], double Cmn[], bool_t verbose)
  {
    /* TEST: void rmxn_mix (int32_t m, int32_t n, double s, double *A, double t, double *B, double *M); */

    if (verbose) { fprintf(stderr, "--- rmxn_mix ---\n"); }
    double s = drandom();
    double t = drandom();
    throw_matrix(m, n, Amn);
    throw_matrix(m, n, Bmn);
    rmxn_mix(m, n, s, Amn, t, Bmn, Cmn);
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { double zij = s*Amn[n*i + j] + t*Bmn[n*i + j];
            rn_check_eq(Cmn[n*i + j],zij, &i, &j, "rmxn_mix error(1)");
          }
      }
  }
    
void test_rmxn_rel_diff(int32_t m, int32_t n, double Amn[], double Bmn[], double Cmn[], bool_t verbose)
  {
    /* TEST: void rmxn_rel_diff(int32_t m, int32_t n, double *A, double *B, double *M); */

    if (verbose) { fprintf(stderr, "--- rmxn_rel_diff ---\n"); }
    throw_matrix(m, n, Amn);
    throw_matrix(m, n, Bmn);
    rmxn_rel_diff(m, n, Amn, Bmn, Cmn);
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { double dij = rel_diff(Amn[n*i + j], Bmn[n*i + j]);
            rn_check_eq(Cmn[n*i + j],dij, &i, &j, "rmxn_rel_diff error(1)");
          }
      }
  }
    
void test_rmxn_map_row(int32_t m, int32_t n, double Amn[], double am[], double bn[], double cn[], bool_t verbose)
  {
    /* TEST: void rmxn_map_row (int32_t m, int32_t n, double *x, double *A, double *r); */

    if (verbose) { fprintf(stderr, "--- rmxn_map_row ---\n"); }
    throw_matrix(m, n, Amn);
    rn_throw_cube(m, am);
    rmxn_map_row(m, n, am, Amn, bn);
    rn_zero(n, cn);
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { cn[j] += am[i] * Amn[n*i + j]; }
      }
    double s = rn_dist(n, bn, cn);
    rn_check_eps(s,0.0,0.000000001 * rn_norm(n, cn), NO, NO, "rn_map_row error");
  }

void test_rmxn_map_col(int32_t m, int32_t n, double Amn[], double an[], double bm[], double cm[], bool_t verbose)
  {
    /* TEST: void rmxn_map_col (int32_t m, int32_t n, double *A, double *x, double *r); */

    if (verbose) { fprintf(stderr, "--- rmxn_map_col ---\n"); }
    throw_matrix(m, n, Amn);
    rn_throw_cube(n, an);
    rmxn_map_col(m, n, Amn, an, bm);
    rn_zero(m, cm);
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { cm[i] += Amn[n*i + j] * an[j]; }
      }
    double r = rn_dist(m, bm, cm);
    rn_check_eps(r,0.0,0.000000001 * rn_norm(m, cm), NO, NO, "rn_map_col error");
  }
    
void test_rmxn_mul(int32_t m, int32_t p, int32_t n, double Amp[], double Bpn[], double Cmn[], bool_t verbose)
  {
    /* TEST: void rmxn_mul (int32_t m, int32_t p, int32_t n, double *A, double *B, double *M); */

    if (verbose) { fprintf(stderr, "--- rmxn_mul ---\n"); }
    throw_matrix(m, p, Amp);
    throw_matrix(p, n, Bpn);
    rmxn_mul(m, p, n, Amp, Bpn, Cmn);
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { double sum = 0.0;
            for (int32_t k = 0; k < p; k++) { sum += Amp[p*i + k]*Bpn[n*k + j]; }
            rn_check_eps(Cmn[n*i + j],sum,0.000000001 * fabs(sum), NO, NO, 
              "rmxn_mul error"
            );
          }
      }
  }
  
void test_rmxn_mul_tr(int32_t m, int32_t p, int32_t n, double Amn[], double Bpn[], double Amp[], bool_t verbose)
  {
    /* TEST: void rmxn_mul_tr (int32_t m, int32_t n, int32_t p, double *A, double *B, double *M); */

    if (verbose) { fprintf(stderr, "--- rmxn_mul_tr ---\n"); }
    throw_matrix(m, n, Amn);
    throw_matrix(p, n, Bpn);
    rmxn_mul_tr(m, p, n, Amn, Bpn, Amp);
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < p; j++)
          { double sum = 0.0;
            for (int32_t k = 0; k < n; k++) { sum += Amn[n*i + k]*Bpn[n*j + k]; }
            rn_check_eps(Amp[p*i + j],sum,0.000000001 * fabs(sum), &i, &j, 
              "rmxn_mul_tr error"
            );
          }
      }
  }
  
void test_rmxn_tr_mul(int32_t m, int32_t p, int32_t n, double Amp[], double Amn[], double Bpn[], bool_t verbose)
  {
    /* TEST: void rmxn_tr_mul (int32_t p, int32_t m, int32_t n, double *A, double *B, double *M); */

    if (verbose) { fprintf(stderr, "--- rmxn_tr_mul ---\n"); }
    throw_matrix(m, p, Amp);
    throw_matrix(m, n, Amn);
    rmxn_tr_mul(m, p, n, Amp, Amn, Bpn);
    for (int32_t i = 0; i < p; i++)
      { for (int32_t j = 0; j < n; j++)
          { double sum = 0.0;
            for (int32_t k = 0; k < m; k++) { sum += Amp[k*p + i]*Amn[k*n + j]; }
            rn_check_eps(Bpn[i*n + j],sum,0.000000001 * fabs(sum), &i, &j, 
              "rmxn_tr_mul error"
            );
          }
      }
  }
    
void test_rmxn_det(int32_t m, double Qmm[], bool_t verbose)
  {
    /* TEST: double rmxn_det (int32_t n, double *A); */

    if (verbose) { fprintf(stderr, "--- rmxn_det ---\n"); }
    throw_matrix(m, m, Qmm);
    for (int32_t i = 0; i < m; i++)
      { int32_t k = (i + 1) % m;
        for (int32_t j = 0; j < m; j++)
          { /* Check for linearity */
            double r = drandom();
            Qmm[m*i + j] = r;
            double rr = rmxn_det(m, Qmm);

            double s = drandom();
            Qmm[m*i + j] = s;
            double ss = rmxn_det(m, Qmm);

            double t = drandom();
            Qmm[m*i + j] = r*(1-t) + s*t;
            double tt = rmxn_det(m, Qmm);

            double mag = fabs(rr) + fabs(ss) + fabs(tt);
            rn_check_eps(tt,(rr*(1-t) + ss*t),000000001 * mag, NO, NO, 
              "rmxn_det error(1)"
            );
          }

        /* Row swap test: */
        double r = rmxn_det(m, Qmm);
        for (int32_t j = 0; j < m; j++)
          { double *Aij = &(Qmm[m*i + j]);
            double *Akj = &(Qmm[m*k + j]);
            { double t = *Aij; *Aij = *Akj; *Akj = t; }
          }
        double rr = rmxn_det(m, Qmm);
        double mag = fabs(r) + fabs(rr);
        rn_check_eps(r,(-rr),000000001 * mag, NO, NO, "rmxn_det error(2)");

        /* Col swap test: */
        r = rmxn_det(m, Qmm);
        for (int32_t j = 0; j < m; j++)
          { double *Aji = &(Qmm[m*j + i]);
            double *Ajk = &(Qmm[m*j + k]);
            { double t = *Aji; *Aji = *Ajk; *Ajk = t; }
          }
        rr = rmxn_det(m, Qmm);
        mag = fabs(r) + fabs(rr);
        rn_check_eps(r,(-rr),000000001 * mag, NO, NO, "rmxn_det error(3)");
      }
  }
    
void test_rmxn_inv_cof_adj(int32_t m, double Qmm[], double Rmm[], double Smm[], bool_t verbose)
  {
    /* TEST: double rmxn_cof (int32_t n, double *A, int32_t ix, int32_t jx); */
    /* TEST: void rmxn_adj (int32_t n, double *A, double *M); */
    /* TEST: double rmxn_inv (int32_t n, double *A, double *M); */
    /* TEST: double rmxn_inv_full (int32_t n, double *A, double *M); */

    if (verbose) { fprintf(stderr, "--- rmxn_cof ---\n"); }
    if (verbose) { fprintf(stderr, "!! warning: rmxn_cof not tested\n"); }
    
    if (verbose) { fprintf(stderr, "--- rmxn_adj ---\n"); }
    if (verbose) { fprintf(stderr, "!! warning: rmxn_adj not tested\n"); }
    
    if (verbose) { fprintf(stderr, "--- rmxn_inv ---\n"); }
    throw_matrix(m, m, Qmm);
    double r = rmxn_det(m, Qmm);
    double s = rmxn_inv(m, Qmm, Rmm);
    rmxn_mul(m, m, m, Qmm, Rmm, Smm);
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < m; j++)
          { double val = (i == j ? 1.0 : 0.0);
            rn_check_eps(Smm[m*i + j],val,0.000000001, NO, NO, "rmxn_inv error");
          }
      }
    rn_check_eps(r,s,0.000000001, NO, NO, "rmxn_inv/rmxn_det error");
    
    if (verbose) { fprintf(stderr, "--- rmxn_inv_full ---\n"); }
    s = rmxn_inv_full(m, Qmm, Rmm);
    rmxn_mul(m, m, m, Qmm, Rmm, Smm);
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < m; j++)
          { double val = (i == j ? 1.0 : 0.0);
            rn_check_eps(Smm[m*i + j],val,0.000000001, NO, NO, "rmxn_inv_full error");
          }
      }
    rn_check_eps(r,s,0.000000001, NO, NO, "rmxn_inv-full/rmxn_det error");
  }
    
void test_rmxn_cholesky(int32_t m, double Qmm[], double Rmm[], double Smm[], bool_t verbose)
  {
    /* TEST: void rmxn_cholesky(int32_t n, double *A, double *L); */

    if (verbose) { fprintf(stderr, "--- rmxn_cholesky ---\n"); }
    { double *Lmm = Qmm, *Jmm = Smm;
      throw_LT_matrix(m, Jmm);
      for (int32_t i = 0; i < m; i++) { Jmm[m*i + i] = fabs(Jmm[m*i + i]) + 1.0; }
      rmxn_mul_tr(m, m, m, Jmm, Jmm, Rmm);
      rmxn_cholesky(m, Rmm, Lmm);
      for (int32_t i = 0; i < m; i++)
        { for (int32_t j = 0; j < m; j++)
            { double eps = 0.000000001;
              rn_check_eps(Lmm[m*i + j],Jmm[m*i + j],eps, NO, NO, 
                "rmxn_cholesky error"
              );
            }
        }
    }
  }
    
void test_rmxn_LT_inv_map(int32_t m, double Qmm[], double Rmm[], double am[], double bm[], double cm[], bool_t verbose)
  {
    /* TEST: void rmxn_LT_inv_map_row(int32_t n, double *y, double *L, double *r); */
    /* TEST: void rmxn_LT_inv_map_col(int32_t m, double *L, double *y, double *r); */

    if (verbose) { fprintf(stderr, "--- rmxn_LT_inv_map_row, rmxn_LT_inv_map_col ---\n"); }
    { double *Lmm = Qmm;
      int32_t by_row; 
      throw_LT_matrix(m, Lmm);
      (void)rmxn_inv(m, Lmm, Rmm);
      rn_throw_cube(m, am);
      for (by_row = 0; by_row < 2; by_row++)
        { char *msg;
          if (by_row) 
            { rmxn_LT_inv_map_row(m, am, Lmm, bm);
              rmxn_map_row(m, m, am, Rmm, cm);
              msg = "rn_LT_inv_map_row error";
            }
          else
            { rmxn_LT_inv_map_col(m, Lmm, am, bm);
              rmxn_map_col(m, m, Rmm, am, cm);
              msg = "rn_LT_inv_map_col error";
            }
          double r = rn_dist(m, bm, cm);
          rn_check_eps(r,0.0,0.000000001 * rn_norm(m, cm), NO, NO, msg);
        }
    }
  }
    
void test_rmxn_LT_div(int32_t m, int32_t n, double Qmm[], double Rmm[], double Smm[], double Amn[], double Bmn[], double Cmn[], bool_t verbose)
  {
    /* TEST: void rmxn_LT_pos_div(int32_t m, int32_t n, double *A, double *L, double *M); */
    /* TEST: void rmxn_LT_pre_div(int32_t m, int32_t n, double *L, double *A, double *M); */

    if (verbose) { fprintf(stderr, "--- rmxn_LT_pos_div ---\n"); }
    { double *Lmm = Qmm, *Jmm = Rmm;
      /* We must use "Amn" as if it were "n x m" not "m x n". */
      double *Anm = Amn, *Bnm = Bmn, *Cnm = Cmn; 
      throw_LT_matrix(m, Lmm);
      (void)rmxn_inv(m, Lmm, Jmm);
      throw_matrix(n, m, Anm);
      rmxn_LT_pos_div(n, m, Anm, Lmm, Cnm);
      rmxn_mul(n, m, m, Anm, Jmm, Bnm);
      for (int32_t i = 0; i < n; i++)
        { for (int32_t j = 0; j < m; j++)
            { double eps = 0.000000001*(fabs(Bnm[m*i + j]) + 1.0e-200);
              rn_check_eps(Cnm[m*i + j],Bnm[m*i + j],eps, NO, NO, 
                "rmxn_LT_pos_div error"
              );
            }
        }
    }
    
    if (verbose) { fprintf(stderr, "--- rmxn_LT_pre_div ---\n"); }
    { double *Lmm = Qmm, *Jmm = Rmm;
      throw_LT_matrix(m, Lmm);
      (void)rmxn_inv(m, Lmm, Jmm);
      throw_matrix(m, n, Amn);
      rmxn_LT_pre_div(m, n, Lmm, Amn, Cmn);
      rmxn_mul(m, m, n, Jmm, Amn, Bmn);
      for (int32_t i = 0; i < m; i++)
        { for (int32_t j = 0; j < n; j++)
            { double eps = 0.000000001*fabs(Bmn[n*i + j]);
              rn_check_eps(Cmn[n*i + j],Bmn[n*i + j],eps, NO, NO, 
                "rmxn_LT_pre_div error"
              );
            }
        }
    }
  }
    
    
void test_rmxn_mod_norm(int32_t m, double Qmm[], bool_t verbose)
  {
    /* TEST: double rmxn_mod_norm_sqr (int32_t n, double *A); */

    if (verbose) { fprintf(stderr, "--- rmxn_mod_norm ---\n"); }
    throw_matrix(m, m, Qmm);
    double t = rmxn_mod_norm_sqr(m,Qmm);
    double tt = 0;
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < m; j++)
          { double Aij = Qmm[m*i + j]; 
            double Dij = (i == j ? Aij - 1 : Aij);
            tt += Dij*Dij;
          }
      }
    affirm(tt >= 0, "rmxn_mod_norm_sqr error");
    affirm(fabs(tt - t) < 000000001, "rmxn_mod_norm_sqr error");
  }

void test_rmxn_norms(int32_t m, int32_t n, double Amn[], bool_t verbose)
  {
    /* TEST: double rmxn_norm_sqr(int32_t m, int32_t n, double *A); */
    /* TEST: double rmxn_norm(int32_t m, int32_t n, double *A); */
    /* TEST: double rmxn_mod_norm_sqr (int32_t n, double *A); */

    if (verbose) { fprintf(stderr, "--- rmxn_norm,rmxn_norm_sqr ---\n"); }
    throw_matrix(m,n,Amn);
    double s = rmxn_norm_sqr(m,n,Amn);
    double r = rmxn_norm(m,n,Amn);
    double ss = 0;
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { double Aij = Amn[n*i + j]; 
            ss += Aij*Aij;
          }
      }
    affirm(ss >= 0, "rmxn_norm_sqr error");
    affirm(fabs(ss - s) < 000000001, "rmxn_norm_sqr error");
    double rr = sqrt(ss);
    affirm(fabs(rr - r) < 000000001, "rmxn_norm error");
  }
    
void test_rmxn_print(int32_t m, int32_t n, bool_t verbose)
  {
    /* TEST: void rmxn_print (FILE *f, int32_t m, int32_t n, double *A); */
    /* TEST: void rmxn_gen_print  */

    if (verbose) {  }
    if (verbose)
      { fprintf(stderr, "--- rmxn_print, rmxn_gen_print, rmxn_gen_print2, rmxn_gen_print3, ---\n");

        int32_t n1 = n % 5;
        int32_t n2 = (n + 2) % 5;
        int32_t n3 = (n + 4) % 5;

        double *Amn1 = rmxn_alloc(m, n1);
        double *Amn2 = rmxn_alloc(m, n2);
        double *Amn3 = rmxn_alloc(m, n3);

        throw_matrix (m, n1, Amn1);
        throw_matrix (m, n2, Amn2);
        throw_matrix (m, n3, Amn3);
        
        fprintf(stderr, "  rmxn_print, %d x %d\n", m, n1);
        rmxn_print(stderr, m, n1, Amn1);
        fputc('\n', stderr);
        
        fprintf(stderr, "  rmxn_gen_print, %d x %d\n", m, n1);
        rmxn_gen_print(stderr, m, n1, Amn1, "%6.3f", "  <\n    ", ",\n    ", "\n  >\n", "{ ", " ' ", " }");
        fputc('\n', stderr);
        
        fprintf(stderr, "  rmxn_gen_print2, %d x %d, %d x %d\n", m, n1, m, n2);
        rmxn_gen_print2(stderr, m, n1, Amn1, n2, Amn2, "%6.3f", "  <\n    ", ",\n    ", "\n  >\n", "{ ", " ' ", " }", " | ");
        fputc('\n', stderr);
        
        fprintf(stderr, "  rmxn_gen_print3, %d x %d, %d x %d, %d x %d\n", m, n1, m, n2, m, n3);
        rmxn_gen_print3(stderr, m, n1, Amn1, n2, Amn2, n3, Amn3, "%6.3f", "  <\n    ", ";\n    ", "\n  >\n", "{ ", " ' ", " }", " | ");
        fputc('\n', stderr);
        
        free(Amn1); free(Amn2); free(Amn3);
      }
  }
    
void test_rmxn_canonical_simplex_and_measures(int32_t n, bool_t verbose)
  {
    /* TEST: void rmxn_canonical_simplex(int32_t d, int32_t n, double V[]); */
    /* TEST: double rmxn_canonical_simplex_radius(int32_t d); */
    /* TEST: double rmxn_canonical_simplex_subradius(int32_t d, int32_t k); */
    /* TEST: double rmxn_canonical_simplex_edge(int32_t d); */
    /* TEST: double rmxn_canonical_simplex_height(int32_t d); */
    /* TEST: double rmxn_canonical_simplex_measure(int32_t d); */

    if (verbose) { fprintf(stderr, "--- rmxn_canonical_simplex_{radius,subradius,edge,height,measure} ---\n"); }
    { int32_t d = n-1;
      double V[(d+1)*n];
      rmxn_canonical_simplex(d, n, V);
      double rExp = rmxn_canonical_simplex_radius(d);
      double iExp = rmxn_canonical_simplex_subradius(d, d-1);
      double sExp = rmxn_canonical_simplex_edge(d);
      double hExp = rmxn_canonical_simplex_height(d);
      double mExp = rmxn_canonical_simplex_measure(d);
      check_simplex(d, n, V, rExp, iExp, sExp, hExp, mExp);
    }
  }
    
void test_rmxn_throw_canonical_simplex(int32_t n, bool_t verbose)
  {
    /* TEST: void rmxn_throw_canonical_simplex(int32_t d, double x[]); */
    /* TEST: void rmxn_throw_canonical_simplex_ball(int32_t d, double x[]); */

    if (verbose) { fprintf(stderr, "--- rmxn_throw_canonical_simplex ---\n"); }
    { int32_t d = n-1;
      double x[d+1];
      double tol = 1.0e-12;
      rmxn_throw_canonical_simplex(d, x);
      /* Check sign, unit-sum: */
      double sObs = 0;
      for (int32_t i = 0; i <= d; i++)
        { demand(x[i] >= 0, "rmxn_throw_canonical_simplex error (1)"); 
          sObs += x[i]; 
        }
      double sExp = 1.0;
      rn_check_eps(sObs, sExp, tol, NO, NO, "rmxn_throw_canonical_simplex error (2)");
    }

    if (verbose) { fprintf(stderr, "--- rmxn_throw_canonical_simplex_ball ---\n"); }
    { int32_t d = n-1;
      double x[d+1];
      double tol = 1.0e-12;
      rmxn_throw_canonical_simplex_ball(d, x);
      /* Check unit-sum, radius: */
      double sObs = 0, r2 = 0;
      double c = 1.0/(d+1);
      for (int32_t i = 0; i <= d; i++)
        { sObs += x[i];
          double d = x[i] - c;
          r2 += d*d;
        }
      double sExp = 1.0;
      rn_check_eps(sObs, sExp, tol, NO, NO, "rmxn_throw_canonical_simplex_ball error (1)");
      double rObs = sqrt(r2);
      double rExp = rmxn_canonical_simplex_radius(d);
      demand(rObs <= rExp*(1 + tol), "rmxn_throw_canonical_simplex_ball error (2)");
    }
  }

void test_rmxn_regular_simplex_and_measures(int32_t n, bool_t verbose)
  {
    /* TEST: void rmxn_regular_simplex(int32_t n, double V[]); */
    /* TEST: double rmxn_regular_simplex_radius(int32_t n); */
    /* TEST: double rmxn_regular_simplex_subradius(int32_t n, int32_t k); */
    /* TEST: double rmxn_regular_simplex_edge(int32_t n); */
    /* TEST: double rmxn_regular_simplex_height(int32_t n); */
    /* TEST: double rmxn_regular_simplex_measure(int32_t n); */

    if (verbose) { fprintf(stderr, "--- rmxn_regular_simplex ---\n"); }
    { double V[(n+1)*n];
      rmxn_regular_simplex(n, V);
      double rExp = rmxn_regular_simplex_radius(n);
      double iExp = rmxn_regular_simplex_subradius(n, n-1);
      double sExp = rmxn_regular_simplex_edge(n);
      double hExp = rmxn_regular_simplex_height(n);
      double mExp = rmxn_regular_simplex_measure(n);
      check_simplex(n, n, V, rExp, iExp, sExp, hExp, mExp);
    }
  }
    
void test_rmxn_throw_ortho_spin(int32_t n, bool_t verbose)
  {
    /* TEST: void rmxn_throw_ortho(int32_t n, double M[]); */
    /* TEST: void rmxn_spin_rows(int32_t m, int32_t n, double A[], double M[]); */
    /* TEST: void rmxn_spin_cols(int32_t m, int32_t n, double A[], double M[]); */
    /* TEST: void rmxn_shift_rows(int32_t m, int32_t n, double A[], double v[], double M[]); */
    /* TEST: void rmxn_shift_cols(int32_t m, int32_t n, double A[], double v[], double M[]); */

    if (verbose) { fprintf(stderr, "--- rmxn_throw_ortho ---\n"); }
    { double M[n*n];
      rmxn_throw_ortho(n, M);
      check_ortho_matrix(n, M);
    }

    if (verbose) { fprintf(stderr, "--- rmxn_spin_rows ---\n"); }
    if (verbose) { fprintf(stderr, "!! warning: rmxn_spin_rows not tested\n"); }
    
    if (verbose) { fprintf(stderr, "--- rmxn_shift_rows ---\n"); }
    if (verbose) { fprintf(stderr, "!! warning: rmxn_shift_rows not tested\n"); }

    if (verbose) { fprintf(stderr, "--- rmxn_spin_cols ---\n"); }
    if (verbose) { fprintf(stderr, "!! warning: rmxn_spin_cols not tested\n"); }
    
    if (verbose) { fprintf(stderr, "--- rmxn_shift_cols ---\n"); }
    if (verbose) { fprintf(stderr, "!! warning: rmxn_shift_cols not tested\n"); }
  }

void test_rmxn_max_abs_elem(int32_t m, int32_t n, double Amn[], bool_t verbose)
  {
    /* TEST: double rmxn_max_abs_elem(int32_t m, int32_t n, double *A); */

    if (verbose) { fprintf(stderr, "--- rmxn_max_abs_elem ---\n"); }
    if (verbose) { fprintf(stderr, "!! warning: rmxn_max_abs_elem not tested\n"); }
  }
    
void test_rmxn_perturb_unif(int32_t m, int32_t n, double Amn[], bool_t verbose)
  {
    /* TEST: void rmxn_perturb_unif(int32_t m, int32_t n, double mag, double M[]); */

    if (verbose) { fprintf(stderr, "--- rmxn_perturb_unif ---\n"); }
    if (verbose) { fprintf(stderr, "!! warning: rmxn_perturb_unif not tested\n"); }
  }
    
void test_rmxn_transform_quadratic(int32_t m, int32_t n, double Amn[], bool_t verbose)
  {
    /* TEST: void rmxn_transform_quadratic(int32_t n, double E[], double e[], int32_t m, double U[], double F[], double f[]); */ 

    if (verbose) { fprintf(stderr, "--- rmxn_transform_quadratic ---\n"); }
    if (verbose) { fprintf(stderr, "!! warning: rmxn_transform_quadratic not tested\n"); }
  }

void throw_matrix(int32_t m, int32_t n, double *Amn)
  { for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++) 
          { Amn[n*i + j] = 2.0 * drandom() - 1.0; }
      }
  }

void throw_LT_matrix(int32_t m, double *Lmm)
  { for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < m; j++) 
          { Lmm[m*i + j] = (j <= i ? 2.0 * drandom() - 1.0 : 0.0); }
      }
  }

void check_simplex(int32_t d, int32_t n, double V[], double rExp, double iExp, double sExp, double hExp, double mExp)
  { if (n < 1) { return; }
    demand(d <= n, "simplex dimension {d} too big for space dimension {n}");
    double tol = 1.0e-12;
    double dMax = -INFINITY;
    double dMin = +INFINITY;
    double rMax = -INFINITY;
    double rMin = +INFINITY;
    /* Compute barycenter {b[0..n-1]}: */
    double b[n];
    for (int32_t k = 0; k < n; k++)
      { double s = 0;
        for (int32_t i = 0; i <= d; i++) { s += V[i*n + k]; }
        b[k] = s/(d+1); 
      }
    /* Compute circum-radius range {rMin,rMax}: */
    for (int32_t i = 0; i <= d; i++)
      { /* Compute distance from corner {i} to {b}: */
        double r2 = 0;
        for (int32_t k = 0; k < n; k++) { double rk = V[i*n + k] - b[k]; r2 += rk*rk; }
        double r = sqrt(r2);
        if (r < rMin) { rMin = r; }
        if (r > rMax) { rMax = r; }
        for (int32_t j = 0; j < i; j++)
          { /* Check distance from corner {i} to corner {j}: */
            double d2 = 0;
            for (int32_t k = 0; k < n; k++)
              { double dk = V[i*n + k] - V[j*n + k]; d2 += dk*dk; }
            double d = sqrt(d2);
            if (d < dMin) { dMin = d; }
            if (d > dMax) { dMax = d; }
          }
      }
    if (rMax - rMin > tol*rMax)
      { print_matrix(stderr, d+1, n, V);
        fprintf(stderr, "    rMin =       %22.16e\n", rMin);
        fprintf(stderr, "    rMax =       %22.16e\n", rMax);
        fprintf(stderr, "    tolerance =  %22.16e\n", tol);
        demand(FALSE, "test simplex has irregular radius");
      }
    if (dMax - dMin > tol*dMax)
      { print_matrix(stderr, d+1, n, V);
        fprintf(stderr, "    dMin =       %22.16e\n", dMin);
        fprintf(stderr, "    dMax =       %22.16e\n", dMax);
        fprintf(stderr, "    tolerance =  %22.16e\n", tol);
        demand(FALSE, "test simplex has irregular sides");
      }
    /* Check circumradius: */
    double rObs = 0.5*(rMin+rMax);
    if (fabs(rExp - rObs) > tol*rObs)
      { print_matrix(stderr, d+1, n, V);
        fprintf(stderr, "    rExp =       %22.16e\n", rExp);
        fprintf(stderr, "    rObs =       %22.16e\n", rObs);
        fprintf(stderr, "    tolerance =  %22.16e\n", tol);
        demand(FALSE, "simplex circum-radius mismatch");
      }
    /* Check edge length: */
    double sObs = (d == 0 ? M_SQRT2 : 0.5*(dMin+dMax));
    if (fabs(sExp - sObs) > tol*sObs)
      { print_matrix(stderr, d+1, n, V);
        fprintf(stderr, "    sExp =       %22.16e\n", sExp);
        fprintf(stderr, "    sObs =       %22.16e\n", sObs);
        fprintf(stderr, "    tolerance =  %22.16e\n", tol);
        demand(FALSE, "simplex edge length mismatch");
      }
    /* Check height: */
    double hObs = ((double)d+1)/((double)d)*rObs;
    if (fabs(hExp - hObs) > tol*hObs)
      { print_matrix(stderr, d+1, n, V);
        fprintf(stderr, "    hExp =       %22.16e\n", hExp);
        fprintf(stderr, "    hObs =       %22.16e\n", hObs);
        fprintf(stderr, "    tolerance =  %22.16e\n", tol);
        demand(FALSE, "simplex height mismatch");
      }
    /* Check inradius: */
    double iObs = hObs - rObs;
    if (fabs(iExp - iObs) > tol*iObs)
      { print_matrix(stderr, d+1, n, V);
        fprintf(stderr, "    iExp =       %22.16e\n", iExp);
        fprintf(stderr, "    iObs =       %22.16e\n", iObs);
        fprintf(stderr, "    tolerance =  %22.16e\n", tol);
        demand(FALSE, "simplex in-radius mismatch");
      }
    /* Check measure: */
    { 
      double det; /* Nominal determinant. */
      if (d == n)
        { /* Assume that it is a regular simplex in all {n} coords: */
          double D[n*n]; 
          for (int32_t i = 1; i <= d; i++)
            { /* Subtract corner 0 from corner {i}: */
              for (int32_t k = 0; k < n; k++) { D[(i-1)*n + k] = V[i*n + k] - V[0*n + k]; }
            }
          det = rmxn_det(n, D);
        }
      else
        { /* Require that simplex be the canonical simplex in the first {d+1} coords: */
          for (int32_t i = 0; i <= d; i++) 
            { for (int32_t j = 0; j < n; j++) 
                { double req = (i == j ? 1 : 0);
                  demand(V[i*n + j] == req, "simplex is not canonical");
                }
            }
          /* From formula: */
          det = sqrt(d+1);
        }
      /* Volume is determinant divided by {d!}: */
      double mObs = det;
      for (int32_t i = 1; i <= d; i++) { mObs /= i; }
      
      if (fabs(mExp - mObs) > tol*mObs)
        { print_matrix(stderr, d+1, n, V);
          fprintf(stderr, "    mExp =       %22.16e\n", mExp);
          fprintf(stderr, "    mObs =       %22.16e\n", mObs);
          fprintf(stderr, "    tolerance =  %22.16e\n", tol);
          demand(FALSE, "simplex measure mismatch");
        }
    }
  }

void check_ortho_matrix(int32_t n, double M[])
  { 
    double tol = 1.0e-12;
    for (int32_t i0 = 0; i0 < n; i0++)
      { for (int32_t i1 = 0; i1 <= i0; i1++)
          { /* Compute dot product {sObs} of rows {i0} and {i1}: */
            double sObs = 0;
            for (int32_t j = 0; j < n; j++) { sObs += M[i0*n + j]*M[i1*n + j]; }
            /* Check dot product against expected value {sExp}: */
            double sExp = (i0 == i1 ? 1 : 0);
            rn_check_eps(sExp, sObs, tol, &i0, &i1, "rmxn_throw_ortho error");
          }
      }
  }

void print_matrix(FILE *wr, int32_t m, int32_t n, double *Amn)
  { fprintf(wr, "%d × %d matrix\n", m, n);
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { fprintf(wr, " %22.16e", Amn[i*n + j]); }
        fprintf(wr, "\n");
      }
    fprintf(wr, "\n");
  }
