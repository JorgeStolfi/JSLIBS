/* test_r3 --- test program for r3.h, r3x3.h  */
/* Last edited on 2024-11-20 13:19:24 by stolfi */

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>
#include <jsrandom.h>

#include <r3.h>
#include <r3x3.h>
#include <r3_hedron.h>
#include <r3_bezier.h>
#include <r3_path.h>
#include <r3_motion.h>

#include <rn_test_tools.h>
#include <rmxn_test_tools.h>

#define N 3
#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);

void test_r3(bool_t verbose);
void test_r3x3(bool_t verbose);
void test_r3_hedron(bool_t verbose);
void test_r3_bezier(bool_t verbose);
void test_r3_path(bool_t verbose);
void test_r3_motion(bool_t verbose);

void test_r3_zero(bool_t verbose);
void test_r3_all(bool_t verbose);
void test_r3_axis(bool_t verbose);
void test_r3_throw_cube(bool_t verbose);
void test_r3_throw_ball(bool_t verbose);
void test_r3_throw_normal(bool_t verbose);
void test_r3_throw_dir(bool_t verbose);
void test_r3_throw_ortho(bool_t verbose);
void test_r3_throw_ortho_pair(bool_t verbose);
void test_r3_add(bool_t verbose);
void test_r3_sub(bool_t verbose);
void test_r3_neg(bool_t verbose);
void test_r3_scale(bool_t verbose);
void test_r3_mix(bool_t verbose);
void test_r3_mix_in(bool_t verbose);
void test_r3_weigh(bool_t verbose);
void test_r3_unweigh(bool_t verbose);
void test_r3_rot_axis(bool_t verbose);
void test_r3_norm__r3_norm_sqr__r3_L_inf_norm(bool_t verbose);
void test_r3_dist__r3_dist_sqr__r3_L_inf_dist(bool_t verbose);
void test_r3_dir__r3_L_inf_dir(bool_t verbose);
void test_r3_dot__r3_cos__r3_sin__r3_angle(bool_t verbose);
void test_r3_cross(bool_t verbose);
void test_r3_det(bool_t verbose);
void test_r3_decomp(bool_t verbose);
void test_r3_print(bool_t verbose);
void test_r3_hedron_tetra_vertices(bool_t verbose);
void test_r3_hedron_octa_vertices(bool_t verbose);
void test_r3_hedron_hexa_vertices(bool_t verbose);
void test_r3_hedron_icosa_vertices(bool_t verbose);
void test_r3_hedron_dodeca_vertices(bool_t verbose);

void test_r3x3_size_allocation(bool_t verbose);
void test_r3x3_indexing_adressing(bool_t verbose);
void test_r3x3_zero__r3x3_ident(bool_t verbose);
void test_r3x3_throw(bool_t verbose);
void test_r3x3_get_row__r3x3_get_col__r3x3_set_row__r3x3_set_col(bool_t verbose);
void test_r3x3_map_row__r3x3_map_col(bool_t verbose);
void test_r3x3_scale(bool_t verbose);
void test_r3x3_mul(bool_t verbose);
void test_r3x3_mul_tr(bool_t verbose);
void test_r3x3_tr_mul(bool_t verbose);
void test_r3x3_transp(bool_t verbose);
void test_r3x3_det(bool_t verbose);
void test_r3x3_inv(bool_t verbose);
void test_r3x3_norm__r3x3_norm_sqr__r3x3_mod_norm_sqr(bool_t verbose);
void test_r3x3_normalize(bool_t verbose);
void test_r3x3_diff_sqr(bool_t verbose);
void test_r3x3_u_v_rotation(bool_t verbose);
void test_r3x3_throw_rotation(bool_t verbose);
void test_r3x3_print(bool_t verbose);

void tr3_throw_matrix(r3x3_t *M);
void tr3_print_matrix(char *name, r3x3_t *A);

void tr3_check_num_eps(char *name, double x, double y, double eps, char *msg);
  /* If {x} and {y} differ by more than {eps}, prints {name}, {x}, {y}, and {msg}, and stops. */

void tr3_check_regular_polyhedron(char *func, double R, double L, uint32_t n, r3_t p[], uint32_t deg);
  /* Check that {p[0..n-1]} are the vertices of a regular 
    polyhedron with radius {R}, side {L}, and vertex degree {deg}. */

int32_t main (int32_t argc, char **argv)
  { srand(1993);
    srandom(1993);

    for (uint32_t i = 0;  i < 100; i++) test_r3(i < 3);
    for (uint32_t i = 0;  i < 100; i++) test_r3x3(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_r3(bool_t verbose)
  { if (verbose)
      { fprintf(stderr,
          "sizeof(r3_t) = %lu  %d*sizeof(double) = %lu\n",
          sizeof(r3_t), N, N*sizeof(double)
        );
      }
      
    test_r3_zero(verbose);
    test_r3_all(verbose);
    test_r3_axis(verbose);
    test_r3_throw_cube(verbose);
    test_r3_throw_ball(verbose);
    test_r3_throw_normal(verbose);
    test_r3_throw_dir(verbose);
    test_r3_throw_ortho(verbose);
    test_r3_throw_ortho_pair(verbose);
    test_r3_add(verbose);
    test_r3_sub(verbose);
    test_r3_neg(verbose);
    test_r3_scale(verbose);
    test_r3_mix(verbose);
    test_r3_mix_in(verbose);
    test_r3_weigh(verbose);
    test_r3_unweigh(verbose);
    test_r3_rot_axis(verbose);
    test_r3_norm__r3_norm_sqr__r3_L_inf_norm(verbose);
    test_r3_dist__r3_dist_sqr__r3_L_inf_dist(verbose);
    test_r3_dir__r3_L_inf_dir(verbose);
    test_r3_dot__r3_cos__r3_sin__r3_angle(verbose);
    test_r3_cross(verbose);
    test_r3_det(verbose);
    test_r3_decomp(verbose);
    test_r3_print(verbose);
    
    if (verbose)
      { 
        fprintf(stderr, "!! r3_barycenter NOT TESTED\n");
        fprintf(stderr, "!! r3_bbox NOT TESTED\n");
        fprintf(stderr, "!! r3_circumcenter NOT TESTED\n");
        fprintf(stderr, "!! r3_eq NOT TESTED\n");
        fprintf(stderr, "!! r3_gen_print NOT TESTED\n");
        fprintf(stderr, "!! r3_insphere NOT TESTED\n");
        fprintf(stderr, "!! r3_is_finite NOT TESTED\n");
        fprintf(stderr, "!! r3_throw_normal NOT TESTED\n");
        fprintf(stderr, "!! r3_orient NOT TESTED\n");
      }
  }

void test_r3_zero(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_zero ---\n"); }
    r3_t a;
    r3_zero(&a);
    for (uint32_t i = 0;  i < N; i++)
      { rn_test_tools_check_eq(a.c[i],0.0, NO, NO, "r3_zero error"); }
  }

void test_r3_all(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_all ---\n"); }
    r3_t a;
    r3_all(3.14, &a);
    for (uint32_t i = 0;  i < N; i++)
      { rn_test_tools_check_eq(a.c[i],3.14, NO, NO, "r3_all error"); }
  }

void test_r3_axis(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_axis ---\n"); }
    r3_t a;
    for (uint32_t k = 0;  k < N; k++)
      { r3_axis(k, &a);
        for (uint32_t i = 0;  i < N; i++)
          { rn_test_tools_check_eq(a.c[i],(i == k ? 1.0 : 0.0), NO, NO, "r3_axis error"); }
      }
  }

void test_r3_throw_cube(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_throw_cube ---\n"); }
    r3_t a;
    r3_throw_cube(&a);
    /* Check variation: */
    rn_test_tools_check_all_different(N, a.c, "r3_throw_cube probable error (1)");
    for (uint32_t i = 0;  i < N; i++)
      { /* Check whether there are more than 8 nonzero bits: */
        double vv = a.c[i]*256.0;
        affirm(vv != floor(vv), "r3_throw_cube probable error (3)"); 
        affirm((a.c[i] > -1.0) && (a.c[i] < 1.0), "r3_throw error (2)"); 
      }
  }

void test_r3_throw_ball(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_throw_ball ---\n"); }
    /* Should check uniformity... */
    r3_t a;
    r3_throw_ball(&a);
    /* Check variation: */
    rn_test_tools_check_all_different(N, a.c, "r3_throw_ball probable error (1)");
    /* Check whether the norm is at most 1: */
    double rr = 0;
    for (uint32_t i = 0;  i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    rr = sqrt(rr);
    demand(rr <= 1 + 0.000000001*rr, "r3_throw_ball error (2)");
  }

void test_r3_throw_normal(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_throw_normal ---\n"); }
    /* Should check Gaussian distribution... */
    r3_t a;
    r3_throw_normal(&a);
    /* Check variation: */
    rn_test_tools_check_all_different(N, a.c, "r3_throw_normal probable error (1)");
  }

void test_r3_throw_dir(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_throw_dir ---\n"); }
    /* Should check uniformity... */
    r3_t a;
    r3_throw_dir(&a);
    /* Check variation: */
    rn_test_tools_check_all_different(N, a.c, "r3_throw_dir probable error (1)");
    /* Check whether the norm is 1: */
    double rr = 0;
    for (uint32_t i = 0;  i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    rr = sqrt(rr);
    rn_test_tools_check_eps(1,rr,0.000000001*rr, NO, NO, "r3_throw_dir error (2)");
  }

void test_r3_throw_ortho(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_throw_ortho ---\n"); }
    /* !!! Should check uniformity !!! */
    r3_t a, d;
    uint32_t nrand = 4; /* Number of random trials. */
    for (uint32_t trial = 0;  trial < N+1+nrand; trial++)
      { if (trial < N)
          { uint32_t i = trial;
            r3_axis(i, &a);
          }
        else if (trial == N)
          { /* Test with a zero vector: */
            r3_zero(&a);
          }
        else
          { /* Pick a random nonzero vector: */
            r3_throw_ball(&a);
            r3_scale(0.1 + 5*drandom(), &a, &a);
          }
         
        /* Check same length as {a}: */
        double ma = r3_norm(&a);
        double mo = r3_throw_ortho(&a, &d);
        double md = r3_norm(&d);
        rn_test_tools_check_eps(ma, mo, 0.00000001*ma + 1.0e-15, NO, NO, "r3_throw_ortho error (1)");
        rn_test_tools_check_eps(ma, md, 0.00000001*ma + 1.0e-15, NO, NO, "r3_throw_ortho error (2)");
        /* Check orthogonality: */
        double dot = 0;
        for (uint32_t i = 0;  i < N; i++) 
          { double di = d.c[i], ai = a.c[i]; dot += di*ai; }
        rn_test_tools_check_eps(dot, 0.0, 0.00000001*ma, NO, NO, "r3_throw_ortho error (3)");
        if (ma > 1.0e-12)
          { /* Check variation: */
            rn_test_tools_check_all_different(N, d.c, "r3_throw_ortho probable error (4)");
          }
      }
  }

void test_r3_throw_ortho_pair(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_throw_ortho_pair ---\n"); }
    /* !!! Should check uniformity !!! */
    r3_t a, d, e;
    uint32_t nrand = 4; /* Number of random trials. */
    for (uint32_t trial = 0;  trial < N+1+nrand; trial++)
      { if (trial < N)
          { uint32_t i = trial;
            r3_axis(i, &a);
          }
        else if (trial == N)
          { /* Test with a zero vector: */
            r3_zero(&a);
          }
        else
          { /* Pick a random nonzero vector: */
            r3_throw_ball(&a);
            r3_scale(0.1 + 5*drandom(), &a, &a);
          }
         
        /* Check same length as {a}: */
        double ma = r3_norm(&a);
        double mo = r3_throw_ortho_pair(&a, &d, &e);
        double md = r3_norm(&d);
        double me = r3_norm(&e);
        rn_test_tools_check_eps(ma, mo, 0.00000001*ma + 1.0e-15, NO, NO, "r3_throw_ortho_pair error (1)");
        rn_test_tools_check_eps(ma, md, 0.00000001*ma + 1.0e-15, NO, NO, "r3_throw_ortho_pair error (2)");
        rn_test_tools_check_eps(ma, me, 0.00000001*ma + 1.0e-15, NO, NO, "r3_throw_ortho_pair error (3)");
        /* Check orthogonality: */
        double addot = 0, aedot = 0, dedot = 0;
        for (uint32_t i = 0;  i < N; i++) 
          { double ai = a.c[i], di = d.c[i], ei = e.c[i];
            addot += ai*di; aedot += ai*ei; dedot += di*ei;
          }
        rn_test_tools_check_eps(addot, 0.0, 0.00000001*ma, NO, NO, "r3_throw_ortho_pair error (4)");
        rn_test_tools_check_eps(aedot, 0.0, 0.00000001*ma, NO, NO, "r3_throw_ortho_pair error (5)");
        rn_test_tools_check_eps(dedot, 0.0, 0.00000001*ma, NO, NO, "r3_throw_ortho_pair error (6)");
        if (ma > 1.0e-12)
          { /* Check variation: */
            rn_test_tools_check_all_different(N, d.c, "r3_throw_ortho_pair probable error (7)");
            rn_test_tools_check_all_different(N, e.c, "r3_throw_ortho_pair probable error (8)");
          }
      }
  }

void test_r3_add(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_add ---\n"); }
    r3_t a, b, d;
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_add(&a, &b, &d);
    for (uint32_t i = 0;  i < N; i++)
      { rn_test_tools_check_eq(d.c[i],a.c[i] + b.c[i], NO, NO, "r3_add error"); }
  }
void test_r3_sub(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_sub ---\n"); }
    r3_t a, b, d;
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_sub(&a, &b, &d);
    for (uint32_t i = 0;  i < N; i++)
      { rn_test_tools_check_eq(d.c[i],a.c[i] - b.c[i], NO, NO, "r3_sub error"); }
  }

void test_r3_neg(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_neg ---\n"); }
    r3_t a, d;
    r3_throw_cube(&a);
    r3_neg(&a, &d);
    for (uint32_t i = 0;  i < N; i++)
      { rn_test_tools_check_eq(d.c[i],- a.c[i], NO, NO, "r3_neg error"); }
  }

void test_r3_scale(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_scale ---\n"); }
    r3_t a, d;
    double s = drandom();
    r3_throw_cube(&a);
    r3_scale(s, &a, &d);
    for (uint32_t i = 0;  i < N; i++)
      { double zi = s*a.c[i];
        rn_test_tools_check_eq(d.c[i],zi, NO, NO, "r3_scale error (1)");
      }
  }

void test_r3_mix(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_mix ---\n"); }
    r3_t a, b, d;
    double s = drandom();
    double t = drandom();
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_mix(s, &a, t, &b, &d);
    for (uint32_t i = 0;  i < N; i++)
      { double ddi = s * a.c[i] + t * b.c[i];
        rn_test_tools_check_eq(d.c[i],ddi, NO, NO, "r3_mix error");
      }
  }

void test_r3_mix_in(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_mix_in ---\n"); }
    r3_t a, b, d;
    double s = drandom();
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    d = b;
    r3_mix_in(s, &a, &d);
    for (uint32_t i = 0;  i < N; i++)
      { double ddi = b.c[i] + s * a.c[i];
        rn_test_tools_check_eq(d.c[i],ddi, NO, NO, "r3_mix_in error");
      }
  }

void test_r3_weigh(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_weigh ---\n"); }
    r3_t a, b, d;
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_weigh(&a, &b, &d);
    for (uint32_t i = 0;  i < N; i++)
      { double ddi = a.c[i] * b.c[i];
        rn_test_tools_check_eq(d.c[i],ddi, NO, NO, "r3_weigh error");
      }
  }

void test_r3_unweigh(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_unweigh ---\n"); }
    r3_t a, b, d;
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_unweigh(&a, &b, &d);
    for (uint32_t i = 0;  i < N; i++)
      { double ddi = a.c[i] / b.c[i];
        rn_test_tools_check_eq(d.c[i],ddi, NO, NO, "r3_unweigh error");
      }
  }

void test_r3_rot_axis(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_rot_axis ---\n"); }
    r3_t a, d;
    r3_throw_cube(&a);
    uint32_t i = uint32_abrandom(0, N-1);
    uint32_t j = uint32_abrandom(0, N-2); if (j >= i) { j++; }
    double ang = 2.1*M_PI*drandom();
    r3_rot_axis(&a, i, j, ang, &d);
    rn_test_tools_check_rot_axis(N, a.c, i, j, ang, d.c, "r3_rot_axis error");
  }

void test_r3_norm__r3_norm_sqr__r3_L_inf_norm(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_norm, r3_norm_sqr, r3_L_inf_norm ---\n"); }
    r3_t a;
    r3_throw_cube(&a);
    double r = r3_norm(&a);
    double s = r3_norm_sqr(&a);
    double t = r3_L_inf_norm(&a);
    double ss = 0.0;
    double tt = 0.0;
    for (uint32_t i = 0;  i < N; i++)
      { double ai = fabs(a.c[i]);
        ss += ai*ai; 
        if (ai > tt) { tt = ai; }
      }
    double rr = sqrt(ss);
    rn_test_tools_check_eps(r,rr,0.000000001 * rr, NO, NO, "r3_norm error");
    rn_test_tools_check_eps(s,ss,0.000000001 * ss, NO, NO, "r3_norm_sqr error");
    rn_test_tools_check_eq(t,tt, NO, NO, "r3_L_inf_norm error");
  }

void test_r3_dist__r3_dist_sqr__r3_L_inf_dist(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_dist, r3_dist_sqr, r3_L_inf_dist ---\n"); }
    r3_t a, b;
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    double r = r3_dist(&a, &b);
    double s = r3_dist_sqr(&a, &b);
    double t = r3_L_inf_dist(&a, &b);
    double ss = 0.0;
    double tt = 0.0;
    for (uint32_t i = 0;  i < N; i++)
      { double di = fabs(a.c[i] - b.c[i]);
        ss += di*di; 
        if (di > tt) { tt = di; }
      }
    double rr = sqrt(ss);
    rn_test_tools_check_eps(r,rr,0.000000001 * rr, NO, NO, "r3_dist error");
    rn_test_tools_check_eps(s,ss,0.000000001 * ss, NO, NO, "r3_dist_sqr error");
    rn_test_tools_check_eq(t,tt, NO, NO, "r3_L_inf_dist error");
  }

void test_r3_dir__r3_L_inf_dir(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_dir, r3_L_inf_dir ---\n"); }
    r3_t a, b, d; 
    r3_throw_cube(&a);
    double r = r3_dir(&a, &b);
    double rr = r3_norm(&a);
    rn_test_tools_check_eps(r, rr, 0.000000001*rr, NO, NO, "r3_dir error (1)");
    
    double s = r3_L_inf_dir(&a, &d);
    double ss = r3_L_inf_norm(&a);
    rn_test_tools_check_eq(s, ss, NO, NO, "r3_L_inf_dir error (1)");
        
    for (uint32_t i = 0;  i < N; i++)
      { rn_test_tools_check_eps(b.c[i],a.c[i]/rr, 0.000000001*rr, NO, NO, "r3_dir error (2)");
        rn_test_tools_check_eps(d.c[i],a.c[i]/ss, 0.000000001*ss, NO, NO, "r3_L_inf_dir error (2)");
      }
  }

void test_r3_dot__r3_cos__r3_sin__r3_angle(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_dot, r3_cos, r3_sin, r3_angle ---\n"); }
    r3_t a, b, d; 
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    double r = r3_dot(&a, &b);
    double S = r3_sin(&a, &b);
    double C = r3_cos(&a, &b);
    double A = r3_angle(&a, &b);
    double mag = sqrt(r3_dot(&a,&a)*r3_dot(&b,&b));
    double rr = 0.0;
    for (uint32_t i = 0;  i < N; i++) { rr += a.c[i]*b.c[i]; }
    double CC = rr/(r3_norm(&a)*r3_norm(&b));
    rn_test_tools_check_eps(r,rr,0.000000001 * mag, NO, NO, "r3_dot error (1)");
    rn_test_tools_check_eps(C,CC,0.000000001, NO, NO, "r3_cos error (1)");
    d = a;
    r3_mix_in(-rr/r3_norm_sqr(&b), &b, &d);
    double SS = r3_norm(&d)/r3_norm(&a);
    rn_test_tools_check_eps(S,SS,0.000000001, NO, NO, "r3_sin error (1)");
    double AA = atan2(SS, CC);
    rn_test_tools_check_eps(A,AA,0.000000001, NO, NO, "r3_angle error (1)");
    for (uint32_t i = 0;  i < N; i++)
      { r3_axis(i, &a);
        for (uint32_t j = 0;  j < N; j++)
          { r3_axis(j, &b);
            double r = r3_dot(&a, &b);
            double s = r3_sin(&a, &b);
            double t = r3_cos(&a, &b);
            double rr = (i == j ? 1.0 : 0.0);
            rn_test_tools_check_eq(r,rr, NO, NO, "r3_dot error (2)");
            rn_test_tools_check_eq(t,rr, NO, NO, "r3_dot error (2)");
            rn_test_tools_check_eq(s,1.0 - rr, NO, NO, "r3_dot error (2)");
          }
      }
  }

void test_r3_cross(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_cross ---\n"); }
    r3_t a, b, d, e;
    /* Test on basis vectors: */
    for (uint32_t i = 0;  i < N; i++)
      { uint32_t i0 = (i + 0) % N;
        uint32_t i1 = (i + 1) % N;
        uint32_t i2 = (i + 2) % N;
        r3_axis(i0, &a);
        r3_axis(i1, &b);
        r3_cross(&a, &b, &d);
        r3_axis(i2, &e);
        for (uint32_t p = 0;  p < N; p++)
          { double ep = e.c[p];
            rn_test_tools_check_eq(d.c[p],ep, NO, NO, "r3_cross error (x)");
          }
      }
    /* Test on random vectors: */
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_cross(&a, &b, &d);
    double mag = r3_norm(&a)*r3_norm(&b);
    double r = r3_dot(&a, &d);
    rn_test_tools_check_eps(r, 0.0, 0.00000001*mag*r3_norm(&a), NO, NO, "r3_cross error (1)");
    double s = r3_dot(&b, &d);
    rn_test_tools_check_eps(s, 0.0, 0.00000001*mag*r3_norm(&b), NO, NO, "r3_cross error (2)");
  }

void test_r3_det(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_det ---\n"); }
    r3_t a, b, c, e;
    /* Test on basis vectors: */
    for (uint32_t i = 0;  i < N; i++)
      { uint32_t i0 = (i + 0) % N;
        uint32_t i1 = (i + 1) % N;
        uint32_t i2 = (i + 2) % N;
        r3_axis(i0, &a);
        r3_axis(i1, &b);
        r3_axis(i2, &c);
        double r = r3_det(&a, &b, &c);
        rn_test_tools_check_eq(r,1.0, NO, NO, "r3_det error (2)");
      }
    /* Test on random vectors: */
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_throw_cube(&c);
    double r = r3_det(&a, &b, &c);
    r3_cross(&a, &b, &e);
    double rr = r3_dot(&e, &c);
    double mag = r3_norm(&a)*r3_norm(&b)*r3_norm(&c);
    rn_test_tools_check_eps(r,rr, 0.00000001*mag, NO, NO, "r3_det error (1)");
  }

void test_r3_decomp(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_decomp ---\n"); } 
    r3_t a, b, c, para, perp;
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    double r = r3_decomp(&a, &b, &para, &perp);
    double rr = r3_dot(&a, &b)/r3_norm_sqr(&b);  
    rn_test_tools_check_eps(r,rr, 0.000000001*(fabs(r) + fabs(rr)), NO, NO, "r3_decomp error (1)");
    
    r3_add(&para, &perp, &c);
    double u = r3_dist(&a, &c);
    affirm (u <= 0.000000001*r3_norm(&a), "r3_decomp error (2)");

    double s = r3_dot(&perp, &b);
    rn_test_tools_check_eps(s, 0.0, 0.000000001*r3_norm(&b), NO, NO, "r3_decomp error (3)");

    double t = r3_dot(&para, &perp);
    rn_test_tools_check_eps(t, 0.0, 0.000000001*r3_norm(&a), NO, NO, "r3_decomp error (4)");
  }

void test_r3_print(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_print ---\n"); }
    r3_t a;
    if (verbose)
      { r3_throw_cube (&a);
        fprintf(stderr, "a = ");
        r3_print(stderr, &a);
        fputc('\n', stderr);
      }
  }

void test_r3_hedron_tetra_vertices(bool_t verbose)
  { if (verbose) { fprintf(stderr, "--- r3_hedron_tetra_vertices ---\n"); }
    r3_t pv[4];
    double R = drandom();
    r3_hedron_tetra_vertices(R, 4, pv);
    double tetra_L = R*sqrt(8.0/3.0);
    tr3_check_regular_polyhedron("r3_hedron_tetra_vertices", R, tetra_L, 4, pv, 3);
  }

void test_r3_hedron_octa_vertices(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_hedron_octa_vertices ---\n"); }
    r3_t pv[6];
    double R = drandom();
    r3_hedron_octa_vertices(R, 6, pv);
    double octa_L = R*sqrt(2.0);
    tr3_check_regular_polyhedron("r3_hedron_octa_vertices", R, octa_L, 6, pv, 4);
  }

void test_r3_hedron_hexa_vertices(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_hedron_hexa_vertices ---\n"); }
    r3_t pv[8];
    double R = drandom();
    r3_hedron_hexa_vertices(R, 8, pv);
    double hexa_L = 2*R/sqrt(3.0);
    tr3_check_regular_polyhedron("r3_hedron_hexa_vertices", R, hexa_L, 8, pv, 3);
  }

void test_r3_hedron_icosa_vertices(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_hedron_icosa_vertices ---\n"); }
    r3_t pv[12];
    double R = drandom();
    r3_hedron_icosa_vertices(R, 12, pv);
    double icosa_L = R*sqrt(2 - 2/sqrt(5));
    tr3_check_regular_polyhedron("r3_hedron_icosa_vertices", R, icosa_L, 12, pv, 5);
  }

void test_r3_hedron_dodeca_vertices(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3_hedron_dodeca_vertices ---\n"); }
    r3_t pv[20];
    double R = drandom();
    r3_hedron_dodeca_vertices(R, 20, pv);
    double dodeca_L = R*(sqrt(5)-1)/sqrt(3);
    tr3_check_regular_polyhedron("r3_hedron_dodeca_vertices", R, dodeca_L, 20, pv, 3);
   }

/* TESTING {r3x3.h} FUNCTIONS */

void test_r3x3(bool_t verbose)
  {
    test_r3x3_size_allocation(verbose);
    test_r3x3_indexing_adressing(verbose);
    test_r3x3_zero__r3x3_ident(verbose);
    test_r3x3_throw(verbose);
    test_r3x3_get_row__r3x3_get_col__r3x3_set_row__r3x3_set_col(verbose);
    test_r3x3_map_row__r3x3_map_col(verbose);
    test_r3x3_scale(verbose);
    test_r3x3_mul(verbose);
    test_r3x3_mul_tr(verbose);
    test_r3x3_tr_mul(verbose);
    test_r3x3_transp(verbose);
    test_r3x3_det(verbose);
    test_r3x3_inv(verbose);
    test_r3x3_norm__r3x3_norm_sqr__r3x3_mod_norm_sqr(verbose);
    test_r3x3_normalize(verbose);
    test_r3x3_diff_sqr(verbose);
    test_r3x3_u_v_rotation(verbose);
    test_r3x3_throw_rotation(verbose);
    test_r3x3_print(verbose);

    if (verbose)
      { 
        fprintf(stderr, "!! r3x3_add NOT TESTED\n");
        fprintf(stderr, "!! r3x3_sub NOT TESTED\n");
        fprintf(stderr, "!! r3x3_neg NOT TESTED\n");
        fprintf(stderr, "!! r3x3_mix NOT TESTED\n");
        fprintf(stderr, "!! r3x3_adj NOT TESTED\n");
        fprintf(stderr, "!! r3x3_is_unif_scaling NOT TESTED\n");
        fprintf(stderr, "!! r3x3_from_cols NOT TESTED\n");
        fprintf(stderr, "!! r3x3_from_rows NOT TESTED\n");
        fprintf(stderr, "!! r3x3_gen_print NOT TESTED\n");
        
        fprintf(stderr, "!! r3x3_L_inf_norm NOT TESTED\n");
        fprintf(stderr, "!! r3x3_L_inf_normalize NOT TESTED\n");
      }
  }

void test_r3x3_size_allocation(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_t size and allocation ---\n"); }
    r3x3_t A, B, C;
    if (verbose)
      { fprintf(stderr,
          "sizeof(r3x3_t) = %lu  %d*%d*sizeof(double) = %lu\n",
          sizeof(r3x3_t), N, N, N*N*sizeof(double)
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

void test_r3x3_indexing_adressing(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_t indexing and addressing ---\n"); }
    r3x3_t A;
    for (uint32_t i = 0;  i < N; i++)
      for (uint32_t j = 0;  j < N; j++)
        { double *Aij = &(A.c[i][j]); 
          affirm(Aij == ((double *)&A)+(N*i)+j, "r3x3_t indexing error");
        }
  }

void test_r3x3_zero__r3x3_ident(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_zero, r3x3_ident ---\n"); }
    r3x3_t A, B;
    r3x3_zero(&A);
    r3x3_ident(&B);
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < N; j++)
          { rn_test_tools_check_eq(A.c[i][j],0.0, NO, NO, "r3x3_zero error"); 
            double Iij = (i == j ? 1.0 : 0.0);
            rn_test_tools_check_eq(B.c[i][j],Iij, NO, NO, "r3x3_ident error");
          }
      }
  }

void test_r3x3_throw(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_throw ---\n"); }
    r3x3_t A;
    for (sign_t sgn = +1; sgn >= -1; sgn--)
      { r3x3_throw(&A, sgn);
        rmxn_test_tools_check_all_different(N, N, &(A.c[0][0]), "r3x3_throw probable error");
        if (sgn != 0)
          { double det = r3x3_det(&A);
            demand(det*sgn > 0, "r3x3_throw det sign error");
          }
      }
  }

void test_r3x3_get_row__r3x3_get_col__r3x3_set_row__r3x3_set_col(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_get_row, r3x3_set_row, r3x3_get_col, r3x3_set_col ---\n"); }
    r3x3_t A;
    r3_t a;
    r3x3_throw(&A, 0);
    for (uint32_t dir = 0;  dir < 2; dir++) 
      { /* {dir} is 0 for row, 1 for col. */
        for (uint32_t i = 0;  i < N; i++)
          { /* Check {r3x3_get_row,r3x3_get_col}: */
            r3_throw_cube(&a);
            if (dir == 0) { r3x3_get_row(&A, i, &a); } else { r3x3_get_col(&A, i, &a); }
            for (uint32_t j = 0;  j < N; j++)
              { double vj = (dir == 0 ? A.c[i][j] : A.c[j][i]);
                affirm(vj = a.c[j], "r3x3_get_row/r3x3_get_col error");
              }
            /* Check {r3x3_set_row,r3x3_set_col}: */
            r3_throw_cube(&a);
            if (dir == 0) { r3x3_set_row(&A, i, &a); } else { r3x3_set_col(&A, i, &a); }
            for (uint32_t j = 0;  j < N; j++)
              { double vj = (dir == 0 ? A.c[i][j] : A.c[j][i]);
                affirm(vj = a.c[j], "r3x3_set_row/r3x3_set_col error");
              }
          }
      }
  }

void test_r3x3_map_row__r3x3_map_col(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_map_row, r3x3_map_col ---\n"); }
    r3x3_t A;
    r3_t a, b, c, bb, cc;
    r3x3_throw(&A, 0);
    r3_throw_cube(&a);
    r3x3_map_row(&a, &A, &b);
    r3x3_map_col(&A, &a, &c);
    r3_zero(&bb);
    r3_zero(&cc);
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < N; j++)
          { bb.c[j] += a.c[i] * A.c[i][j];
            cc.c[i] += A.c[i][j] * a.c[j];
          }
      }
    double r = r3_dist(&b, &bb);
    affirm(r < 0.000000001 * r3_norm(&bb), "r3_map_row error");
    double s = r3_dist(&c, &cc);
    affirm(s < 0.000000001 * r3_norm(&cc), "r3_map_col error");
  }

void test_r3x3_scale(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_scale ---\n"); }
    r3x3_t A, C;
    r3x3_throw(&A, 0);
    double r = drandom();
    r3x3_scale(r, &A, &C);
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < N; j++)
          { double sel = r * A.c[i][j];
            rn_test_tools_check_eps(C.c[i][j],sel,0.000000001 * fabs(sel), NO, NO,
              "r3x3_scale error"
            );
          }
      }
  }

void test_r3x3_mul(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_mul ---\n"); }
    r3x3_t A, B, C;
    r3x3_throw(&A, 0);
    r3x3_throw(&B, 0);
    r3x3_mul(&A, &B, &C);
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < N; j++)
          { double sum = 0.0;
            for (uint32_t k = 0;  k < N; k++) { sum += A.c[i][k]*B.c[k][j]; }
            rn_test_tools_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
              "r3x3_mul error"
            );
          }
      }
  }

void test_r3x3_mul_tr(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_mul_tr ---\n"); }
    r3x3_t A, B, C;
    r3x3_throw(&A, 0);
    r3x3_throw(&B, 0);
    r3x3_mul_tr(&A, &B, &C);
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < N; j++)
          { double sum = 0.0;
            for (uint32_t k = 0;  k < N; k++) { sum += A.c[i][k]*B.c[j][k]; }
            rn_test_tools_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
              "r3x3_mul_tr error"
            );
          }
      }
  }

void test_r3x3_tr_mul(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_tr_mul ---\n"); }
    r3x3_t A, B, C;
    r3x3_throw(&A, 0);
    r3x3_throw(&B, 0);
    r3x3_tr_mul(&A, &B, &C);
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < N; j++)
          { double sum = 0.0;
            for (uint32_t k = 0;  k < N; k++) { sum += A.c[k][i]*B.c[k][j]; }
            rn_test_tools_check_eps(C.c[i][j], sum, 0.000000001*fabs(sum), NO, NO,
              "r3x3_tr_mul error"
            );
          }
      }
  }

void test_r3x3_transp(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_transp ---\n"); }
    r3x3_t A, B;
    r3x3_throw(&A, 0);
    r3x3_transp(&A, &B);
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < N; j++)
          { rn_test_tools_check_eq(B.c[i][j],A.c[j][i], NO, NO, "r3x3_transp error (1)"); }
      }
    /* In-place transpose: */
    B = A;
    r3x3_transp(&B, &B);
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < N; j++)
          { rn_test_tools_check_eq(B.c[i][j],A.c[j][i], NO, NO, "r3x3_transp error (2)"); }
      }
  }

void test_r3x3_det(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_det ---\n"); }
    r3x3_t A;
    r3x3_throw(&A, 0);
    for (uint32_t i = 0;  i < N; i++)
      { uint32_t k = (i + 1) % N;
        for (uint32_t j = 0;  j < N; j++)
          { /* Check for linearity */
            double r = drandom();
            A.c[i][j] = r;
            double rr = r3x3_det(&A);

            double s = drandom();
            A.c[i][j] = s;
            double ss = r3x3_det(&A);

            double t = drandom();
            A.c[i][j] = r*(1-t) + s*t;
            double tt = r3x3_det(&A);
            double mag = fabs(rr) + fabs(ss) + fabs(tt);
            rn_test_tools_check_eps(tt, (rr*(1-t)+ss*t), 0.000000001*mag, NO, NO,
              "r3x3_det error (1)"
            );
          }

        /* Row swap test: */
        double r = r3x3_det(&A);
        for (uint32_t j = 0;  j < N; j++)
          { double t = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = t; }
        double rr = r3x3_det(&A);
        double rmag = fabs(r) + fabs(rr);
        rn_test_tools_check_eps(r, -rr, 0.000000001*rmag, NO, NO, "r3x3_det error (2)");

        /* Col swap test: */
        double s = r3x3_det(&A);
        for (uint32_t j = 0;  j < N; j++)
          { double t = A.c[j][i]; A.c[j][i] = A.c[j][k]; A.c[j][k] = t; }
        double ss = r3x3_det(&A);
        double smag = fabs(s) + fabs(ss);
        rn_test_tools_check_eps(s, -ss, 0.000000001*smag, NO, NO, "r3x3_det error (3)");
      }
  }

void test_r3x3_inv(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_inv ---\n"); }
    r3x3_t A, B, C;
    r3x3_throw(&A, 0);
    r3x3_inv(&A, &B);
    r3x3_mul(&A, &B, &C);
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < N; j++)
          { double val = (i == j ? 1.0 : 0.0);
            affirm(fabs(C.c[i][j] - val) < 0.000000001, "r3x3_inv error");
          }
      }
  }

void test_r3x3_norm__r3x3_norm_sqr__r3x3_mod_norm_sqr(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_norm, r3x3_norm_sqr, r3x3_mod_norm_sqr ---\n"); }
    r3x3_t A;
    r3x3_throw(&A, 0);
    double s = r3x3_norm_sqr(&A);
    double r = r3x3_norm(&A);
    double t = r3x3_mod_norm_sqr(&A);
    double ss = 0, tt = 0;
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < N; j++)
          { double Aij = A.c[i][j];
            double Dij = (i == j ? Aij - 1 : Aij);
            ss += Aij*Aij;
            tt += Dij*Dij;
          }
      }
    affirm(ss >= 0, "r3x3_norm_sqr error");
    affirm(fabs(ss - s) < 0.000000001, "r3x3_norm_sqr error");
    double rr = sqrt(ss);
    affirm(fabs(rr - r) < 0.000000001, "r3x3_norm error");
    affirm(tt >= 0, "r3x3_mod_norm_sqr error");
    affirm(fabs(tt - t) < 0.000000001, "r3x3_mod_norm_sqr error");
  }

void test_r3x3_normalize(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_normalize ---\n"); }
    r3x3_t A, B;
    r3x3_throw(&A, 0);
    B = A;
    double s = r3x3_norm(&B);
    double ss = r3x3_normalize(&B);
    affirm(fabs(ss - s) < 0.000000001, "r3x3_normalize result error");
    double t = r3x3_norm(&B);
    double tt = 1.0;
    affirm(fabs(tt - t) < 0.000000001, "r3x3_normalize norm error");
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < N; j++)
          { double Aij = A.c[i][j];
            double Bij = B.c[i][j];
            affirm(fabs(Bij*ss - Aij) < 0.000000001, "r3x3_normalize elem error");
          }
      }
  }

void test_r3x3_diff_sqr(bool_t verbose)
  {
    bool_t debug = FALSE;
    if (verbose) { fprintf(stderr, "--- r3x3_diff_sqr ---\n"); }
    
    r3x3_t A, B, R;
    r3x3_throw(&A, 0);
    r3x3_throw(&B, 0);
    r3x3_throw(&R, 0);
    uint32_t iz = uint32_abrandom(0, 2);
    uint32_t jz = uint32_abrandom(0, 2);
    R.c[iz][jz] = 0;
    if (debug)
      { tr3_print_matrix("A", &A);
        tr3_print_matrix("B", &B);
        tr3_print_matrix("R", &R);
      }
    double dabs2, drel2;
    r3x3_diff_sqr(&A, &B, &R, &dabs2, &drel2);
    double cabs2 = 0, crel2 = 0;
    for (uint32_t j = 0;  j < N; j++)
      { for (uint32_t i = 0;  i < N; i++)
          { double rij = R.c[i][j];
            if (rij != 0.0)
              { double d = A.c[i][j] - B.c[i][j];
                cabs2 += d*d;
                crel2 += (d/rij)*(d/rij);
              }
          }
      }
    tr3_check_num_eps("drel2", drel2, crel2, 0.0000001, "r3x3_diff_sqr failed");
    tr3_check_num_eps("dabs2", dabs2, cabs2, 0.0000001, "r3x3_diff_sqr failed");
  }

void test_r3x3_u_v_rotation(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_u_v_rotation ---\n"); }
    r3x3_t A;
    r3_t a, b, c;
    r3_throw_dir(&a);
    r3_throw_dir(&b);
    r3x3_u_v_rotation(&a, &b, &A);
    r3x3_map_row(&a, &A, &c);
    for (uint32_t i = 0;  i < N; i++)
      { affirm(fabs(b.c[i] - c.c[i]) < 0.000000001, "r3x3_u_v_rotation error"); }
  }

void test_r3x3_throw_rotation(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_throw_rotation ---\n"); }
    r3x3_t A;
    r3x3_throw_rotation(&A);
    rmxn_test_tools_check_all_different(N, N, &(A.c[0][0]), "r3x3_throw_totation probable error");
    /* A rotation's inverse is its transpose: */
    r3x3_t B;
    r3x3_mul_tr(&A, &A, &B);
    for (uint32_t i = 0;  i < N; i++)
      { for (uint32_t j = 0;  j < N; j++)
          { double Bij = B.c[i][j];
            double Bij_exp = (i == j ? 1.0 : 0.0);
            rn_test_tools_check_eps(Bij,Bij_exp, 1.0e-14, NO, NO, "r3x3_throw_rotation error (1)"); 
          }
      }
    /* The determinant must be {+1}: */
    double det = r3x3_det(&A);
    rn_test_tools_check_eps(det,1.0, 1.0e-14, NO, NO, "r3x3_throw_rotation error (2)"); 
  }

void test_r3x3_print(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r3x3_print ---\n"); }
    r3x3_t A;
    if (verbose)
      { r3x3_throw(&A, 0);
        fprintf(stderr, "A = ");
        r3x3_print(stderr, &A);
        fputc('\n', stderr);
      }
  }

/* TESTS OF {r3_hedron} FUNCTIONS */

void test_r3_hedron(bool_t verbose)
  {
    test_r3_hedron_tetra_vertices(verbose);
    test_r3_hedron_octa_vertices(verbose);
    test_r3_hedron_hexa_vertices(verbose);
    test_r3_hedron_icosa_vertices(verbose);
    test_r3_hedron_dodeca_vertices(verbose);
    
    if (verbose)
      { 
        fprintf(stderr, "!! r3_hedron_cylinder NOT TESTED\n");
      }
  }  

/* TESTS OF {r3_bezier} FUNCTIONS */

void test_r3_bezier(bool_t verbose)
  {
    if (verbose)
      { 
        fprintf(stderr, "!! r3_bezier_eval NOT TESTED\n");
        fprintf(stderr, "!! r3_bezier_length_estimate NOT TESTED\n");
        fprintf(stderr, "!! r3_bezier_split NOT TESTED\n");
      }
  }  

/* TESTS OF {r3_path} FUNCTIONS */

void test_r3_path(bool_t verbose)
  {
    if (verbose)
      { 
        fprintf(stderr, "!! r3_path_bezier_from_states NOT TESTED\n");
        fprintf(stderr, "!! r3_path_interpolate_some NOT TESTED\n");
        fprintf(stderr, "!! r3_path_length_estimate NOT TESTED\n");
        fprintf(stderr, "!! r3_path_state_debug NOT TESTED\n");
        fprintf(stderr, "!! r3_path_state_from_r3_motion_state NOT TESTED\n");
        fprintf(stderr, "!! r3_path_state_from_r3_motion_state NOT TESTED\n");
        fprintf(stderr, "!! r3_path_state_gen_print NOT TESTED\n");
        fprintf(stderr, "!! r3_path_state_print NOT TESTED\n");
        fprintf(stderr, "!! r3_path_state_to_r3_motion_state NOT TESTED\n");
      }
  }  

/* TESTS OF {r3_motion} FUNCTIONS */

void test_r3_motion(bool_t verbose)
  {
    if (verbose)
      { 
        fprintf(stderr, "!! r3_motion_circle NOT TESTED\n");
        fprintf(stderr, "!! r3_motion_helix NOT TESTED\n");
        fprintf(stderr, "!! r3_motion_sample_uniform NOT TESTED\n");
        fprintf(stderr, "!! r3_motion_state_canonical NOT TESTED\n");
        fprintf(stderr, "!! r3_motion_state_compose NOT TESTED\n");
        fprintf(stderr, "!! r3_motion_state_debug NOT TESTED\n");
        fprintf(stderr, "!! r3_motion_state_gen_print NOT TESTED\n");
        fprintf(stderr, "!! r3_motion_state_print NOT TESTED\n");      
        
      }
  }  

/* TESTING TOOLS */
  
void tr3_print_matrix(char *name, r3x3_t *A)
  { fprintf(stderr, "%s = ", name);
    r3x3_gen_print(stderr, A, "%+10.6f", NULL,NULL,NULL, NULL,NULL,NULL);
    fprintf(stderr, "\n");
  }

void tr3_check_regular_polyhedron(char *func, double R, double L, uint32_t n, r3_t p[], uint32_t deg)
  {
    /* Not a complete test... */
    R = fabs(R); /* We can't distinguish {R} from {-R}. */
    /* Find the smallest nonzero vertex-vertex distance {dmin}: */
    double dmin = +INF;
    for (uint32_t i = 0;  i < n; i++)
      { for (uint32_t j = 0;  j < i; j++)
          { dmin = fmin(dmin, r3_dist(&(p[i]), &(p[j]))); }
      }
    rn_test_tools_check_eps(dmin, L, 0.000001*R, NULL, NULL, "polyhedron has wrong side");
    /* Check each vertex: */
    for (uint32_t i = 0;  i < n; i++)
      { /* Check distance from origin: */
        double Ri = r3_norm(&(p[i]));
        rn_test_tools_check_eps(Ri, R, 0.000001*R, &i, NULL, "vertex has wrong radius");
        /* Check number of nearest neighbors: */
        uint32_t degi = 0;
        for (uint32_t j = 0;  j < n; j++)
          { double dij = r3_dist(&(p[i]), &(p[j]));
            if (fabs(dij - L) < 0.0001*R) { degi++; }
          }
        rn_test_tools_check_eq(degi, deg, &i, NULL, "vertex has wrong degree");
      }
  }

void tr3_check_num_eps(char *name, double x, double y, double eps, char *msg)
  { double diff = fabs(x - y);
    if (diff > eps)
      { fprintf(stderr, " ** %s: %+20.16e %+20.16e", name, x, y);
        fprintf(stderr, " diff = %+20.16e  max = %+20.16e - %s\n", diff, eps, msg);
        exit(1);
      }
  }
