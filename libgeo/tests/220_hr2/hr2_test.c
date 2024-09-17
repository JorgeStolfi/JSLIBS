/* hr2test --- test program for hr2.h  */
/* Last edited on 2024-09-17 17:28:38 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <flt.h>
#include <jsrandom.h>
#include <jsstring.h>
#include <affirm.h>

#include <rn.h>
#include <r3x3.h>
#include <r3.h>
#include <r3_extra.h>
#include <r2x2.h>
#include <r2.h>
#include <hr2.h>

#include <hr2_test_tools.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_hr2(bool_t verbose);

void test_hr2_from_r2__r2_from_hr2(bool_t verbose);
void test_hr2_pt_pt_diff(bool_t verbose);
void test_hr2_side(bool_t verbose);
void test_hr2_orient(bool_t verbose);
void test_hr2_join(bool_t verbose);
void test_hr2_meet(bool_t verbose);
void test_hr2_point_point_dir(bool_t verbose);
void test_hr2_line_dir(bool_t verbose);
void test_hr2_line_normal(bool_t verbose);
void test_hr2_L_inf_normalize_line__hr2_L_inf_normalize_point(bool_t verbose);
void test_hr2_dist__hr2_dist_sqr(bool_t verbose);

int32_t main (int32_t argc, char **argv)
  {
    int32_t i;
    srand(1993);
    srandom(1933);

    for (i = 0; i < 100; i++) test_hr2(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_hr2(bool_t verbose)
  {
    if (verbose)
      { fprintf(stderr,
          "sizeof(hr2_point_t) = %lu  %d*sizeof(double) = %lu\n",
          sizeof(hr2_point_t), NH, NH*sizeof(double)
        );
      }

    test_hr2_from_r2__r2_from_hr2(verbose);
    test_hr2_pt_pt_diff(verbose);
    test_hr2_side(verbose);
    test_hr2_orient(verbose);
    test_hr2_join(verbose);
    test_hr2_meet(verbose);
    test_hr2_point_point_dir(verbose);
    test_hr2_line_dir(verbose);
    test_hr2_line_normal(verbose);
    test_hr2_L_inf_normalize_line__hr2_L_inf_normalize_point(verbose);
    test_hr2_dist__hr2_dist_sqr(verbose);
  }



void test_hr2_from_r2__r2_from_hr2(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_from_r2 ---\n"); }
    r2_t pc; r2_throw_cube(&pc);
    hr2_point_t p = hr2_from_r2(&pc);
    affirm(p.c.c[0] == 1.0, "hr2_from_r2 error(1)");

    if (verbose) { fprintf(stderr, "--- r2_from_hr2 ---\n"); }
    r3_throw_cube(&(p.c));
    pc = r2_from_hr2(&p);
    hr2_point_t q = hr2_from_r2(&pc);
    { double tol = fabs(p.c.c[0])*1.0e-12;
      for (int32_t i = 1; i <= NC; i++)
        { double di = q.c.c[i]*p.c.c[0] - p.c.c[i];
          affirm(fabs(di) < tol, "r2_from_hr2 error(1)");
        }
    }
  }

void test_hr2_pt_pt_diff(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pt_pt_diff ---\n"); }
    { /* Check zero distance: */
      hr2_point_t p; r3_throw_cube(&(p.c));
      double dpp = hr2_pt_pt_diff(&p, &p);
      check_eq(dpp, 0.0, "hr2_pt_pt_diff(p,p) error(1)");

      /* Check symmetry: */
      r3_throw_cube(&(p.c));
      hr2_point_t q; r3_throw_cube(&(q.c));
      double dpq = hr2_pt_pt_diff(&p, &q);
      double dqp = hr2_pt_pt_diff(&q, &p);
      check_eq(dpq, dqp, "hr2_pt_pt_diff error(2)");

      /* Check range {[0 _ PI]}: */
      affirm(dpq >= 0.0, "hr2_pt_pt_diff error(sign)");
      affirm(dpq <= 1.000000001*M_PI, "hr2_pt_pt_diff error(max)");

      /* Generate two points {p,q} with known distance: */
      double alfa = 0.5*M_PI*drandom();
      double beta = 2.0*M_PI*drandom();
      p.c.c[0] = sin(alfa)*cos(beta);
      p.c.c[1] = sin(alfa)*sin(beta);
      p.c.c[2] = cos(alfa);

      q.c.c[0] = - p.c.c[0];
      q.c.c[1] = - p.c.c[1];
      q.c.c[2] = + p.c.c[2];
      /* Check their distance: */
      double dex = 2*alfa; /* Expected. */
      if (dex > M_PI) { dex = 2*M_PI - dex; }
      double dob = hr2_pt_pt_diff(&p, &q); /* Actual. */
      hr2_test_check_eps(dob, dex, 1.0e-8, "hr2_pt_pt_diff error(3)");
      /* Check invariance under rotations: */
      for (int32_t i = 0; i < NH; i++)
        { int32_t j = (i + 1) % NH; /* Another axis. */
          /* Rotate {p,q} by a random angle in {R^3} parallel to plane {i,j}: */
          double ang = 2*M_PI*drandom();
          double ca = cos(ang), sa = sin(ang);
          hr2_point_t pp = p, qq = q; /* Rotated points. */

          pp.c.c[i] = p.c.c[i]*ca + p.c.c[j]*sa;
          pp.c.c[j] = p.c.c[j]*ca - p.c.c[i]*sa;

          qq.c.c[i] = q.c.c[i]*ca + q.c.c[j]*sa;
          qq.c.c[j] = q.c.c[j]*ca - q.c.c[i]*sa;

          /* Check whether distance is preserved: */
          double drt = hr2_pt_pt_diff(&pp, &qq);
          hr2_test_check_eps(drt, dex, 1.0e-8, "hr2_pt_pt_diff error(4)");
        }
    }
  }

void test_hr2_side(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_side ---\n"); }
    hr2_point_t p; r3_throw_cube(&(p.c));
    hr2_line_t L; r3_throw_cube(&(L.f));
    { double dd = r3_dot(&(p.c), &(L.f));
      sign_t sgn = hr2_side(&p, &L);
      affirm(((dd == 0.0) && (sgn == 0)) || (dd*sgn > 0), "hr2_side error(1)");
    }
  }

void test_hr2_orient(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_orient ---\n"); }
    hr2_point_t p; r3_throw_cube(&(p.c));
    hr2_point_t q; r3_throw_cube(&(q.c));
    hr2_point_t r; r3_throw_cube(&(r.c));
    { double dd = r3_det(&(p.c), &(q.c), &(r.c));
      sign_t sgn1 = hr2_orient(&p, &q, &r);
      affirm(((dd == 0.0) && (sgn1 == 0)) || (dd*sgn1 > 0), "hr2_orient error(1)");
      sign_t sgn2 = hr2_orient(&p, &p, &q);
      affirm(sgn2 == 0, "hr2_orient error(2)");
    }
  }

void test_hr2_join(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_join ---\n"); }
    hr2_point_t p; r3_throw_cube(&(p.c));
    hr2_point_t q; r3_throw_cube(&(q.c));
    hr2_point_t r; r3_throw_cube(&(r.c)); /* Random test point */
    hr2_line_t L = hr2_join(&p, &q);
    
    double tp = r3_norm(&(L.f))*r3_norm(&(p.c))*1.0e-12;
    double dp = r3_dot(&(p.c), &(L.f));
    affirm(fabs(dp) < tp, "hr2_join error(1)");

    double tq = r3_norm(&(L.f))*r3_norm(&(q.c))*1.0e-12;
    double dq = r3_dot(&(q.c), &(L.f));
    affirm(fabs(dq) < tq, "hr2_join error(2)");

    sign_t sgnL = hr2_side(&r, &L);
    sign_t sgnO = hr2_orient(&p, &q, &r);
    affirm(sgnL == sgnO, "hr2_join error(3)");
  }

void test_hr2_meet(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_meet ---\n"); }
    hr2_line_t L; r3_throw_cube(&(L.f));
    hr2_line_t M; r3_throw_cube(&(M.f));
    hr2_point_t r = hr2_meet(&L, &M);

    /* Consitency with {hr2_join}: */
    hr2_point_t p; p.c = L.f;
    hr2_point_t q; q.c = M.f;
    hr2_line_t N = hr2_join(&p, &q);
    for (int32_t i = 0; i < NH; i++)
      { check_eq(r.c.c[i], N.f.c[i], "hr2_meet error(1)"); }
  }

void test_hr2_point_point_dir(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_point_point_dir ---\n"); }
    hr2_point_t p; r3_throw_cube(&(p.c)); p.c.c[0] = fabs(p.c.c[0]);
    hr2_point_t q; r3_throw_cube(&(q.c)); q.c.c[0] = fabs(q.c.c[0]);
    { r2_t pc = r2_from_hr2(&p);
      r2_t qc = r2_from_hr2(&q);
      r2_t upq = hr2_point_point_dir(&p, &q);
      r2_t vpq;
      r2_sub(&qc, &pc, &vpq);
      r2_dir(&vpq, &vpq);
      for (int32_t i = 0; i < NC; i++)
        { hr2_test_check_eps(upq.c[i], vpq.c[i], 1.0e-12, "hr2_point_point_dir error"); }
    }
  }

void test_hr2_line_dir(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_line_dir ---\n"); }
    hr2_point_t p; r3_throw_cube(&(p.c)); p.c.c[0] = fabs(p.c.c[0]) + 0.00001;
    hr2_point_t q; r3_throw_cube(&(q.c)); q.c.c[0] = fabs(q.c.c[0]) + 0.00001;
    hr2_line_t L = hr2_join(&p, &q);
    { r2_t dL = hr2_line_dir(&L);
      r2_t eL = hr2_point_point_dir(&p, &q);;
      double tol = 1.0e-12;
      for (int32_t i = 0; i < NC; i++)
        { hr2_test_check_eps(dL.c[i], eL.c[i], tol, "hr2_line_dir error"); }
    }
  }

void test_hr2_line_normal(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_line_normal ---\n"); }
    hr2_line_t L; r3_throw_cube(&(L.f));
    { r2_t nL = hr2_line_normal(&L);
      r2_t mL = (r2_t){{ L.f.c[1], L.f.c[2] }};
      double mLmag = r2_dir(&mL, &mL);
      assert(mLmag != 0);
      double tol = 1.0e-12;
      for (int32_t i = 0; i < NC; i++)
        { hr2_test_check_eps(nL.c[i], mL.c[i], tol, "hr2_line_normal error"); }
    }
  }

void test_hr2_L_inf_normalize_line__hr2_L_inf_normalize_point(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "!! {hr2_L_inf_normalize_line} NOT TESTED\n"); }
    if (verbose) { fprintf(stderr, "!! {hr2_L_inf_normalize_point} NOT TESTED\n"); } 
  }

void test_hr2_dist__hr2_dist_sqr(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr2_dist, hr2_dist_sqr ---\n"); }
    double tol = 2.0e-15;
    r2_t p = (r2_t){{ 100*(drandom() - 0.5), 100*(drandom() - 0.5), }}; 
    r2_t q = (r2_t){{ 100*(drandom() - 0.5), 100*(drandom() - 0.5), }}; 
    double dist_sqr_exp = r2_dist_sqr(&p, &q);
    hr2_point_t hp = hr2_from_r2(&p);
    hr2_point_t hq = hr2_from_r2(&q);
    double dist_sqr_cmp = hr2_dist_sqr(&hp, &hq);
    demand(fabs(dist_sqr_cmp/dist_sqr_exp - 1.0) < tol, "{hr2_dist_sqr} failed");
    double dist_exp = r2_dist(&p, &q);
    double dist_cmp = hr2_dist(&hp, &hq);
    demand(fabs(dist_cmp/dist_exp - 1.0) < tol, "{hr2_dist} failed");
  }
