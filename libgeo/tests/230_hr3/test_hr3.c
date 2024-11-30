/* Test program for hr3.h  */
/* Last edited on 2024-11-20 18:17:41 by stolfi */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <jsrandom.h>
#include <jsstring.h>
#include <affirm.h>

#include <rn.h>
#include <r3x3.h>
#include <r3.h>
#include <r2x2.h>
#include <r2.h>

#include <hr3.h>
#include <hr3_test_tools.h>

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_hr3(bool_t verbose);

void test_hr4_to_from_r3(bool_t verbose);
void test_hr3_pt_pt_diff(bool_t verbose);
void test_hr3_side(bool_t verbose);
void test_hr3_orient(bool_t verbose);
void test_hr3_plane_from_three_points(bool_t verbose);
void test_hr3_point_from_three_planes(bool_t verbose);
void test_hr3_point_point_dir(bool_t verbose);
void test_hr3_plane_normal(bool_t verbose);


int32_t main (int32_t argc, char **argv)
  {
    srand(1993);
    srandom(1933);

    for (uint32_t i = 0;  i < 100; i++) test_hr3(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_hr3(bool_t verbose)
  {
    if (verbose)
      { fprintf(stderr,
          "sizeof(hr3_point_t) = %lu  %d*sizeof(double) = %lu\n",
          sizeof(hr3_point_t), NH, NH*sizeof(double)
        );
      }

    test_hr4_to_from_r3(verbose);
    test_hr3_pt_pt_diff(verbose);
    test_hr3_side(verbose);
    test_hr3_orient(verbose);
    test_hr3_plane_from_three_points(verbose);
    test_hr3_point_from_three_planes(verbose);
    test_hr3_point_point_dir(verbose);
    test_hr3_plane_normal(verbose);
  }

void test_hr4_to_from_r3(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_from_r3 ---\n"); }
    r3_t pc; r3_throw_cube(&pc);
    hr3_point_t p = hr3_from_r3(&pc);
    affirm(p.c.c[0] == 1.0, "hr3_from_r3 error(1)");
    
    if (verbose) { fprintf(stderr, "--- r3_from_hr3 ---\n"); }
    r4_throw_cube(&(p.c));
    pc = r3_from_hr3(&p);
    hr3_point_t q = hr3_from_r3(&pc);
    { double tol = fabs(p.c.c[0])*1.0e-12;
      for (uint32_t i = 1;  i <= NC; i++)
        { double di = q.c.c[i]*p.c.c[0] - p.c.c[i];
          affirm(fabs(di) < tol, "r3_from_hr3 error(1)");
        }
    }
  }

void test_hr3_pt_pt_diff(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_pt_pt_diff ---\n"); }
    { /* Check zero distance: */
      hr3_point_t p; r4_throw_cube(&(p.c));
      double dpp = hr3_pt_pt_diff(&p, &p); 
      check_eq(dpp, 0.0, "hr3_pt_pt_diff(p,p) error(1)");
      
      /* Check symmetry: */
      r4_throw_cube(&(p.c));
      hr3_point_t q; r4_throw_cube(&(q.c));
      double dpq = hr3_pt_pt_diff(&p, &q); 
      double dqp = hr3_pt_pt_diff(&q, &p); 
      check_eq(dpq, dqp, "hr3_pt_pt_diff error(2)");
      
      /* Check range {[0 _ PI]}: */
      affirm(dpq >= 0.0, "hr3_pt_pt_diff error(sign)");
      affirm(dpq <= 1.000000001*M_PI, "hr3_pt_pt_diff error(max)");
      
      /* Generate two points {p,q} with known distance: */
      double alfa = 0.5*M_PI*drandom();
      double beta = 2.0*M_PI*drandom();
      p.c.c[0] = cos(alfa);
      p.c.c[1] = sin(alfa);
      p.c.c[2] = cos(beta);
      p.c.c[3] = sin(beta);
      
      q.c.c[0] = - p.c.c[1];
      q.c.c[1] = + p.c.c[0];
      q.c.c[2] = - p.c.c[3];
      q.c.c[3] = + p.c.c[2];
      r4_add(&(p.c), &(q.c), &(q.c));
      /* Check their distance: */
      double dex = M_PI/4; /* Expected. */
      double dob = hr3_pt_pt_diff(&p, &q); /* Actual. */
      hr3_test_check_eps(dob, dex, 1.0e-8, "hr3_pt_pt_diff error(3)");
      /* Check invariance under rotations: */
      for (uint32_t i = 0;  i < NH; i++)
        { uint32_t j = (i + 1) % NH; /* Another axis. */
          /* Rotate {p,q} by a random angle in {R^3} parallel to plane {i,j}: */
          double ang = 2*M_PI*drandom();
          double ca = cos(ang), sa = sin(ang);
          hr3_point_t pp = p, qq = q; /* Rotated points. */

          pp.c.c[i] = p.c.c[i]*ca + p.c.c[j]*sa;
          pp.c.c[j] = p.c.c[j]*ca - p.c.c[i]*sa;

          qq.c.c[i] = q.c.c[i]*ca + q.c.c[j]*sa;
          qq.c.c[j] = q.c.c[j]*ca - q.c.c[i]*sa;
          
          /* Check whether distance is preserved: */
          double drt = hr3_pt_pt_diff(&pp, &qq);
          hr3_test_check_eps(drt, dex, 1.0e-8, "hr3_pt_pt_diff error(4)");
        }
    }
  }

void test_hr3_side(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_side ---\n"); }
    hr3_point_t p; r4_throw_cube(&(p.c));
    hr3_plane_t L; r4_throw_cube(&(L.f));
    { double dd = r4_dot(&(p.c), &(L.f));
      sign_t sgn = hr3_side(&p, &L);
      affirm(((dd == 0.0) && (sgn == 0)) || (dd*sgn > 0), "hr3_side error(1)");
    }
  }

void test_hr3_orient(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_orient ---\n"); }
    hr3_point_t p; r4_throw_cube(&(p.c));
    hr3_point_t q; r4_throw_cube(&(q.c));
    hr3_point_t r; r4_throw_cube(&(r.c));
    hr3_point_t s; r4_throw_cube(&(s.c));
    { double dd = r4_det(&(p.c), &(q.c), &(r.c), &(s.c));
      sign_t sgn1 = hr3_orient(&p, &q, &r, &s);
      affirm(((dd == 0.0) && (sgn1 == 0)) || (dd*sgn1 > 0), "hr3_orient error(1)");
      sign_t sgn2 = hr3_orient(&p, &p, &q, &r);
      affirm(sgn2 == 0, "hr3_orient error(2)");
    }
  }

void test_hr3_plane_from_three_points(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_plane_from_three_points ---\n"); }
    hr3_point_t p; r4_throw_cube(&(p.c));
    hr3_point_t q; r4_throw_cube(&(q.c));
    hr3_point_t r; r4_throw_cube(&(r.c));
    hr3_point_t s; r4_throw_cube(&(s.c)); /* Random test point */
    hr3_plane_t L = hr3_plane_from_three_points(&p, &q, &r);
    
    double tp = r4_norm(&(L.f))*r4_norm(&(p.c))*1.0e-12;
    double dp = r4_dot(&(p.c), &(L.f));
    affirm(fabs(dp) < tp, "hr3_plane_from_three_points error(1)");
    
    double tq = r4_norm(&(L.f))*r4_norm(&(q.c))*1.0e-12;
    double dq = r4_dot(&(q.c), &(L.f));
    affirm(fabs(dq) < tq, "hr3_plane_from_three_points error(2)");
    
    sign_t sgnL = hr3_side(&s, &L);
    sign_t sgnO = hr3_orient(&p, &q, &r, &s);
    affirm(sgnL == sgnO, "hr3_plane_from_three_points error(3)");
  }
    
void test_hr3_point_from_three_planes(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_point_from_three_planes ---\n"); }
    hr3_plane_t J; r4_throw_cube(&(J.f)); 
    hr3_plane_t K; r4_throw_cube(&(K.f)); 
    hr3_plane_t L; r4_throw_cube(&(L.f)); 
    hr3_point_t pm = hr3_point_from_three_planes(&J, &K, &L);
    
    /* Consitency with {hr3_plane_from_three_points}: */
    hr3_point_t p; p.c = L.f;
    hr3_point_t q; q.c = K.f;
    hr3_point_t r; r.c = J.f;
    hr3_plane_t Lj = hr3_plane_from_three_points(&p, &q, &r);
    for (uint32_t i = 0;  i < NH; i++)
      { check_eq(pm.c.c[i], Lj.f.c[i], "hr3_point_from_three_planes error(1)"); }
  }
    
void test_hr3_point_point_dir(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_point_point_dir ---\n"); }
    hr3_point_t p; r4_throw_cube(&(p.c)); p.c.c[0] = fabs(p.c.c[0]);
    hr3_point_t q; r4_throw_cube(&(q.c)); q.c.c[0] = fabs(q.c.c[0]);
    { r3_t pc = r3_from_hr3(&p);
      r3_t qc = r3_from_hr3(&q);
      r3_t upq = hr3_point_point_dir(&p, &q);
      r3_t vpq;
      r3_sub(&qc, &pc, &vpq);
      r3_dir(&vpq, &vpq);
      for (uint32_t i = 0;  i < NC; i++)
        { hr3_test_check_eps(upq.c[i], vpq.c[i], 1.0e-12, "hr3_point_point_dir error"); }
    }
  }
    
void test_hr3_plane_normal(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_plane_normal ---\n"); }
    hr3_plane_t L; r4_throw_cube(&(L.f));
    { r3_t nL = hr3_plane_normal(&L);
      r3_t mL = (r3_t){{ L.f.c[1], L.f.c[2], L.f.c[3] }};
      double mLmag = r3_dir(&mL, &mL);
      assert(mLmag != 0);
      double tol = 1.0e-12;
      for (uint32_t i = 0;  i < NC; i++)
        { hr3_test_check_eps(nL.c[i], mL.c[i], tol, "hr3_plane_normal error"); }
    }
  }

