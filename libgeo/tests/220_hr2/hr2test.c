/* hr2test --- test program for hr2.h  */
/* Last edited on 2023-10-09 19:39:06 by stolfi */

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
#include <hr2test_tools.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_hr2(bool_t verbose);
void test_hr2_pmap(bool_t verbose);
void test_hr2_pmap_aff(bool_t verbose);

void test_hr2_to_from_r2(bool_t verbose);
void test_hr2_pt_pt_diff(bool_t verbose);
void test_hr2_side(bool_t verbose);
void test_hr2_orient(bool_t verbose);
void test_hr2_join(bool_t verbose);
void test_hr2_meet(bool_t verbose);
void test_hr2_point_point_dir(bool_t verbose);
void test_hr2_line_dir(bool_t verbose);
void test_hr2_line_normal(bool_t verbose);

void test_hr2_pmap_congruence_from_point_and_dir(bool_t verbose, bool_t flip);
void test_hr2_pmap_similarity_from_two_points(bool_t verbose, bool_t flip);
void test_hr2_pmap_from_four_points(bool_t verbose);
void test_hr2_pmap_gen_print(bool_t verbose);
void test_hr2_pmap_identity(bool_t verbose);
void test_hr2_pmap_inv(bool_t verbose);
void test_hr2_pmap_compose(bool_t verbose);
void test_hr2_pmap_inv_compose(bool_t verbose);
void test_hr2_pmap_is_affine(bool_t verbose);
void test_hr2_pmap_is_identity(bool_t verbose);
void test_hr2_pmap_line(bool_t inv, bool_t verbose);
void test_hr2_pmap_point(bool_t inv, bool_t verbose);
void test_hr2_pmap_print(bool_t verbose);
void test_hr2_pmap_r2_point(bool_t verbose);
void test_hr2_pmap_aff_from_mat_and_disp(bool_t verbose);
void test_hr2_pmap_aff_from_three_points(bool_t verbose);
void test_hr2_pmap_diff_sqr(bool_t verbose);
void test_hr2_pmap_mismatch_sqr(bool_t verbose);
void test_hr2_pmap_aff_discr_sqr(bool_t verbose);
void test_hr2_pmap_deform_sqr(bool_t verbose);
void test_hr2_pmap_rotation(bool_t verbose);
void test_hr2_pmap_rotation_and_scaling(bool_t verbose);
void test_hr2_pmap_scaling(bool_t verbose);
void test_hr2_pmap_translation(bool_t verbose);


int32_t main (int32_t argc, char **argv)
  {
    int32_t i;
    srand(1993);
    srandom(1933);

    for (i = 0; i < 100; i++) test_hr2(i < 3);
    for (i = 0; i < 100; i++) test_hr2_pmap(i < 3);
    for (i = 0; i < 100; i++) test_hr2_pmap_aff(i < 3);
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

    test_hr2_to_from_r2(verbose);
    test_hr2_pt_pt_diff(verbose);
    test_hr2_side(verbose);
    test_hr2_orient(verbose);
    test_hr2_join(verbose);
    test_hr2_meet(verbose);
    test_hr2_point_point_dir(verbose);
    test_hr2_line_dir(verbose);
    test_hr2_line_normal(verbose);
  }

void test_hr2_pmap(bool_t verbose)
  {

    /* Size: */
    if (verbose)
      { fprintf(stderr,
          "sizeof(hr2_pmap_t) = %lu  2*%d*%d*sizeof(double) = %lu\n",
          sizeof(hr2_pmap_t), NH, NH, 2*NH*NH*sizeof(double)
        );
      }

    test_hr2_pmap_point(FALSE, verbose);
    test_hr2_pmap_point(TRUE, verbose);
    test_hr2_pmap_line(FALSE, verbose);
    test_hr2_pmap_line(TRUE, verbose);
    test_hr2_pmap_r2_point(verbose);

    test_hr2_pmap_is_identity(verbose);
    test_hr2_pmap_identity(verbose);
    test_hr2_pmap_inv(verbose);
    test_hr2_pmap_compose(verbose);
    test_hr2_pmap_inv_compose(verbose);

    test_hr2_pmap_from_four_points(verbose);

    test_hr2_pmap_diff_sqr(verbose);
    test_hr2_pmap_deform_sqr(verbose);
    test_hr2_pmap_mismatch_sqr(verbose);

    test_hr2_pmap_print(verbose);
    test_hr2_pmap_gen_print(verbose);
  }

void test_hr2_pmap_aff(bool_t verbose)
  {
    test_hr2_pmap_is_affine(verbose);
    test_hr2_pmap_aff_discr_sqr(verbose);

    test_hr2_pmap_aff_from_mat_and_disp(verbose);

    test_hr2_pmap_translation(verbose);
    test_hr2_pmap_rotation(verbose);
    test_hr2_pmap_scaling(verbose);
    test_hr2_pmap_rotation_and_scaling(verbose);
    test_hr2_pmap_congruence_from_point_and_dir(verbose, FALSE);
    test_hr2_pmap_congruence_from_point_and_dir(verbose, TRUE);
    test_hr2_pmap_similarity_from_two_points(verbose, FALSE);
    test_hr2_pmap_similarity_from_two_points(verbose, TRUE);
    test_hr2_pmap_aff_from_three_points(verbose);

    /* TO BE COMPLETED !!! */

    if (verbose)
      {
        /* fprintf(stderr, "!! r2_aff_map_from_XXX NOT TESTED\n"); */
      }
  }

void test_hr2_pmap_diff_sqr(bool_t verbose)
  { if (verbose) { fprintf(stderr, "--- hr2_pmap_diff_sqr ---\n"); }
    hr2_pmap_t M, N;
    /* Throw two maps and normalizes their matrices: */
    for (int32_t k = 0; k < 2; k++)
      { N = M;
        h2tt_throw_pmap(&M);
        r3x3_normalize(&(M.dir));
        r3x3_normalize(&(M.inv));
      }
    /* Compute exepected diff sqr {d2_exp}: */
    double d2_exp = 0;
    for (int32_t i = 0; i < NH; i++)
      { for (int32_t j = 0; j < NH; j++)
          { double d = M.dir.c[i][j] - N.dir.c[i][j];
            d2_exp += d*d;
            d = M.inv.c[i][j] - N.inv.c[i][j];
            d2_exp += d*d;
          }
      }
    /* Scale the matrices by arbitrary amounts: */
    for (int32_t i = 0; i < NH; i++)
      { for (int32_t j = 0; j < NH; j++)
          { M.dir.c[i][j] *= 2;
            M.inv.c[i][j] *= 4;
            N.inv.c[i][j] *= 0.25;
            N.dir.c[i][j] *= 0.0625;
          }
      }
    /* Compute {hr2_pmap_diff_sqr} and compare: */
    double d2_cmp = hr2_pmap_diff_sqr(&M, &N);
    h2tt_check_eps(d2_cmp, d2_exp, 1.0e-13, "hr2_pmap_diff_sqr failed");
  }

void test_hr2_pmap_mismatch_sqr(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_mismatch_sqr ---\n"); }

    bool_t debug = FALSE;

    /* Choose a bunch of points: */
    int32_t np = 7;
    r2_t p1[np], p2[np];
    for (int32_t ip = 0; ip < np; ip++)
      { r2_throw_cube(&(p1[np]));
        r2_throw_cube(&(p2[np]));
      }
    /* Choose a projective map: */
    hr2_pmap_t M; h2tt_throw_pmap(&M);

    if (debug) { hr2_pmap_print(stderr, &M, "  M =\n", "\n"); }

    /* Compute the expected mismatch {m2_exp}: */
    double sum2 = 0.0;
    hr2_pmap_t N = hr2_pmap_inv(&M);
    for (int32_t ip = 0; ip < np; ip++)
      { r2_t *p1k = &(p1[ip]);
        r2_t *p2k = &(p2[ip]);
        r2_t q1k = hr2_pmap_r2_point(p1k, &M);
        r2_t q2k = hr2_pmap_r2_point(p2k, &N);
        double d2 = r2_dist_sqr(&q1k, &q2k);
        if (debug)
          { r2_gen_print(stderr, p1k, "%+8.5f", "  ( ", " ", " )");
            r2_gen_print(stderr, &q1k, "%+12.8f", " -> ( ", " ", " )");
            fprintf(stderr, " |%.6f| ", d2);
            r2_gen_print(stderr, &q2k, "%+12.8f", "( ", " ", " ) <- ");
            r2_gen_print(stderr, p2k, "%+8.5f", "( ", " ", " )\n");
          }
        sum2 += d2;
      }
    double m2_exp = sum2/np;

    /* Compute the mismatch by the library: */
    double m2_cmp = hr2_pmap_mismatch_sqr(&M, np, p1, p2);
    h2tt_check_eps(m2_cmp, m2_exp, 1.0e-13, "hr2_pmap_diff_sqr");
  }

void test_hr2_pmap_aff_discr_sqr(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_aff_discr_sqr ---\n"); }

    bool_t debug = FALSE;

    hr2_pmap_t M, N;
     h2tt_throw_aff_map(&M);
     h2tt_throw_aff_map(&N);
    if (debug)
      { h2tt_print_pmap("M", &M);
        h2tt_print_pmap("N", &N);
      }
    double mis2 = hr2_pmap_aff_discr_sqr(&M, &N);
    /* The integrand |(M - N)(u(t))|^2 is a sinusoidal function with freqs 0,1,2. */
    /* Can be integrated numerically with 5 or more samples. */
    int32_t nang = 7;
    double sum_d2 = 0;
    for (int32_t i = 0; i < nang; i++)
      { double ang = 2*M_PI*((double)i)/((double)nang);
        r2_t p = (r2_t){{ cos(ang), sin(ang) }};
        r2_t q = hr2_pmap_r2_point(&p, &M);
        r2_t r = hr2_pmap_r2_point(&p, &N);
        double d2 = r2_dist_sqr(&q, &r);
        sum_d2 += d2;
      }
    double cis2 = sum_d2/nang;
    h2tt_check_num_eps("mis", mis2, cis2, 0.0000001, "hr2_pmap_aff_discr_sqr failed");
  }

void test_hr2_pmap_deform_sqr(bool_t verbose)
  {
    fprintf(stderr, "!! {hr2_pmap_deform_sqr} NOT TESTED\n");
  }

void test_hr2_to_from_r2(bool_t verbose)
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
      h2tt_check_eps(dob, dex, 1.0e-8, "hr2_pt_pt_diff error(3)");
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
          h2tt_check_eps(drt, dex, 1.0e-8, "hr2_pt_pt_diff error(4)");
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
        { h2tt_check_eps(upq.c[i], vpq.c[i], 1.0e-12, "hr2_point_point_dir error"); }
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
        { h2tt_check_eps(dL.c[i], eL.c[i], tol, "hr2_line_dir error"); }
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
        { h2tt_check_eps(nL.c[i], mL.c[i], tol, "hr2_line_normal error"); }
    }
  }

void test_hr2_pmap_aff_from_mat_and_disp(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_aff_from_mat_disp ---\n"); }

    r2_t disp; r2_throw_cube(&disp);
    r2x2_t mat;
    for (int32_t i = 0; i < 2; i++)
      { r2_t s; r2_throw_cube(&s);
        mat.c[i][0] = s.c[0];
        mat.c[i][1] = s.c[1];
      }
    hr2_pmap_t A = hr2_pmap_aff_from_mat_and_disp(&mat, &disp);
    for (int32_t k = 0; k < 5; k++)
      { /* Should take the origin to {disp}* */
        r2_t o = (r2_t){{ 0.0, 0.0 }};
        h2tt_check_pmap_r2_point("o", &o, &A, FALSE, &disp, "hr2_pmap_aff_from_mat_disp failed");
        /* Test with a few other points: */
        for (int32_t kp = 0; kp < 3; kp++)
          { r2_t p; r2_throw_cube(&p);
            r2_t q; r2x2_map_row(&p, &mat, &q);
            r2_add(&disp, &q, &q);
            h2tt_check_pmap_r2_point("p", &p, &A, FALSE, &q, "hr2_pmap_aff_from_mat_disp failed");
          }
      }
  }

void test_hr2_pmap_scaling(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_scaling ---\n"); }

    r2_t scale;
    do { r2_throw_cube(&scale); } while ((scale.c[0] == 0) || (scale.c[1] == 0));

    hr2_pmap_t A = hr2_pmap_scaling(&scale);

    /* Check whether the map works: */
    r2_t oc = (r2_t){{ 0.0, 0.0 }};
    r2_t pc = (r2_t){{ 1.0, 0.0 }};
    r2_t qc = (r2_t){{ 0.0, 1.0 }};
    r2_t uc = (r2_t){{ 1.0, 1.0 }};

    r2_t pcM = (r2_t){{ scale.c[0], 0.0 }};
    r2_t qcM = (r2_t){{ 0.0, scale.c[1] }};
    r2_t ucM = (r2_t){{ scale.c[0], scale.c[1] }};

    h2tt_check_pmap_r2_point("o", &oc, &A, FALSE, &oc, "hr2_pmap_scaling failed");
    h2tt_check_pmap_r2_point("p", &pc, &A, FALSE, &pcM, "hr2_pmap_scaling failed");
    h2tt_check_pmap_r2_point("q", &qc, &A, FALSE, &qcM, "hr2_pmap_scaling failed");
    h2tt_check_pmap_r2_point("u", &uc, &A, FALSE, &ucM, "hr2_pmap_scaling failed");
  }

void test_hr2_pmap_congruence_from_point_and_dir(bool_t verbose, bool_t flip)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_congruence_from_point_and_dir ---\n"); }

    bool_t debug = FALSE;

    r2_t oc = (r2_t){{ 0.0, 0.0 }};
    r2_t pc = (r2_t){{ 1.0, 0.0 }};
    r2_t qc = (r2_t){{ 0.0, 1.0 }};
    r2_t rc = (r2_t){{ 1.0, 1.0 }};

    r2_t ocM; r2_throw_cube(&ocM);
    r2_t udM; r2_throw_dir(&udM);

    /* Determine the vector {v} that is going to be the image of vector {(0,1)}: */
    r2_t vdM = (r2_t){{ -udM.c[1], +udM.c[0] }};
    if (flip) { r2_neg(&vdM, &vdM); }

    if (debug)
      { fprintf(stderr, "  flip = %c\n", "FT"[flip]);
        r2_gen_print(stderr, &ocM, "%12.8f", "  ocM = ( ", " ", " )\n");
        r2_gen_print(stderr, &udM, "%12.8f", "  udM = ( ", " ", " )\n");
        r2_gen_print(stderr, &vdM, "%12.8f", "  vdM = ( ", " ", " )\n");
      }

    /* Compute the expected images of points {pc,qc,rc}: */
    r2_t pcM; r2_add(&ocM, &udM, &pcM);
    r2_t qcM; r2_add(&ocM, &vdM, &qcM);
    r2_t rcM; r2_add(&pcM, &vdM, &rcM);

    hr2_pmap_t A = hr2_pmap_congruence_from_point_and_dir(&ocM, &udM, flip);

    /* Check whether the map works: */
    h2tt_check_pmap_r2_point("o", &oc, &A, FALSE, &ocM, "hr2_pmap_congruence_from_point_and_dir failed");
    h2tt_check_pmap_r2_point("p", &pc, &A, FALSE, &pcM, "hr2_pmap_congruence_from_point_and_dir failed");
    h2tt_check_pmap_r2_point("q", &qc, &A, FALSE, &qcM, "hr2_pmap_congruence_from_point_and_dir failed");
    h2tt_check_pmap_r2_point("r", &rc, &A, FALSE, &rcM, "hr2_pmap_congruence_from_point_and_dir failed");
  }

void test_hr2_pmap_similarity_from_two_points(bool_t verbose, bool_t flip)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_similarity_from_two_points ---\n"); }

    bool_t debug = FALSE;

    r2_t oc = (r2_t){{ 0.0, 0.0 }};
    r2_t pc = (r2_t){{ 1.0, 0.0 }};
    r2_t qc = (r2_t){{ 0.0, 1.0 }};
    r2_t rc = (r2_t){{ 1.0, 1.0 }};

    r2_t ocM; r2_throw_cube(&ocM);
    r2_t pcM; r2_throw_cube(&pcM);

    /* Determine the vector {vdM} that is going to be the image of vector {(0,1)}: */
    r2_t udM; r2_sub(&pcM, &ocM, &udM);
    r2_t vdM = (r2_t){{ -udM.c[1], +udM.c[0] }};
    if (flip) { r2_neg(&vdM, &vdM); }

    if (debug)
      { fprintf(stderr, "  flip = %c\n", "FT"[flip]);
        r2_gen_print(stderr, &ocM, "%12.8f", " ocM = ( ", " ", " )\n");
        r2_gen_print(stderr, &pcM, "%12.8f", " pcM = ( ", " ", " )\n");
        r2_gen_print(stderr, &udM, "%12.8f", " udM =  ( ", " ", " )\n");
        r2_gen_print(stderr, &vdM, "%12.8f", " vdM =  ( ", " ", " )\n");
      }

    /* Compute the expected images of points {pc,qc,rc}: */
    r2_t qcM; r2_add(&ocM, &vdM, &qcM);
    r2_t rcM; r2_add(&pcM, &vdM, &rcM);

    hr2_pmap_t A = hr2_pmap_similarity_from_two_points(&ocM, &pcM, flip);

    /* Check whether the map works: */
    h2tt_check_pmap_r2_point("o", &oc, &A, FALSE, &ocM, "hr2_pmap_similarity_from_two_points failed");
    h2tt_check_pmap_r2_point("p", &pc, &A, FALSE, &pcM, "hr2_pmap_similarity_from_two_points failed");
    h2tt_check_pmap_r2_point("q", &qc, &A, FALSE, &qcM, "hr2_pmap_similarity_from_two_points failed");
    h2tt_check_pmap_r2_point("r", &rc, &A, FALSE, &rcM, "hr2_pmap_similarity_from_two_points failed");
  }

void test_hr2_pmap_from_four_points(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_from_four_points ---\n"); }

    bool_t debug = FALSE;

    auto void print_pair(hr2_point_t *p, hr2_point_t *pM);

    hr2_point_t p = (hr2_point_t){{{ 1.0, 0.0, 0.0 }}};
    hr2_point_t q = (hr2_point_t){{{ 0.0, 1.0, 0.0 }}};
    hr2_point_t r = (hr2_point_t){{{ 0.0, 0.0, 1.0 }}};
    hr2_point_t u;

    for (int32_t kt = 0; kt < 2*(1 << NH); kt++)
      { hr2_point_t pM_exp, qM_exp, rM_exp, uM_exp;
        u = (hr2_point_t){{{ 1.0, 1.0, 1.0 }}};
        if (kt < (1 << NH))
          { if (kt & 1) { u.c.c[0] = -u.c.c[0]; }
            if (kt & 2) { u.c.c[1] = -u.c.c[1]; }
            if (kt & 4) { u.c.c[2] = -u.c.c[2]; }
            pM_exp = p; qM_exp = q; rM_exp = r; uM_exp = u;
          }
        else
          { r3_throw_cube(&(pM_exp.c));
            r3_throw_cube(&(qM_exp.c));
            r3_throw_cube(&(rM_exp.c));
            r3_throw_cube(&(uM_exp.c));
          }

        hr2_pmap_t M = hr2_pmap_from_four_points(&pM_exp, &qM_exp, &rM_exp, &uM_exp);

        if (kt >= (1 << NH))
          { /* Choose {u = [±1,±1,±1]} based on inverse map of {uM_exp}: */
            u = hr2_pmap_inv_point(&uM_exp, &M);
            for (int32_t i = 0; i < NH; i++)
              { if (fabs(u.c.c[i]) < 1.0e-8)
                  { u.c.c[i] = 0.0; }
                else if (fabs(fabs(u.c.c[i]) - 1) < 1.0e-8)
                  { u.c.c[i] = (u.c.c[i] > 0 ? +1 : -1); }
                else
                  { fprintf(stderr, "** bad u[%d] = %+14.11f\n", i, u.c.c[i]); }
              }
          }

        if (debug)
          { fprintf(stderr, "  kt = %d goal:\n", kt);
            print_pair(&p, &pM_exp);
            print_pair(&q, &qM_exp);
            print_pair(&r, &rM_exp);
            print_pair(&u, &uM_exp);
            fprintf(stderr, "  actual:\n");
            hr2_point_t pM_cmp = hr2_pmap_point(&p, &M); print_pair(&p, &pM_cmp);
            hr2_point_t qM_cmp = hr2_pmap_point(&q, &M); print_pair(&q, &qM_cmp);
            hr2_point_t rM_cmp = hr2_pmap_point(&r, &M); print_pair(&r, &rM_cmp);
            hr2_point_t uM_cmp = hr2_pmap_point(&u, &M); print_pair(&u, &uM_cmp);
          }

        /* Check whether the map works: */
        h2tt_check_pmap_point("p", &p, FALSE, &M, FALSE, &pM_exp, "hr2_pmap_from_four_points failed");
        h2tt_check_pmap_point("q", &q, FALSE, &M, FALSE, &qM_exp, "hr2_pmap_from_four_points failed");
        h2tt_check_pmap_point("r", &r, FALSE, &M, FALSE, &rM_exp, "hr2_pmap_from_four_points failed");
        h2tt_check_pmap_point("u", &u, FALSE, &M, FALSE, &uM_exp, "hr2_pmap_from_four_points failed");
      }
    return;

    void print_pair(hr2_point_t *p, hr2_point_t *pM)
      { r3_gen_print(stderr, &(p->c), "%+4.1f", "    [ ", " ", " ]");
        r3_gen_print(stderr, &(pM->c), "%+14.10f", " -> [ ", " ", " ]\n");
      }
  }

void test_hr2_pmap_r2_point(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_r2_point, hr2_pmap_inv_r2_point,  ---\n"); }
    hr2_pmap_t M;  h2tt_throw_aff_map(&M);
    r2_t pc; r2_throw_cube(&pc);
    /* Compute {qc_exp} without using {hr2_pmap_r2_point}: */
    hr2_point_t ph = hr2_from_r2(&pc);
    hr2_point_t qh_exp = hr2_pmap_point(&ph, &M);
    r2_t qc_exp = r2_from_hr2(&qh_exp);
    h2tt_check_pmap_r2_point("p", &pc,     &M, FALSE, &qc_exp, "hr2_pmap_r2_point failed");
    h2tt_check_pmap_r2_point("q", &qc_exp, &M, TRUE,  &pc,     "hr2_pmap_inv_r2_point failed");
  }

void test_hr2_pmap_inv(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_inv ---\n"); }
    hr2_pmap_t M;  h2tt_throw_aff_map(&M);
    hr2_pmap_t N = hr2_pmap_inv(&M);
    for (int32_t k = 0; k < 5; k++)
      { hr2_point_t p = hr2_point_throw();
        hr2_point_t q = hr2_pmap_point(&p, &M);
        h2tt_check_pmap_point("q", &q, FALSE, &N, FALSE, &p, "hr2_pmap_inv failed");
      }
 }

void test_hr2_pmap_compose(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_compose ---\n"); }
    hr2_pmap_t M;  h2tt_throw_aff_map(&M);
    hr2_pmap_t N;  h2tt_throw_aff_map(&N);
    hr2_pmap_t P = hr2_pmap_compose(&M, &N);
    for (int32_t k = 0; k < 5; k++)
      { hr2_point_t p = hr2_point_throw();
        hr2_point_t q = hr2_pmap_point(&p, &M);
        hr2_point_t r = hr2_pmap_point(&q, &N);
        h2tt_check_pmap_point("p", &p, FALSE, &P, FALSE, &r, "hr2_pmap_compose failed");
      }
   }

void test_hr2_pmap_inv_compose(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_inv_compose ---\n"); }
    hr2_pmap_t M;  h2tt_throw_aff_map(&M);
    hr2_pmap_t Minv = hr2_pmap_inv(&M);
    hr2_pmap_t N;  h2tt_throw_aff_map(&N);
    hr2_pmap_t P = hr2_pmap_inv_compose(&M, &N);
    for (int32_t k = 0; k < 5; k++)
      { hr2_point_t p = hr2_point_throw();
        hr2_point_t q = hr2_pmap_point(&p, &Minv);
        hr2_point_t r = hr2_pmap_point(&q, &N);
        h2tt_check_pmap_point("p", &p, FALSE, &P, FALSE, &r, "hr2_pmap_inv_compose failed");
      }
   }

void test_hr2_pmap_translation(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_translation ---\n"); }
    r2_t r; r2_throw_cube(&r);
    hr2_pmap_t M = hr2_pmap_translation(&r);
    for (int32_t k = 0; k < 5; k++)
      { r2_t p; r2_throw_cube(&p);
        r2_t q; r2_add(&r, &p, &q);
        h2tt_check_pmap_r2_point("p", &p, &M, FALSE, &q, "hr2_pmap_translation failed");
      }
  }

void test_hr2_pmap_rotation(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_rotation ---\n"); }
    double ang = dabrandom(-1.58, +1.58);
    hr2_pmap_t M = hr2_pmap_rotation(ang);
    double ca = cos(ang);
    double sa = sin(ang);
    for (int32_t k = 0; k < 5; k++)
      { r2_t p; r2_throw_cube(&p);
        r2_t q;
        q.c[0] = + ca*p.c[0] - sa*p.c[1];
        q.c[1] = + sa*p.c[0] + ca*p.c[1];
        h2tt_check_pmap_r2_point("p", &p, &M, FALSE, &q, "hr2_pmap_rotation failed");
      }
  }

void test_hr2_pmap_rotation_and_scaling(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_rotation_and_scaling ---\n"); }
    double ang = dabrandom(-1.58, +1.58);
    double scale = dabrandom(0.5, 2.0);
    hr2_pmap_t M = hr2_pmap_rotation_and_scaling(ang, scale);
    double ca = cos(ang);
    double sa = sin(ang);
    for (int32_t k = 0; k < 5; k++)
      { r2_t p; r2_throw_cube(&p);
        r2_t q;
        q.c[0] = + ca*scale*p.c[0] - sa*scale*p.c[1];
        q.c[1] = + sa*scale*p.c[0] + ca*scale*p.c[1];
        h2tt_check_pmap_r2_point("p", &p, &M, FALSE, &q, "hr2_pmap_rotation_and_scaling failed");
      }
  }

void test_hr2_pmap_aff_from_three_points(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_aff_from_three_points ---\n"); }
    r2_t o; r2_throw_cube(&o);
    r2_t p; r2_throw_cube(&p);
    r2_t q; r2_throw_cube(&q);
    hr2_pmap_t M = hr2_pmap_aff_from_three_points(&o, &p, &q);
    /* Check whether the map works: */
    r2_t oo = (r2_t){{ 0.0, 0.0 }};
    r2_t pp = (r2_t){{ 1.0, 0.0 }};
    r2_t qq = (r2_t){{ 0.0, 1.0 }};

    h2tt_check_pmap_r2_point("o", &oo, &M, FALSE, &o, "hr2_pmap_aff_from_three_points failed");
    h2tt_check_pmap_r2_point("p", &pp, &M, FALSE, &p, "hr2_pmap_aff_from_three_points failed");
    h2tt_check_pmap_r2_point("q", &qq, &M, FALSE, &q, "hr2_pmap_aff_from_three_points failed");
  }

void test_hr2_pmap_point(bool_t inv, bool_t verbose)
  { char *tag = (inv ? "_inv" : "");
    if (verbose) { fprintf(stderr, "--- hr2_pmap%s_point ---\n", tag); }
    hr2_point_t p = hr2_point_throw();
    hr2_pmap_t M; h2tt_throw_pmap(&M);

    /* Compute expected images of {p} by {M.dir} and {M.inv}: */
    hr2_point_t q_exp;
    r3x3_map_row(&(p.c), (inv ? &(M.inv) : &(M.dir)), &(q_exp.c));

    /* Compare with library func results: */
    h2tt_check_pmap_point("p", &p, FALSE, &M, inv, &q_exp, txtcat3("hr2_pmap", tag,  "_point failed"));
  }

void test_hr2_pmap_line(bool_t inv, bool_t verbose)
  { 
    char *tag = (inv ? "_inv" : "");
    if (verbose) { fprintf(stderr, "--- hr2_pmap%s_line ---\n", tag); }
    hr2_line_t A = hr2_line_throw();
    hr2_pmap_t M; h2tt_throw_pmap(&M);

    /* Compute expected images of {p} by {M.dir} and {M.inv}: */
    hr2_line_t B_exp;
    r3x3_map_col((inv ? &(M.dir) : &(M.inv)), &(A.f), &(B_exp.f));

    /* Compare with library func results: */
    h2tt_check_pmap_line("A", &A, FALSE, &M, inv, &B_exp, txtcat3("hr2_pmap", tag, "_line failed"));

    /* Pick a point on the line: */
    hr2_point_t p;
    (void)r3_pick_ortho (&(A.f), &(p.c));
    hr2_point_t q = (inv ? hr2_pmap_inv_point(&p, &M) : hr2_pmap_point(&p, &M));
    double dot = r3_dot(&(B_exp.f), &(q.c));
    affirm(fabs(dot) < 1.0e-11, txtcat3(txtcat3("hr2_pmap", tag, "_line"), " incompat with ", txtcat3("hr2_pmap", tag, "_point")));
  }

void test_hr2_pmap_is_identity(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_is_identity} NOT TESTED\n"); }

void test_hr2_pmap_identity(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_identity} NOT TESTED\n"); }

void test_hr2_pmap_print(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_print} NOT TESTED\n"); }

void test_hr2_pmap_gen_print(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_gen_print} NOT TESTED\n"); }

void test_hr2_pmap_is_affine(bool_t verbose)
  { fprintf(stderr, "!! {test_hr2_pmap_is_affine} NOT TESTED\n"); }

