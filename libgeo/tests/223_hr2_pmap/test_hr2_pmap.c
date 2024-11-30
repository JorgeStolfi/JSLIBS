/* Test program for hr2_pmap.h  */
/* Last edited on 2024-11-24 08:58:10 by stolfi */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <jsrandom.h>
#include <jsstring.h>
#include <jsprintf.h>
#include <affirm.h>

#include <rn.h>
#include <rmxn.h>
#include <r3x3.h>
#include <r3.h>
#include <r2x2.h>
#include <r2.h>
#include <hr2.h>
#include <hr2_pmap.h>
#include <hr2_pmap_translation.h>
#include <hr2_pmap_congruence.h>
#include <hr2_pmap_similarity.h>
#include <hr2_pmap_affine.h>

#include <hr2_test_tools.h>
#include <hr2_pmap_test_tools.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */
  
#define Pr fprintf
#define Er stderr
  /* Shorter names. */

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_hr2_pmap(bool_t verbose);
void test_hr2_pmap_affine(bool_t verbose);

void test_hr2_pmap_sign__hr2_pmap_invert_sign__hr2_pmap_set_sign(bool_t verbose);

void test_hr2_pmap_inv(bool_t verbose);
void test_hr2_pmap_compose(bool_t verbose);
void test_hr2_pmap_inv_compose(bool_t verbose);
void test_hr2_pmap_line(bool_t inv, bool_t verbose);
void test_hr2_pmap_point(bool_t inv, bool_t verbose);
void test_hr2_pmap_r2_point__hr2_pmap_inv_r2_point(bool_t verbose);
void test_hr2_pmap_diff_sqr(bool_t verbose);
void test_hr2_pmap_mismatch_sqr(bool_t verbose);
void test_hr2_pmap_affine_discr_sqr(bool_t verbose);
void test_hr2_pmap_deform_sqr(bool_t verbose);

void test_hr2_pmap_identity(bool_t verbose);
void test_hr2_pmap_mirror(bool_t verbose);
void test_hr2_pmap_rotation(bool_t verbose);
void test_hr2_pmap_rotation_and_scaling(bool_t verbose);
void test_hr2_pmap_scaling(bool_t verbose);
void test_hr2_pmap_translation_from_disp(bool_t verbose);
void test_hr2_pmap_congruence_from_point_and_dir(bool_t verbose);
void test_hr2_pmap_similarity_from_two_points(bool_t verbose);
void test_hr2_pmap_affine_from_mat_and_disp(bool_t verbose);
void test_hr2_pmap_affine_from_three_points(bool_t verbose);
void test_hr2_pmap_from_parms(bool_t verbose);
void test_hr2_pmap_from_four_points(bool_t verbose);
void test_hr2_pmap_from_four_r2_points(bool_t verbose);
void test_hr2_pmap_r2_from_sign_class(bool_t verbose);

void test_hr2_pmap_is_valid(bool_t verbose);
void test_hr2_pmap_is_identity(bool_t verbose);
void test_hr2_pmap_is_translation(bool_t verbose);
void test_hr2_pmap_is_congruence(bool_t verbose);
void test_hr2_pmap_is_similarity(bool_t verbose);
void test_hr2_pmap_is_affine(bool_t verbose);

void test_hr2_pmap_print(bool_t verbose);
void test_hr2_pmap_gen_print(bool_t verbose);

void thr2_random_scale_and_print(char *Rname, hr2_pmap_t *R, uint32_t indent, bool_t verbose);
  /* Scales the matrices {R.dir} and {R.inv} by random homogeneous factors.
    if {verbose} is true, also prints {R} befoe and after scaling, with the given {name},
    indented by {indent} spaces. */
    
int32_t main (int32_t argc, char **argv)
  {
    srand(1993);
    srandom(1933);

    for (uint32_t i = 0;  i < 100; i++) test_hr2_pmap(i < 3);
    for (uint32_t i = 0;  i < 100; i++) test_hr2_pmap_affine(i < 3);
    fclose(Er);
    fclose(stdout);
    return (0);
  }


void test_hr2_pmap(bool_t verbose)
  {
    /* Size: */
    if (verbose)
      { Pr(Er,
          "sizeof(hr2_pmap_t) = %lu  2*%d*%d*sizeof(double) = %lu\n",
          sizeof(hr2_pmap_t), NH, NH, 2*NH*NH*sizeof(double)
        );
      }

    test_hr2_pmap_sign__hr2_pmap_invert_sign__hr2_pmap_set_sign(verbose);

    /* Testing in order of increasing strictness: */
    test_hr2_pmap_is_valid(verbose);
    test_hr2_pmap_is_affine(verbose);
    test_hr2_pmap_is_similarity(verbose);
    test_hr2_pmap_is_congruence(verbose);
    test_hr2_pmap_is_translation(verbose);
    test_hr2_pmap_is_identity(verbose);

    test_hr2_pmap_point(FALSE, verbose);
    test_hr2_pmap_point(TRUE, verbose);
    test_hr2_pmap_line(FALSE, verbose);
    test_hr2_pmap_line(TRUE, verbose);
    test_hr2_pmap_r2_point__hr2_pmap_inv_r2_point(verbose);

    test_hr2_pmap_identity(verbose);
    test_hr2_pmap_inv(verbose);
    test_hr2_pmap_compose(verbose);
    test_hr2_pmap_inv_compose(verbose);

    test_hr2_pmap_from_four_points(verbose);
    test_hr2_pmap_from_four_r2_points(verbose);
    test_hr2_pmap_r2_from_sign_class(verbose);

    test_hr2_pmap_diff_sqr(verbose);
    test_hr2_pmap_deform_sqr(verbose);
    test_hr2_pmap_mismatch_sqr(verbose);

    test_hr2_pmap_print(verbose);
    test_hr2_pmap_gen_print(verbose);

    if (verbose)
      {
      }
  }

void test_hr2_pmap_affine(bool_t verbose)
  {
    test_hr2_pmap_affine_discr_sqr(verbose);

    test_hr2_pmap_affine_from_mat_and_disp(verbose);
    test_hr2_pmap_affine_from_three_points(verbose);
    test_hr2_pmap_from_parms(verbose);

    test_hr2_pmap_identity(verbose);
    test_hr2_pmap_mirror(verbose);
    test_hr2_pmap_translation_from_disp(verbose);
    test_hr2_pmap_rotation(verbose);
    test_hr2_pmap_scaling(verbose);
    test_hr2_pmap_rotation_and_scaling(verbose);
    test_hr2_pmap_congruence_from_point_and_dir(verbose);
    test_hr2_pmap_congruence_from_point_and_dir(verbose);
    test_hr2_pmap_similarity_from_two_points(verbose);
    test_hr2_pmap_similarity_from_two_points(verbose);

    /* TO BE COMPLETED !!! */

    if (verbose)
      {
      }
  }

void test_hr2_pmap_diff_sqr(bool_t verbose)
  { if (verbose) { Pr(Er, "=== testing hr2_pmap_diff_sqr ==============================\n"); }
    hr2_pmap_t M, N;
    /* Throw two maps and normalize their matrices: */
    for (uint32_t k = 0;  k < 2; k++)
      { N = M;
        hr2_pmap_throw(&M);
        r3x3_normalize(&(M.dir));
        r3x3_normalize(&(M.inv));
        hr2_pmap_test_tools_check("M (normalized)", &M, 1.0e-8, "{r3x3_normalize} failed (1)");
      }
    /* Compute exepected diff sqr {d2_exp}: */
    double d2_exp = 0;
    for (uint32_t i = 0;  i < NH; i++)
      { for (uint32_t j = 0;  j < NH; j++)
          { double d = M.dir.c[i][j] - N.dir.c[i][j];
            d2_exp += d*d;
            d = M.inv.c[i][j] - N.inv.c[i][j];
            d2_exp += d*d;
          }
      }
    /* Scale the matrices by arbitrary amounts: */
    for (uint32_t i = 0;  i < NH; i++)
      { for (uint32_t j = 0;  j < NH; j++)
          { M.dir.c[i][j] *= 2;
            M.inv.c[i][j] *= 4;
            N.inv.c[i][j] *= 0.25;
            N.dir.c[i][j] *= 0.0625;
          }
      }
    /* Compute {hr2_pmap_diff_sqr} and compare: */
    double d2_cmp = hr2_pmap_diff_sqr(&M, &N);
    hr2_test_tools_check_eps(d2_cmp, d2_exp, 1.0e-13, "{hr2_pmap_diff_sqr} failed (1)");
  }

void test_hr2_pmap_mismatch_sqr(bool_t verbose)
  {
    bool_t debug = FALSE;
    if (debug) { verbose = TRUE; }
    if (verbose) { Pr(Er, "=== testing hr2_pmap_mismatch_sqr ==============================\n"); }

    /* Choose a bunch of points and weights: */
    uint32_t np = 7;
    r2_t p1[np], p2[np];
    double w[np];
    for (uint32_t ip = 0;  ip < np; ip++)
      { r2_throw_cube(&(p1[np]));
        r2_throw_cube(&(p2[np]));
        w[ip] = 0.000001 + 0.999999*drandom();
      }
    /* Choose a projective map: */
    hr2_pmap_t M; hr2_pmap_throw(&M);

    if (verbose) { hr2_pmap_test_tools_print(2, "M", &M); }

    /* Compute the expected mismatch {m2_exp}: */
    double sum_wD2 = 0.0;
    double sum_w = 0.0;
    hr2_pmap_t N = hr2_pmap_inv(&M);
    hr2_pmap_test_tools_check("N", &N, 1.0e-8, "{hr2_pmap_inv} failed (1)");
    for (uint32_t ip = 0;  ip < np; ip++)
      { r2_t *p1k = &(p1[ip]);
        r2_t *p2k = &(p2[ip]);
        r2_t q1k = hr2_pmap_r2_point(p1k, &M);
        r2_t q2k = hr2_pmap_r2_point(p2k, &N);
        r2_t q2k_bis = hr2_pmap_inv_r2_point(p2k, &M);
        demand(r2_dist(&q2k, &q2k_bis) < 1.0e-13, "bug in inverse mapping");
        double D2k_dir = r2_dist_sqr(&q1k, p2k);
        double D2k_inv = r2_dist_sqr(p1k, &q2k);
        double wk = w[ip];
        if (verbose)
          { r2_gen_print(Er, p1k, "%+8.5f", "  ( ", " ", " )");
            r2_gen_print(Er, &q1k, "%+12.8f", " -> ( ", " ", " )\n");
            r2_gen_print(Er, &q2k, "%+12.8f", "( ", " ", " )");
            r2_gen_print(Er, p2k, "%+8.5f", " <- ( ", " ", " )\n");
            Pr(Er, " D2_dir = %.6f  D2_inv = %.6f w = %.10f\n", D2k_dir, D2k_inv, wk);
          }
        sum_wD2 += wk*(D2k_dir + D2k_inv);
        sum_w += wk;
      }
    double m2_exp = 0.5*sum_wD2/sum_w;

    /* Compute the mismatch by the library: */
    double m2_cmp = hr2_pmap_mismatch_sqr(&M, np, p1, p2, w);
    hr2_test_tools_check_eps(m2_cmp, m2_exp, 1.0e-13, "hr2_pmap_diff_sqr");
  }

void test_hr2_pmap_affine_discr_sqr(bool_t verbose)
  {
    bool_t debug = FALSE;
    if(debug) { verbose = TRUE; }
    if (verbose) { Pr(Er, "=== testing hr2_pmap_affine_discr_sqr ==============================\n"); }


    hr2_pmap_t M, N;
    hr2_pmap_test_tools_throw_aff_map(&M); hr2_pmap_test_tools_throw_aff_map(&N);
    if (verbose)
      { hr2_pmap_test_tools_print(2, "M", &M); hr2_pmap_test_tools_print(2, "N", &N); }
    
    double mis2 = hr2_pmap_affine_discr_sqr(&M, &N);
    /* The integrand |(M - N)(u(t))|^2 is a sinusoidal function with freqs 0,1,2. */
    /* Can be integrated numerically with 5 or more samples. */
    uint32_t nang = 7;
    double sum_d2 = 0;
    for (uint32_t i = 0;  i < nang; i++)
      { double ang = 2*M_PI*((double)i)/((double)nang);
        r2_t p = (r2_t){{ cos(ang), sin(ang) }};
        r2_t q = hr2_pmap_r2_point(&p, &M);
        r2_t r = hr2_pmap_r2_point(&p, &N);
        double d2 = r2_dist_sqr(&q, &r);
        sum_d2 += d2;
      }
    double cis2 = sum_d2/nang;
    hr2_test_tools_check_num_eps("mis", mis2, cis2, 0.0000001, "{hr2_pmap_affine_discr_sqr} failed (1)");
  }

void test_hr2_pmap_deform_sqr(bool_t verbose)
  { if (verbose) { Pr(Er, "!! {hr2_pmap_deform_sqr} NOT TESTED\n"); } }

void test_hr2_pmap_affine_from_mat_and_disp(bool_t verbose)
  {
    if (verbose) { Pr(Er, "=== testing hr2_pmap_affine_from_mat_disp ==============================\n"); }

    r2_t disp; r2_throw_cube(&disp);
    r2x2_t mat;
    for (uint32_t i = 0;  i < 2; i++)
      { r2_t s; r2_throw_cube(&s);
        mat.c[i][0] = s.c[0];
        mat.c[i][1] = s.c[1];
      }
    hr2_pmap_t A = hr2_pmap_affine_from_mat_and_disp(&mat, &disp);
    hr2_pmap_test_tools_check("A", &A, 1.0e-8, "{hr2_pmap_affine_from_mat_and_disp} failed (1)"); 
    
    for (uint32_t k = 0;  k < 5; k++)
      { /* Should take the origin to {disp}* */
        r2_t o = (r2_t){{ 0.0, 0.0 }};
        hr2_pmap_test_tools_check_map_r2_point("o", &o, &A, FALSE, &disp, "{hr2_pmap_affine_from_mat_disp} failed (1)");
        /* Test with a few other points: */
        for (uint32_t kp = 0;  kp < 3; kp++)
          { r2_t p; r2_throw_cube(&p);
            r2_t q; r2x2_map_row(&p, &mat, &q);
            r2_add(&disp, &q, &q);
            hr2_pmap_test_tools_check_map_r2_point("p", &p, &A, FALSE, &q, "{hr2_pmap_affine_from_mat_disp} failed (2)");
          }
      }
  }

void test_hr2_pmap_from_parms(bool_t verbose)
  { if (verbose) { Pr(Er, "!! {hr2_pmap_from_parms} NOT TESTED\n"); } }

void test_hr2_pmap_scaling(bool_t verbose)
  {
    if (verbose) { Pr(Er, "=== testing hr2_pmap_scaling ==============================\n"); }

    r2_t scale;
    do { r2_throw_cube(&scale); } while ((scale.c[0] == 0) || (scale.c[1] == 0));

    hr2_pmap_t A = hr2_pmap_scaling(&scale);
    hr2_pmap_test_tools_check("A", &A, 1.0e-8, "{hr2_pmap_scaling} failed (1)");

    /* Check whether the map works: */
    r2_t oc = (r2_t){{ 0.0, 0.0 }};
    r2_t pc = (r2_t){{ 1.0, 0.0 }};
    r2_t qc = (r2_t){{ 0.0, 1.0 }};
    r2_t uc = (r2_t){{ 1.0, 1.0 }};

    r2_t pcM = (r2_t){{ scale.c[0], 0.0 }};
    r2_t qcM = (r2_t){{ 0.0, scale.c[1] }};
    r2_t ucM = (r2_t){{ scale.c[0], scale.c[1] }};

    hr2_pmap_test_tools_check_map_r2_point("o", &oc, &A, FALSE, &oc,  "{hr2_pmap_scaling} failed (1)");
    hr2_pmap_test_tools_check_map_r2_point("p", &pc, &A, FALSE, &pcM, "{hr2_pmap_scaling} failed (2)");
    hr2_pmap_test_tools_check_map_r2_point("q", &qc, &A, FALSE, &qcM, "{hr2_pmap_scaling} failed (3)");
    hr2_pmap_test_tools_check_map_r2_point("u", &uc, &A, FALSE, &ucM, "{hr2_pmap_scaling} failed (4)");
  }

void test_hr2_pmap_congruence_from_point_and_dir(bool_t verbose)
  {
    bool_t debug = FALSE;
    if (debug) { verbose = TRUE; }
    if (verbose) { Pr(Er, "=== testing hr2_pmap_congruence_from_point_and_dir ==============================\n"); }

    r2_t oc = (r2_t){{ 0.0, 0.0 }};
    r2_t pc = (r2_t){{ 1.0, 0.0 }};
    r2_t qc = (r2_t){{ 0.0, 1.0 }};
    r2_t rc = (r2_t){{ 1.0, 1.0 }};

    for (sign_t sgn = -1; sgn <= +1; sgn += 2)
      { r2_t ocM; r2_throw_cube(&ocM);
        r2_t udM; r2_throw_dir(&udM);
        /* Determine the vector {v} that is going to be the image of vector {(0,1)}: */
        r2_t vdM = (r2_t){{ -udM.c[1], +udM.c[0] }};
        if (sgn < 0) { r2_neg(&vdM, &vdM); }

        if (verbose)
          { Pr(Er, "  sgn = %+d\n", sgn);
            r2_gen_print(Er, &ocM, "%12.8f", "  ocM = ( ", " ", " )\n");
            r2_gen_print(Er, &udM, "%12.8f", "  udM = ( ", " ", " )\n");
            r2_gen_print(Er, &vdM, "%12.8f", "  vdM = ( ", " ", " )\n");
          }

        /* Compute the expected images of points {pc,qc,rc}: */
        r2_t pcM; r2_add(&ocM, &udM, &pcM);
        r2_t qcM; r2_add(&ocM, &vdM, &qcM);
        r2_t rcM; r2_add(&pcM, &vdM, &rcM);

        hr2_pmap_t A = hr2_pmap_congruence_from_point_and_dir(&ocM, &udM, sgn);
        hr2_pmap_test_tools_check("A", &A, 1.0e-8, "{hr2_pmap_congruence_from_point_and_dir} failed (1)");

        /* Check whether the map works: */
        hr2_pmap_test_tools_check_map_r2_point("o", &oc, &A, FALSE, &ocM, "{hr2_pmap_congruence_from_point_and_dir} failed (2)");
        hr2_pmap_test_tools_check_map_r2_point("p", &pc, &A, FALSE, &pcM, "{hr2_pmap_congruence_from_point_and_dir} failed (3)");
        hr2_pmap_test_tools_check_map_r2_point("q", &qc, &A, FALSE, &qcM, "{hr2_pmap_congruence_from_point_and_dir} failed (4)");
        hr2_pmap_test_tools_check_map_r2_point("r", &rc, &A, FALSE, &rcM, "{hr2_pmap_congruence_from_point_and_dir} failed (5)");
     }
  }

void test_hr2_pmap_similarity_from_two_points(bool_t verbose)
  {
    bool_t debug = FALSE;
    if (debug) { verbose = TRUE; }
    if (verbose) { Pr(Er, "=== testing hr2_pmap_similarity_from_two_points ==============================\n"); }


    r2_t oc = (r2_t){{ 0.0, 0.0 }};
    r2_t pc = (r2_t){{ 1.0, 0.0 }};
    r2_t qc = (r2_t){{ 0.0, 1.0 }};
    r2_t rc = (r2_t){{ 1.0, 1.0 }};

    for (sign_t sgn = -1; sgn <= +1; sgn += 2)
      { r2_t ocM; r2_throw_cube(&ocM);
        r2_t pcM; r2_throw_cube(&pcM);

        /* Determine the vector {vdM} that is going to be the image of vector {(0,1)}: */
        r2_t udM; r2_sub(&pcM, &ocM, &udM);
        r2_t vdM = (r2_t){{ -udM.c[1], +udM.c[0] }};
        if (sgn < 0) { r2_neg(&vdM, &vdM); }

        if (verbose)
          { Pr(Er, "  sgn = %+d\n", sgn);
            r2_gen_print(Er, &ocM, "%12.8f", " ocM = ( ", " ", " )\n");
            r2_gen_print(Er, &pcM, "%12.8f", " pcM = ( ", " ", " )\n");
            r2_gen_print(Er, &udM, "%12.8f", " udM =  ( ", " ", " )\n");
            r2_gen_print(Er, &vdM, "%12.8f", " vdM =  ( ", " ", " )\n");
          }

        /* Compute the expected images of points {pc,qc,rc}: */
        r2_t qcM; r2_add(&ocM, &vdM, &qcM);
        r2_t rcM; r2_add(&pcM, &vdM, &rcM);

        hr2_pmap_t A = hr2_pmap_similarity_from_two_points(&ocM, &pcM, sgn);
        hr2_pmap_test_tools_check("A", &A, 1.0e-8, "{hr2_pmap_similarity_from_two_points} failed (1)");

        /* Check whether the map works: */
        hr2_pmap_test_tools_check_map_r2_point("o", &oc, &A, FALSE, &ocM, "{hr2_pmap_similarity_from_two_points} failed (2)");
        hr2_pmap_test_tools_check_map_r2_point("p", &pc, &A, FALSE, &pcM, "{hr2_pmap_similarity_from_two_points} failed (3)");
        hr2_pmap_test_tools_check_map_r2_point("q", &qc, &A, FALSE, &qcM, "{hr2_pmap_similarity_from_two_points} failed (4)");
        hr2_pmap_test_tools_check_map_r2_point("r", &rc, &A, FALSE, &rcM, "{hr2_pmap_similarity_from_two_points} failed (5)");
      }
  }

void test_hr2_pmap_from_four_points(bool_t verbose)
  {
    bool_t debug = FALSE;
    if (debug) { verbose = TRUE; }
    if (verbose) { Pr(Er, "=== testing hr2_pmap_from_four_points ==============================\n"); }


    auto void print_pair(hr2_point_t *p, hr2_point_t *pM);

    /* Source points: */
    hr2_point_t p[4];
    p[0] = (hr2_point_t){{{ 1.0, 0.0, 0.0 }}};
    p[1] = (hr2_point_t){{{ 0.0, 1.0, 0.0 }}};
    p[2] = (hr2_point_t){{{ 0.0, 0.0, 1.0 }}};
    /* Will set {p[3]} later. */

    for (uint32_t kt = 0;  kt < 2*(1 << NH); kt++)
      { if (verbose) { Pr(Er, "  --- trial %d ---\n", kt); }
        bool_t ident = (kt < (1 << NH));  /* Try to get the identity map: */
        hr2_point_t pM_exp[4];
        if (ident)
          { /* Try all choices of {p[3] = [±1,±1,±1]} */
            /* and set {pM_exp[0..4]} so that the map is the identity: */
            hr2_point_t u = (hr2_point_t){{{ 1.0, 1.0, 1.0 }}};
            if (kt & 1) { u.c.c[0] = -u.c.c[0]; }
            if (kt & 2) { u.c.c[1] = -u.c.c[1]; }
            if (kt & 4) { u.c.c[2] = -u.c.c[2]; }
            p[3] = u;
            for (uint32_t kp = 0;  kp < 4; kp++) { pM_exp[kp] = p[kp]; }
          }
        else
          { /* Fix {p[3] = [1,1,1]} and throw {pM_exp} random: */
            p[3] = (hr2_point_t){{{ 1.0, 1.0, 1.0 }}};
            for (uint32_t kp = 0;  kp < 4; kp++) { r3_throw_cube(&(pM_exp[kp].c)); }
          }

        hr2_pmap_t M = hr2_pmap_from_four_points(&(pM_exp[0]), &(pM_exp[1]), &(pM_exp[2]), &(pM_exp[3]));
        hr2_pmap_test_tools_check("M", &M, 1.0e-8, "{hr2_pmap_from_four_points} failed (1)");
        
        if (verbose)
          { hr2_pmap_gen_print (Er, &M, "%+10.5f", "  map:\n", "    ", "  ", "\n", "[ ", " ", " ]", "\n"); }

        if (! ident)
          { /* Set {p[3] = [±1,±1,±1]} based on inverse map of {pM_exp[3]}: */
            hr2_point_t u = hr2_pmap_inv_point(&(pM_exp[3]), &M);
            for (uint32_t i = 0;  i < NH; i++)
              { if (fabs(fabs(u.c.c[i]) - 1) < 1.0e-8)
                  { u.c.c[i] = (u.c.c[i] > 0 ? +1 : -1); }
                else
                  { Pr(Er, "** bad u[%d] = %+14.11f\n", i, u.c.c[i]); }
              }
            p[3] = u;
          }

        if (verbose)
          { Pr(Er, "  trial kt = %d goal:\n", kt);
            for (uint32_t kp = 0;  kp < 4; kp++) 
              { print_pair(&(p[kp]), &(pM_exp[kp])); }
            Pr(Er, "  actual:\n");
            for (uint32_t kp = 0;  kp < 4; kp++) 
              { hr2_point_t pM_cmp = hr2_pmap_point(&(p[kp]), &M); 
                print_pair(&(p[kp]), &pM_cmp);
              }
            Pr(Er, "\n");
          }

        /* Check whether the map works: */
        for (uint32_t kp = 0;  kp < 4; kp++) 
          { char *lab = (char*[4]){ "p[0]", "p[1]", "p[2]", "p[3]" }[kp];
            hr2_pmap_test_tools_check_map_point(lab, &(p[kp]), FALSE, FALSE, &M, FALSE, &(pM_exp[kp]), "{hr2_pmap_from_four_points} failed (2)");
          }
      }
    return;

    void print_pair(hr2_point_t *p, hr2_point_t *pM)
      { r3_gen_print(Er, &(p->c), "%+6.3f", "    [ ", " ", " ]");
        r3_gen_print(Er, &(pM->c), "%+14.10f", " -> [ ", " ", " ]\n");
      }
  }

void test_hr2_pmap_from_four_r2_points(bool_t verbose)
  {
    bool_t debug = FALSE;
    if (debug) { verbose = TRUE; }
    if (verbose) { Pr(Er, "=== testing hr2_pmap_from_four_r2_points ==============================\n"); }

    auto void print_pair(r2_t *p, r2_t *pM);

    r2_t p[4];
    p[0] = (r2_t){{ 0.0, 0.0 }};
    p[1] = (r2_t){{ 1.0, 0.0 }};
    p[2] = (r2_t){{ 0.0, 1.0 }};
    /* Will set {p[3]} later. */

    for (uint32_t kt = 0;  kt < 8; kt++)
      { if (verbose) { Pr(Er, "  --- trial %d ---\n", kt); }
            
        bool_t ident = (kt < 4);  /* Try to get the identity map: */
        r2_t pM_exp[4];
        if (ident)
          { /* Try all choices {p[3]} {(+1,+1)}, {(-1,+1)}, {(+1,-1)}, or {(1/3,1/3)}, */
            /* and set {pM_exp[0..3] = p[0..3]} so that the map is the identity: */
            if (kt == 0)
              { p[3] = (r2_t){{ +1.0, +1.0 }}; }
            else if (kt == 1)
              { p[3] = (r2_t){{ -1.0, +1.0 }}; }
            else if (kt == 2)
              { p[3] = (r2_t){{ +1.0, -1.0 }}; }
            else if (kt == 3)
              { p[3] = (r2_t){{ +1/3.0, +1/3.0 }}; }
            for (uint32_t kp = 0;  kp < 4; kp++) { pM_exp[kp] = p[kp]; }
          }
        else
          { /* Fix {p[3] = (1,1)} and throw {pM_exp} random: */
            p[3] = (r2_t){{ 1.0, 1.0 }};
            for (uint32_t kp = 0;  kp < 4; kp++) { r2_throw_cube(&(pM_exp[kp])); }
          }

        hr2_pmap_t M = hr2_pmap_from_four_r2_points(&(pM_exp[0]), &(pM_exp[1]), &(pM_exp[2]), &(pM_exp[3]));
        hr2_pmap_test_tools_check("M", &M, 1.0e-8, "{hr2_pmap_from_four_r2_points} failed (1)");
        
        if (verbose)
          { hr2_pmap_gen_print (Er, &M, "%+10.5f", "  map:\n", "    ", "  ", "\n", "[ ", " ", " ]", "\n"); }

        if (! ident)
          { /* Set {p[3]} to {(+1,+1)}, {(-1,+1)}, {(+1,-1)}, or {(1/3,1/3)} */
            /* based on inverse map of {pM_exp[3]}: */
            r2_t u = hr2_pmap_inv_r2_point(&(pM_exp[3]), &M);
            double tw = 1.0 - u.c[0] - u.c[1];
            double tx = u.c[0];
            double ty = u.c[1];
            if (verbose)
              { r2_gen_print(Er, &u, "%+20.16f", "  u = ( ", " ", " )\n");
                Pr(Er, "  t = [ %+14.11f %+14.11f %+14.11f ]\n", tw, tx, ty);
              }
            if (((tw > 0) && (tx > 0) && (ty > 0)) || ((tw < 0) && (tx < 0) && (ty < 0)))
              { p[3] = (r2_t){{ +1/3.0, +1/3.0 }}; }
            else if (tx*ty > 0)
              { p[3] = (r2_t){{ +1.0, +1.0 }}; }
            else if (tw*tx > 0)
              { p[3] = (r2_t){{ +1.0, -1.0 }}; }
            else if (tw*ty > 0)
              { p[3] = (r2_t){{ -1.0, +1.0 }}; }
            else 
              { Pr(Er, "** bad u\n"); }
          }

        if (verbose)
          { Pr(Er, "  goal:\n");
            for (uint32_t kp = 0;  kp < 4; kp++) 
              { print_pair(&(p[kp]), &(pM_exp[kp])); }
            Pr(Er, "  actual:\n");
            for (uint32_t kp = 0;  kp < 4; kp++) 
              { r2_t pM_cmp = hr2_pmap_r2_point(&(p[kp]), &M);
                print_pair(&(p[kp]), &(pM_cmp));
              }
            Pr(Er, "\n");
          }

        /* Check whether the map works: */
        for (uint32_t kp = 0;  kp < 4; kp++) 
          { char *lab = (char*[4]){ "p[0]", "p[1]", "p[2]", "p[3]" }[kp];
            hr2_pmap_test_tools_check_map_r2_point(lab, &(p[kp]), &M, FALSE, &(pM_exp[kp]), "{hr2_pmap_from_four_r2_points} failed (2)");
          }
      }
    return;

    void print_pair(r2_t *p, r2_t *pM)
      { r2_gen_print(Er, p, "%+7.4f", "    ( ", " ", " )");
        r2_gen_print(Er, pM, "%+14.10f", " -> ( ", " ", " )\n");
      }
  }

void test_hr2_pmap_r2_from_sign_class(bool_t verbose)
  {
    if (verbose) { Pr(Er, "=== testing hr2_pmap_r2_from_sign_class ==============================\n"); }
    
    /* The four test points: */
    r2_t p[4];
    p[0] = (r2_t){{ 0, 0 }};
    p[1] = (r2_t){{ 1, 0 }};
    p[2] = (r2_t){{ 0, 1 }};
    p[3] = (r2_t){{ 1.0/3.0, 1.0/3.0 }};
    char *lab[4] = { "(00)", "(10)", "(01)", "(tt)" };
    
    /* Image of {p[3]} for each {class}: */
    r2_t u[4];
    u[0] = p[3];
    u[1] = (r2_t){{ +1, -1 }};
    u[2] = (r2_t){{ -1, +1 }};
    u[3] = (r2_t){{ +1, +1 }};
    
    for (uint32_t class = 0;  class <= 3; class++)
      { hr2_pmap_t M = hr2_pmap_r2_from_sign_class(class);
        for (uint32_t kp = 0;  kp < 4; kp++)
          { r2_t q;
            if (kp != 3)
              { q = p[kp]; }
            else
              { q = u[class]; }
            hr2_pmap_test_tools_check_map_r2_point(lab[kp], &(p[kp]), &M, FALSE, &(q), "{hr2_pmap_r2_from_sign_class} failed (1)");
          }
      }
  }

void test_hr2_pmap_r2_point__hr2_pmap_inv_r2_point(bool_t verbose)
  {
    if (verbose) { Pr(Er, "=== testing hr2_pmap_r2_point, hr2_pmap_inv_r2_point,  ==============================\n"); }
    hr2_pmap_t M;  hr2_pmap_test_tools_throw_aff_map(&M);
    r2_t pc; r2_throw_cube(&pc);
    /* Compute {qc_exp} without using {hr2_pmap_r2_point}: */
    hr2_point_t ph = hr2_from_r2(&pc);
    hr2_point_t qh_exp = hr2_pmap_point(&ph, &M);
    r2_t qc_exp = r2_from_hr2(&qh_exp);
    hr2_pmap_test_tools_check_map_r2_point("p", &pc,     &M, FALSE, &qc_exp, "{hr2_pmap_r2_point} failed (1)");
    hr2_pmap_test_tools_check_map_r2_point("q", &qc_exp, &M, TRUE,  &pc,     "{hr2_pmap_inv_r2_point} failed (2)");
  }

void test_hr2_pmap_inv(bool_t verbose)
  {
    if (verbose) { Pr(Er, "=== testing hr2_pmap_inv ==============================\n"); }
    hr2_pmap_t M;  hr2_pmap_test_tools_throw_aff_map(&M);
    hr2_pmap_t N = hr2_pmap_inv(&M); 
    hr2_pmap_test_tools_check("N", &N, 1.0e-8, "{hr2_pmap_inv} failed (1)");
    for (uint32_t k = 0;  k < 5; k++)
      { hr2_point_t p = hr2_point_throw();
        hr2_point_t q = hr2_pmap_point(&p, &M);
        hr2_pmap_test_tools_check_map_point("q", &q, FALSE, FALSE, &N, FALSE, &p, "{hr2_pmap_inv} failed (2)");
      }
 }

void test_hr2_pmap_compose(bool_t verbose)
  {
    if (verbose) { Pr(Er, "=== testing hr2_pmap_compose ==============================\n"); }
    hr2_pmap_t M;  hr2_pmap_test_tools_throw_aff_map(&M);
    hr2_pmap_t N;  hr2_pmap_test_tools_throw_aff_map(&N);
    hr2_pmap_t P = hr2_pmap_compose(&M, &N);
    for (uint32_t k = 0;  k < 5; k++)
      { hr2_point_t p = hr2_point_throw();
        hr2_point_t q = hr2_pmap_point(&p, &M);
        hr2_point_t r = hr2_pmap_point(&q, &N);
        hr2_pmap_test_tools_check_map_point("p", &p, FALSE, FALSE, &P, FALSE, &r, "{hr2_pmap_compose} failed (1)");
      }
   }

void test_hr2_pmap_inv_compose(bool_t verbose)
  {
    bool_t debug = TRUE;
    if (debug) { verbose = TRUE; }
    
    if (verbose) { Pr(Er, "=== testing hr2_pmap_inv_compose ==============================\n"); }
    
    if (verbose) { Pr(Er, "  generating a random affine map {M} ...\n"); }
    hr2_pmap_t M;  hr2_pmap_test_tools_throw_aff_map(&M);
    hr2_pmap_test_tools_check("M", &M, 1.0e-8, "{hr2_pmap_test_tools_throw_aff_map} failed (1)");
    
    if (verbose) { Pr(Er, "  getting the inverse {Minv} of map {M} ...\n"); }
    hr2_pmap_t Minv = hr2_pmap_inv(&M); 
    hr2_pmap_test_tools_check("Minv", &Minv, 1.0e-8, "{hr2_pmap_inv} failed");
    
    if (verbose) { Pr(Er, "  generating another random affine map {N} ...\n"); }
    hr2_pmap_t N;  hr2_pmap_test_tools_throw_aff_map(&N);
    hr2_pmap_test_tools_check("N", &N, 1.0e-8, "{hr2_pmap_test_tools_throw_aff_map} failed (2)");
    
    if (verbose) { Pr(Er, "  computing the composite {P = Minv N} ...\n"); }
    hr2_pmap_t P = hr2_pmap_inv_compose(&M, &N); 
    hr2_pmap_test_tools_check("P", &P, 1.0e-8, "{hr2_pmap_inv_compose} failed");
    
    if (verbose) { Pr(Er, "  testing {P} with some random points ...\n"); }
    for (uint32_t k = 0;  k < 5; k++)
      { hr2_point_t p = hr2_point_throw();
        hr2_point_t q = hr2_pmap_point(&p, &Minv);
        hr2_point_t r = hr2_pmap_point(&q, &N);
        hr2_pmap_test_tools_check_map_point("p", &p, FALSE, FALSE, &P, FALSE, &r, "{hr2_pmap_inv_compose} failed (2)");
      }
   }

void test_hr2_pmap_translation_from_disp(bool_t verbose)
  {
    if (verbose) { Pr(Er, "=== testing hr2_pmap_translation_from_disp ==============================\n"); }
    r2_t disp; r2_throw_cube(&disp);
    hr2_pmap_t M = hr2_pmap_translation_from_disp(&disp); 
    hr2_pmap_test_tools_check("M", &M, 1.0e-8, "{hr2_pmap_translation_from_disp} failed (1)");
    for (uint32_t k = 0;  k < 5; k++)
      { r2_t p; r2_throw_cube(&p);
        r2_t q; r2_add(&p, &disp, &q);
        hr2_pmap_test_tools_check_map_r2_point("p", &p, &M, FALSE, &q, "{hr2_pmap_translation_from_disp} failed (2)");
      }
  }

void test_hr2_pmap_rotation(bool_t verbose)
  {
    if (verbose) { Pr(Er, "=== testing hr2_pmap_rotation ==============================\n"); }
    double ang = dabrandom(-1.58, +1.58);
    hr2_pmap_t M = hr2_pmap_rotation(ang); 
    hr2_pmap_test_tools_check("M", &M, 1.0e-8, "{hr2_pmap_rotation} failed (1)");
    double ca = cos(ang);
    double sa = sin(ang);
    for (uint32_t k = 0;  k < 5; k++)
      { r2_t p; r2_throw_cube(&p);
        r2_t q;
        q.c[0] = + ca*p.c[0] - sa*p.c[1];
        q.c[1] = + sa*p.c[0] + ca*p.c[1];
        hr2_pmap_test_tools_check_map_r2_point("p", &p, &M, FALSE, &q, "{hr2_pmap_rotation} failed (2)");
      }
  }

void test_hr2_pmap_rotation_and_scaling(bool_t verbose)
  {
    if (verbose) { Pr(Er, "=== testing hr2_pmap_rotation_and_scaling ==============================\n"); }
    double ang = dabrandom(-1.58, +1.58);
    double scale = dabrandom(0.5, 2.0);
    hr2_pmap_t M = hr2_pmap_rotation_and_scaling(ang, scale); 
    hr2_pmap_test_tools_check("M", &M, 1.0e-8, "{hr2_pmap_rotation_and_scaling} failed (1)");
    double ca = cos(ang);
    double sa = sin(ang);
    for (uint32_t k = 0;  k < 5; k++)
      { r2_t p; r2_throw_cube(&p);
        r2_t q;
        q.c[0] = + ca*scale*p.c[0] - sa*scale*p.c[1];
        q.c[1] = + sa*scale*p.c[0] + ca*scale*p.c[1];
        hr2_pmap_test_tools_check_map_r2_point("p", &p, &M, FALSE, &q, "{hr2_pmap_rotation_and_scaling} failed (2)");
      }
  }

void test_hr2_pmap_affine_from_three_points(bool_t verbose)
  {
    if (verbose) { Pr(Er, "=== testing hr2_pmap_affine_from_three_points ==============================\n"); }
    r2_t o; r2_throw_cube(&o);
    r2_t p; r2_throw_cube(&p);
    r2_t q; r2_throw_cube(&q);
    hr2_pmap_t M = hr2_pmap_affine_from_three_points(&o, &p, &q); 
    hr2_pmap_test_tools_check("M", &M, 1.0e-8, "{hr2_pmap_affine_from_three_points} failed (1)");
    /* Check whether the map works: */
    r2_t oo = (r2_t){{ 0.0, 0.0 }};
    r2_t pp = (r2_t){{ 1.0, 0.0 }};
    r2_t qq = (r2_t){{ 0.0, 1.0 }};

    hr2_pmap_test_tools_check_map_r2_point("o", &oo, &M, FALSE, &o, "{hr2_pmap_affine_from_three_points} failed (2)");
    hr2_pmap_test_tools_check_map_r2_point("p", &pp, &M, FALSE, &p, "{hr2_pmap_affine_from_three_points} failed (3)");
    hr2_pmap_test_tools_check_map_r2_point("q", &qq, &M, FALSE, &q, "{hr2_pmap_affine_from_three_points} failed (4)");
  }

void test_hr2_pmap_point(bool_t inv, bool_t verbose)
  { char *tag = (inv ? "_inv" : "");
    if (verbose) { Pr(Er, "=== testing hr2_pmap%s_point ==============================\n", tag); }
    hr2_point_t p = hr2_point_throw();
    hr2_pmap_t M; hr2_pmap_throw(&M);

    /* Compute expected images of {p} by {M.dir} and {M.inv}: */
    hr2_point_t q_exp;
    r3x3_map_row(&(p.c), (inv ? &(M.inv) : &(M.dir)), &(q_exp.c));

    /* Compare with library func results: */
    hr2_pmap_test_tools_check_map_point("p", &p, FALSE, FALSE, &M, inv, &q_exp, txtcat3("hr2_pmap", tag,  "_point failed (1)"));
  }

void test_hr2_pmap_line(bool_t inv, bool_t verbose)
  { 
    char *tag = (inv ? "_inv" : "");
    if (verbose) { Pr(Er, "=== testing hr2_pmap%s_line ==============================\n", tag); }
    hr2_line_t A = hr2_line_throw();
    hr2_pmap_t M; hr2_pmap_throw(&M);

    /* Compute expected images of {p} by {M.dir} and {M.inv}: */
    hr2_line_t B_exp;
    r3x3_map_col((inv ? &(M.dir) : &(M.inv)), &(A.f), &(B_exp.f));

    /* Compare with library func results: */
    hr2_pmap_test_tools_check_map_line("A", &A, FALSE, FALSE, &M, inv, &B_exp, txtcat3("hr2_pmap", tag, "_line failed (1)"));

    /* Pick a point on the line: */
    hr2_point_t p;
    (void)r3_throw_ortho(&(A.f), &(p.c));
    hr2_point_t q = (inv ? hr2_pmap_inv_point(&p, &M) : hr2_pmap_point(&p, &M));
    double dot = r3_dot(&(B_exp.f), &(q.c));
    affirm(fabs(dot) < 1.0e-11, txtcat3(txtcat3("hr2_pmap", tag, "_line"), " incompat with ", txtcat3("hr2_pmap", tag, "_point")));
  }

void test_hr2_pmap_identity(bool_t verbose)
  { if (verbose) { Pr(Er, "!! {hr2_pmap_identity} NOT TESTED\n"); } }

void test_hr2_pmap_mirror(bool_t verbose)
  { if (verbose) { Pr(Er, "!! {hr2_pmap_mirror} NOT TESTED\n"); } }

void test_hr2_pmap_print(bool_t verbose)
  { if (verbose) { Pr(Er, "!! {hr2_pmap_print} NOT TESTED\n"); } }

void test_hr2_pmap_gen_print(bool_t verbose)
  { if (verbose) { Pr(Er, "!! {hr2_pmap_gen_print} NOT TESTED\n"); } }

void test_hr2_pmap_is_valid(bool_t verbose)
  { 
    bool_t debug = FALSE;
    if (debug) { verbose = TRUE; }
    if (verbose) { Pr(Er, "=== testing hr2_pmap_is_valid ==============================\n"); }

    double tol = 1.0e-10;  

    /* Create a barely invalid map {M}: */
    double detMaxA = 0.1*tol; /* Significantly smaller than {tol}. */
    hr2_pmap_t MA = hr2_pmap_test_tools_throw_almost_singular(detMaxA);

    /* Create a barely valid map {N}: */
    double detMinB = 10.0*tol; /* Significantly bigger than {tol}. */
    hr2_pmap_t MB = hr2_pmap_test_tools_throw_non_singular(detMinB);

    if (verbose) { Pr(Er, "  testing {MA}, should be invalid...\n"); }
    demand(! hr2_pmap_is_valid(&MA, tol), "{hr2_pmap_is_valid(MA)} true, should be false" );

    if (verbose) { Pr(Er, "  testing {MB}, should be valid...\n"); }
    demand(hr2_pmap_is_valid(&MB, tol),  "{hr2_pmap_is_valid(MB)} false, should be true" );
  }

void test_hr2_pmap_sign__hr2_pmap_invert_sign__hr2_pmap_set_sign(bool_t verbose)
  { if (verbose) { Pr(Er, "=== testing hr2_pmap_sign, hr2_pmap_invert_sign, hr2_pmap_set_sign ==============================\n"); }
    
    hr2_pmap_t M;
    double detM;
    do { 
      hr2_pmap_throw(&M);
      detM = r3x3_det(&(M.dir));
    } while (detM < 0.001);
    demand(hr2_pmap_sign(&M) == +1, "{hr2_pmap_sign} failed");
    for (sign_t sgn = +1; sgn >= -1; sgn -= 2)
      { hr2_pmap_t N = M;
        hr2_pmap_set_sign(&N, sgn);
        double detN = r3x3_det(&(N.dir));
        demand(sgn*detN > 0, "{hr2_pmap_set_sign} failed");
        sign_t sgnN = hr2_pmap_sign(&N);
        demand(sgnN == sgn, "{hr2_pmap_sign} failed");
        hr2_pmap_t P = N;
        hr2_pmap_invert_sign(&P);
        sign_t sgnP = hr2_pmap_sign(&P);
        demand(sgnP == -sgn, "{hr2_pmap_invert_sign} failed");
        hr2_pmap_invert_sign(&P);
        for (uint32_t i = 0;  i < 3; i++)
          { for (uint32_t j = 0;  j < 3; j++)
             { demand(P.dir.c[i][j] == N.dir.c[i][j], "{hr2_invert_sign^2} failed (1)");
               demand(P.inv.c[i][j] == N.inv.c[i][j], "{hr2_invert_sign^2} failed (2)");
             }
          }
      }
  }

void test_hr2_pmap_is_identity(bool_t verbose)
  { 
    bool_t debug = FALSE;
    if (debug) { verbose = TRUE; }
    if (verbose) { Pr(Er, "=== testing hr2_pmap_is_identity ==============================\n"); }
    
    double tol = 1.0e-14;
 
    /* Perturbation small enough to keep {M} "practically identity" and {N} definitely not identity: */
    double eps = 0.1*tol;
    r3x3_t P = (r3x3_t){{ { eps, eps, eps }, { eps, eps, eps }, { eps, eps, eps } }};

    auto bool_t ok(hr2_pmap_t *R);
      /* Returns {hr2_pmap_is_identity(R, tol)}. */
   
    if (verbose) { Pr(Er, "  creating a map {M} that is the identity\n"); }
    hr2_pmap_t M; 
    r3x3_ident(&(M.dir));
    r3x3_ident(&(M.inv));
    thr2_random_scale_and_print("M", &M, 2, verbose);

    if (verbose) { Pr(Er, "  testing {M} and perturbed versions, should be ok ...\n"); }
    hr2_pmap_test_tools_check_perturbed("M", &M, &ok, &P, TRUE, "hr2_pmap_is_identity", debug);
    
    if (verbose) { Pr(Er, "  creating a map {N} that is a translation by a non-trivial vector\n"); }
    hr2_pmap_t N; 
    r3x3_ident(&(N.dir));
    N.dir.c[0][1] = 10*tol; N.dir.c[0][2] = 20*tol;
    r3x3_inv(&(N.dir), &(N.inv)); 
    thr2_random_scale_and_print("N", &N, 2, verbose);
 
    if (verbose) { Pr(Er, "  testing {N} and perturbed versions, should not be ok ...\n"); }
    hr2_pmap_test_tools_check_perturbed("N", &N, &ok, &P, FALSE, "hr2_pmap_is_identity", debug);
    
    return;
    
    bool_t ok(hr2_pmap_t *R)
      { return hr2_pmap_is_identity(R, tol); }
  }

void test_hr2_pmap_is_translation(bool_t verbose)
  { 
    bool_t debug = FALSE;
    if (debug) { verbose = TRUE; }
    if (verbose) { Pr(Er, "=== testing hr2_pmap_is_translation ==============================\n"); }
    
    double tol = 1.0e-14;
    
    auto hr2_pmap_t throw_congruence(char *name, double ang);
      /* Throws a map {M} that is a translation by a random significant vector
        after a rotation by angle {ang}. The matrices {M.dir} and {M.inv} are scaled
        random homogeneous factors. */

    auto bool_t ok(hr2_pmap_t *R);
      /* Returns {hr2_pmap_is_translation(R, tol). */

    /* Perturbation small enough to keep {M} "practically translation", and {N} "not translation": */
    double eps = 0.1*tol;
    r3x3_t P = (r3x3_t){{ { eps, eps, eps }, { eps, eps, eps }, { eps, eps, eps } }};
    
    if (verbose) { Pr(Er, "  creating {M}, effectively a non-trivial translation\n"); }
    double ang_M = 0.05*tol;
    hr2_pmap_t M = throw_congruence("M", ang_M);

    if (verbose) { Pr(Er, "  testing {M} and perturbed versions, should be ok ...\n"); }
    hr2_pmap_test_tools_check_perturbed("M", &M, &ok, &P, TRUE, "hr2_pmap_is_translation", debug);
    
    if (verbose) { Pr(Er, "  creating {N}, a translation after a small but signif rotation\n"); }
    double ang_N = 3.0*tol;
    hr2_pmap_t N = throw_congruence("N", ang_N);
 
    if (verbose) { Pr(Er, "  testing {N} and perturbed versions, should not be ok ...\n"); }
    hr2_pmap_test_tools_check_perturbed("N", &N, &ok, &P, FALSE, "hr2_pmap_is_identity", debug);

    return;

    hr2_pmap_t throw_congruence(char *name, double ang)
      { double tmag = dabrandom(5.0,10.0);
        r2_t tvec; r2_throw_dir(&tvec); r2_scale(tmag, &tvec, &tvec);
        hr2_pmap_t R = hr2_pmap_from_parms
          ( /*persp*/ NULL, /*shear */ 0, /*skew */ 0, 
            /*scale*/ 1.0, /*ang*/ ang, /*disp*/ &tvec
          );
        thr2_random_scale_and_print(name, &R, 2, verbose);
        return R;
      }
   
    bool_t ok(hr2_pmap_t *R)
      { return hr2_pmap_is_translation(R, tol); }
  }

void test_hr2_pmap_is_congruence(bool_t verbose)
  { bool_t debug = TRUE;
    if (debug) { verbose = TRUE; }
    
    if (verbose) { Pr(Er, "=== testing hr2_pmap_is_congruence ===\n"); }
    
    double tol = 2.0e-13;
    if (verbose) { Pr(Er, "  tol = %24.16e\n", tol); }

    /* Perturbation small enough to keep {M} "practically ok", and {N} "not ok": */
    double eps = 0.01*tol;
    r3x3_t P = (r3x3_t){{ { eps, 0.0, 0.0 }, { eps, eps, eps }, { eps, eps, eps } }};
    
    auto hr2_pmap_t throw_similarity(char *name, double scale);
      /* Throws a map {M} that is a congruence with a random
        displacement and angle after a uniform scaling by {scale}. The
        matrices {M.dir} and {M.inv} are scaled random homogeneous
        factors. */

    auto bool_t ok(hr2_pmap_t *R);
      /* Returns {hr2_pmap_is_congruence(R,tol)}. */
    
    if (verbose)  { Pr(Er, "  creating {M}, definitely congruence not translation\n"); }
    double sgn_sc_M = (drandom() < 0.5 ? -1 : +1);
    double scale_M = 1.0 + sgn_sc_M*dabrandom(0.03, 0.05)*tol;
    hr2_pmap_t M = throw_similarity("M", scale_M);

    if (verbose) { Pr(Er, "  testing {M} and perturbed versions, should be ok ...\n"); }
    hr2_pmap_test_tools_check_perturbed("M", &M, &ok, &P, TRUE, "hr2_pmap_is_congruence", debug);

    if (verbose) { Pr(Er, "  creating {N}, that is almost but definitely not a congruence\n"); }
    double sgn_sc_N = (drandom() < 0.5 ? -1 : +1);
    double scale_N = 1.0 + sgn_sc_N*dabrandom(3.0, 5.0)*tol;
    hr2_pmap_t N = throw_similarity("N", scale_N);
 
    if (verbose) { Pr(Er, "  testing {N} and perturbed versions, should not be ok ...\n"); }
    hr2_pmap_test_tools_check_perturbed("N", &N, &ok, &P, FALSE, "hr2_pmap_is_congruence", debug);

    return;
    
    hr2_pmap_t throw_similarity(char *name, double scale)
      { 
        double tmag = dabrandom(5.0, 10.0);
        r2_t tvec; r2_throw_dir(&tvec); r2_scale(tmag, &tvec, &tvec);
        double ang = 2*M_PI*(0.1 + 0.8*drandom()); /* Rotation angle. */
        hr2_pmap_t R = hr2_pmap_from_parms
          ( /*persp*/ NULL, /*shear */ 0, /*skew */ 0, 
            /*scale*/ scale, /*ang*/ ang, /*disp*/ &tvec
          );
        thr2_random_scale_and_print(name, &R, 2, verbose);
        return R;
      }
      
    bool_t ok(hr2_pmap_t *R)
      { return hr2_pmap_is_congruence(R, tol); }

  }

void test_hr2_pmap_is_similarity(bool_t verbose)
  { bool_t debug = TRUE;
    if (debug) { verbose = TRUE; }
    
    if (verbose) { Pr(Er, "=== testing hr2_pmap_is_similarity ===\n"); }
    
    double tol = 2.0e-13;
    if (verbose) { Pr(Er, "  tol = %24.16e\n", tol); }

    /* Perturbation small enough to keep {M} "practically ok", and {N} "not ok": */
    double eps = 0.01*tol;
    r3x3_t P = (r3x3_t){{ { eps, 0.0, 0.0 }, { eps, eps, eps }, { eps, eps, eps } }};
    
    auto hr2_pmap_t throw_affine(char *name, double shear, double skew, double scale);
      /* Generates an affine map {M} that is a shear and non-unif
        scaling with parameters {shear} and {skew}, follwed by uniform
        scaling by {scale}, rotation by a non-trivial random angle, and
        displacement by a non-trivial distance. The matrices {M.dir} and
        {M.inv} are scaled random homogeneous factors. */

    double sgn_sc = (drandom() < 0.5 ? -1 : +1);
    double scale_gen = 1.0 + sgn_sc*dabrandom(3.0, 5.0)*tol;  /* Used to create the maps.  */
    double scale_check = (drandom() < 0.5 ? scale_gen : NAN);  /* Used when checking.  */

    auto bool_t ok(hr2_pmap_t *R);
      /* Returns {hr2_pmap_is_similarity(R,scale_check,tol)}. */
    
    if (verbose)  { Pr(Er, "  creating {M}, definitely similarity not congruence\n"); }
    double shear_M = (drandom() < 0.5 ? -1 : +1)*dabrandom(0.03, 0.05)*tol;
    double skew_M = exp((drandom() < 0.5 ? -1 : +1)*log(1 + dabrandom(0.03, 0.05)*tol));
    hr2_pmap_t M = throw_affine("M", shear_M, skew_M, scale_gen);

    if (verbose) { Pr(Er, "  testing {M} and perturbed versions, should be ok ...\n"); }
    hr2_pmap_test_tools_check_perturbed("M", &M, &ok, &P, TRUE, "hr2_pmap_is_similarity", debug);

    if (verbose) { Pr(Er, "  creating {N}, affinity that is almost but definitely not a similarity\n"); }
    double shear_N = (drandom() < 0.5 ? -1 : +1)*dabrandom(3.0, 5.0)*tol;
    double skew_N = exp((drandom() < 0.5 ? -1 : +1)*log(1 + dabrandom(3.0, 5.0)*tol));
    hr2_pmap_t N = throw_affine("N", shear_N, skew_N, scale_gen);
 
    if (verbose) { Pr(Er, "  testing {N} and perturbed versions, should not be ok ...\n"); }
    hr2_pmap_test_tools_check_perturbed("N", &N, &ok, &P, FALSE, "hr2_pmap_is_similarity", debug);

    return;

    hr2_pmap_t throw_affine(char *name, double shear, double skew, double scale)
      { 
        double tmag = dabrandom(5.0, 10.0);
        r2_t tvec; r2_throw_dir(&tvec); r2_scale(tmag, &tvec, &tvec);
        double ang = 2*M_PI*(0.1 + 0.8*drandom()); /* Rotation angle. */
        hr2_pmap_t R = hr2_pmap_from_parms
          ( /*persp*/ NULL, /*shear */ shear, /*skew */ skew, 
            /*scale*/ scale, /*ang*/ ang, /*disp*/ &tvec
          );
        thr2_random_scale_and_print(name, &R, 2, verbose);

        return R;
      }

    bool_t ok(hr2_pmap_t *R)
      { return hr2_pmap_is_similarity(R, scale_check, tol); }
  }

void test_hr2_pmap_is_affine(bool_t verbose)
  { 
    bool_t debug = TRUE;
    if (debug) { verbose = TRUE; }
    if (verbose) { Pr(Er, "=== testing hr2_pmap_is_affine ===\n"); }
    
    double tol = 1.0e-14;
    if (verbose) { Pr(Er, "  tol = %24.16e\n", tol); }
   
    /* Perturbation small enough to keep {M} "practically ok", and {N} "not ok": */
    double eps = 0.05*tol;
    r3x3_t P = (r3x3_t){{ { 0.0, 0.0, 0.0 }, { eps, 0.0, 0.0 }, { eps, 0.0, 0.0 } }};
         
    auto hr2_pmap_t throw_persp(char *name, double pmag);
      /* Generates a projective map that is the specified amount of perspective
        distortion about the origin, followed by a random affine map
        with non-trivial random shear, skew, scaling, rotation angle,
        displacement. The
        matrices {M.dir} and {M.inv} are scaled random homogeneous
        factors. */

    auto bool_t ok(hr2_pmap_t *R);
      /* Returns {hr2_pmap_is_affine(R,tol)}. */
    
    if (verbose) { Pr(Er, "  creating {M}, essentially affine not similarity\n"); }
    double pmag_M = dabrandom(0.03, 0.05)*tol;
    hr2_pmap_t M = throw_persp("M", pmag_M);
 
    if (verbose) { Pr(Er, "  testing {N} and perturbed versions, should not be ok ...\n"); }
    hr2_pmap_test_tools_check_perturbed("M", &M, &ok, &P, TRUE, "hr2_pmap_is_affine", debug);

    if (verbose) { Pr(Er, "  creating {N}, almost but definitely nor affine\n");  }
    double pmag_N = dabrandom(30.0, 50.0)*tol;
    hr2_pmap_t N = throw_persp("N", pmag_N);
 
    if (verbose) { Pr(Er, "  testing {N} and perturbed versions, should not be ok ...\n"); }
    hr2_pmap_test_tools_check_perturbed("N", &N, &ok, &P, FALSE, "hr2_pmap_is_affine", debug);

    return;
    
    hr2_pmap_t throw_persp(char *name, double pmag)
      { double tmag = dabrandom(5.0, 10.0);
        r2_t tvec; r2_throw_dir(&tvec); r2_scale(tmag, &tvec, &tvec);
        double scale = (drandom() < 0.5 ? dabrandom(0.1, 0.4) :  dabrandom(3.0, 5.0));
        double ang = 2*M_PI*(0.1 + 0.8*drandom()); /* Rotation angle. */
        double shear = (drandom() < 0.5 ? -1 : +1)*dabrandom(0.2, 0.9);
        double skew = exp((drandom() < 0.5 ? -1 : +1)*log(1 + dabrandom(0.2, 2.9)));
        r2_t pvec; r2_throw_dir(&pvec); r2_scale(pmag, &pvec, &pvec);
        hr2_pmap_t R = hr2_pmap_from_parms
          ( /*persp*/ &pvec, /*shear */ shear, /*skew */ skew, 
            /*scale*/ scale, /*ang*/ ang, /*disp*/ &tvec
          );
        thr2_random_scale_and_print(name, &R, 2, verbose);
        return R;
      }
    bool_t ok(hr2_pmap_t *R)
      { return hr2_pmap_is_affine(R, tol); }
  }
  
void thr2_random_scale_and_print(char *Rname, hr2_pmap_t *R, uint32_t indent, bool_t verbose)
  { if (verbose)
      { char *Rname_raw = jsprintf("%s (raw)", Rname);
        hr2_pmap_test_tools_print(indent, Rname_raw, R);
        free(Rname_raw);
      }
    r3x3_scale(0.1+5*drandom(), &(R->dir), &(R->dir));
    r3x3_scale(0.1+5*drandom(), &(R->inv), &(R->inv));
    if (verbose)
      { char *Rname_scaled = jsprintf("%s (scaled)", Rname);
        hr2_pmap_test_tools_print(indent, Rname_scaled, R);
        free(Rname_scaled);
      }
  }
