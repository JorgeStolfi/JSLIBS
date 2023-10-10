/* hr3test --- test program for hr3.h  */
/* Last edited on 2023-10-09 21:29:09 by stolfi */

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
#include <r4_extra.h>
#include <r2x2.h>
#include <r2.h>

#include <hr3.h>
#include <hr3test_tools.h>

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_hr3(bool_t verbose);
void test_hr3_pmap(bool_t verbose);
void test_hr3_pmap_aff(bool_t verbose);

void test_hr4_to_from_r3(bool_t verbose);
void test_hr3_pt_pt_diff(bool_t verbose);
void test_hr3_side(bool_t verbose);
void test_hr3_orient(bool_t verbose);
void test_hr3_plane_from_three_points(bool_t verbose);
void test_hr3_point_from_three_planes(bool_t verbose);
void test_hr3_point_point_dir(bool_t verbose);
void test_hr3_plane_normal(bool_t verbose);

void test_hr3_pmap_gen_print(bool_t verbose);
void test_hr3_pmap_identity(bool_t verbose);
void test_hr3_pmap_inv(bool_t verbose);
void test_hr3_pmap_compose(bool_t verbose);
void test_hr3_pmap_inv_comp(bool_t verbose);
void test_hr3_pmap_is_affine(bool_t verbose);
void test_hr3_pmap_is_identity(bool_t verbose);
void test_hr3_pmap_plane(bool_t inv, bool_t verbose);
void test_hr3_pmap_point(bool_t inv, bool_t verbose);
void test_hr3_pmap_print(bool_t verbose);
void test_hr3_pmap_r3_point(bool_t verbose);

void test_hr3_pmap_aff_from_mat_and_disp(bool_t verbose);
void test_hr3_pmap_aff_from_four_points(bool_t verbose);
void test_hr3_pmap_u_v_rotation(bool_t verbose);
void test_hr3_pmap_translation(bool_t verbose);
void test_hr3_pmap_scaling(bool_t verbose);
void test_hr3_pmap_from_five_points(bool_t verbose);

int32_t main (int32_t argc, char **argv)
  {
    int32_t i;
    srand(1993);
    srandom(1933);

    for (i = 0; i < 100; i++) test_hr3(i < 3);
    for (i = 0; i < 100; i++) test_hr3_pmap(i < 3);
    for (i = 0; i < 100; i++) test_hr3_pmap_aff(i < 3);
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

void test_hr3_pmap(bool_t verbose)
  {

    /* Size: */
    if (verbose)
      { fprintf(stderr,
          "sizeof(hr3_pmap_t) = %lu  2*%d*%d*sizeof(double) = %lu\n",
          sizeof(hr3_pmap_t), NH, NH, 2*NH*NH*sizeof(double)
        );
      }
    
    test_hr3_pmap_point(FALSE, verbose);
    test_hr3_pmap_point(TRUE, verbose);
    test_hr3_pmap_plane(FALSE, verbose);
    test_hr3_pmap_plane(TRUE, verbose);
    test_hr3_pmap_r3_point(verbose);

    test_hr3_pmap_is_identity(verbose);
    test_hr3_pmap_identity(verbose);
    test_hr3_pmap_inv(verbose);
    test_hr3_pmap_compose(verbose);
    test_hr3_pmap_inv_comp(verbose);
    
    test_hr3_pmap_from_five_points(verbose);

    test_hr3_pmap_print(verbose);
    test_hr3_pmap_gen_print(verbose);
  }  

void test_hr3_pmap_aff(bool_t verbose)
  {
    test_hr3_pmap_is_affine(verbose);

    test_hr3_pmap_aff_from_mat_and_disp(verbose);

    test_hr3_pmap_translation(verbose);
    test_hr3_pmap_u_v_rotation(verbose);
    test_hr3_pmap_scaling(verbose);
    test_hr3_pmap_aff_from_four_points(verbose);

    /* TO BE COMPLETED !!! */
    
    if (verbose)
      { 
        /* fprintf(stderr, "!! r3_aff_map_from_XXX NOT TESTED\n"); */
      }
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
      for (int32_t i = 1; i <= NC; i++)
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
      h3tt_check_eps(dob, dex, 1.0e-8, "hr3_pt_pt_diff error(3)");
      /* Check invariance under rotations: */
      for (int32_t i = 0; i < NH; i++)
        { int32_t j = (i + 1) % NH; /* Another axis. */
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
          h3tt_check_eps(drt, dex, 1.0e-8, "hr3_pt_pt_diff error(4)");
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
    for (int32_t i = 0; i < NH; i++)
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
      for (int32_t i = 0; i < NC; i++)
        { h3tt_check_eps(upq.c[i], vpq.c[i], 1.0e-12, "hr3_point_point_dir error"); }
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
      for (int32_t i = 0; i < NC; i++)
        { h3tt_check_eps(nL.c[i], mL.c[i], tol, "hr3_plane_normal error"); }
    }
  }

void test_hr3_pmap_from_five_points(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_pmap_from_five_points ---\n"); }
    
    bool_t debug = FALSE;
        
    auto void print_pair(hr3_point_t *p, hr3_point_t *pM);
    
    hr3_point_t p = (hr3_point_t){{{ 1.0, 0.0, 0.0, 0.0 }}};
    hr3_point_t q = (hr3_point_t){{{ 0.0, 1.0, 0.0, 0.0 }}};
    hr3_point_t r = (hr3_point_t){{{ 0.0, 0.0, 1.0, 0.0 }}};
    hr3_point_t s = (hr3_point_t){{{ 0.0, 0.0, 0.0, 1.0 }}};
    hr3_point_t u;

    for (int32_t kt = 0; kt < 2*(1<<NH); kt++)
      { hr3_point_t pM_exp, qM_exp, rM_exp, sM_exp, uM_exp;
        u = (hr3_point_t){{{ 1.0, 1.0, 1.0 }}};
        if (kt < (1 << NH))
          { if (kt & 1) { u.c.c[0] = -u.c.c[0]; }
            if (kt & 2) { u.c.c[1] = -u.c.c[1]; }
            if (kt & 4) { u.c.c[2] = -u.c.c[2]; }
            if (kt & 8) { u.c.c[3] = -u.c.c[3]; }
            pM_exp = p; qM_exp = q; rM_exp = r; sM_exp = s; uM_exp = u;
          }
        else
          { r4_throw_cube(&(pM_exp.c));
            r4_throw_cube(&(qM_exp.c));
            r4_throw_cube(&(rM_exp.c));
            r4_throw_cube(&(sM_exp.c));
            r4_throw_cube(&(uM_exp.c));
          }

        hr3_pmap_t M = hr3_pmap_from_five_points(&pM_exp, &qM_exp, &rM_exp, &sM_exp, &uM_exp);
        
        if (kt >= (1 << NH))
          { /* Choose {u = [±1,±1,±1]} based on inverse map of {uM_exp}: */
            u = hr3_pmap_inv_point(&uM_exp, &M);
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
            print_pair(&s, &sM_exp);
            print_pair(&u, &uM_exp);
            fprintf(stderr, "  actual:\n");
            hr3_point_t pM_cmp = hr3_pmap_point(&p, &M); print_pair(&p, &pM_cmp);
            hr3_point_t qM_cmp = hr3_pmap_point(&q, &M); print_pair(&q, &qM_cmp);
            hr3_point_t rM_cmp = hr3_pmap_point(&r, &M); print_pair(&r, &rM_cmp);
            hr3_point_t sM_cmp = hr3_pmap_point(&s, &M); print_pair(&s, &sM_cmp);
            hr3_point_t uM_cmp = hr3_pmap_point(&u, &M); print_pair(&u, &uM_cmp);
          }

        /* Check whether the map works: */
        h3tt_check_pmap_point("p", &p, FALSE, &M, FALSE, &pM_exp, "hr3_pmap_from_five_points failed");
        h3tt_check_pmap_point("q", &q, FALSE, &M, FALSE, &qM_exp, "hr3_pmap_from_five_points failed");
        h3tt_check_pmap_point("r", &r, FALSE, &M, FALSE, &rM_exp, "hr3_pmap_from_five_points failed");
        h3tt_check_pmap_point("s", &s, FALSE, &M, FALSE, &sM_exp, "hr3_pmap_from_five_points failed");
        h3tt_check_pmap_point("u", &u, FALSE, &M, FALSE, &uM_exp, "hr3_pmap_from_five_points failed");
      }
    return;
            
    void print_pair(hr3_point_t *p, hr3_point_t *pM)
      { r4_gen_print(stderr, &(p->c), "%+4.1f", "    [ ", " ", " ]");
        r4_gen_print(stderr, &(pM->c), "%+14.10f", " -> [ ", " ", " ]\n");
      }
  }
    
void test_hr3_pmap_r3_point(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr3_pmap_r3_point, hr3_pmap_inv_r3_point,  ---\n"); }
    hr3_pmap_t M;  h3tt_throw_aff_map(&M);
    r3_t pc; r3_throw_cube(&pc);
    /* Compute {qc_exp} without using {hr3_pmap_r3_point}: */
    hr3_point_t ph = hr3_from_r3(&pc);
    hr3_point_t qh_exp = hr3_pmap_point(&ph, &M);
    r3_t qc_exp = r3_from_hr3(&qh_exp);
    h3tt_check_pmap_r3_point("p", &pc,     &M, FALSE, &qc_exp, "hr3_pmap_r3_point failed");
    h3tt_check_pmap_r3_point("q", &qc_exp, &M, TRUE,  &pc,     "hr3_pmap_inv_r3_point failed");
  }

void test_hr3_pmap_inv(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr3_pmap_inv ---\n"); }
    hr3_pmap_t M;  h3tt_throw_aff_map(&M);
    hr3_pmap_t N = hr3_pmap_inv(&M);
    for (int32_t k = 0; k < 5; k++)
      { hr3_point_t p = hr3_point_throw();
        hr3_point_t q = hr3_pmap_point(&p, &M);
        h3tt_check_pmap_point("q", &q, FALSE, &N, FALSE, &p, "hr3_pmap_inv failed");
      }
 }

void test_hr3_pmap_compose(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr3_pmap_compose ---\n"); }
    hr3_pmap_t M;  h3tt_throw_aff_map(&M);
    hr3_pmap_t N;  h3tt_throw_aff_map(&N);
    hr3_pmap_t P = hr3_pmap_compose(&M, &N);
    for (int32_t k = 0; k < 5; k++)
      { hr3_point_t p = hr3_point_throw();
        hr3_point_t q = hr3_pmap_point(&p, &M);
        hr3_point_t r = hr3_pmap_point(&q, &N);
        h3tt_check_pmap_point("p", &p, FALSE, &P, FALSE, &r, "hr3_pmap_compose failed");
      }
   }

void test_hr3_pmap_inv_comp(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr3_pmap_inv_comp ---\n"); }
    hr3_pmap_t M;  h3tt_throw_aff_map(&M);
    hr3_pmap_t Minv = hr3_pmap_inv(&M);
    hr3_pmap_t N;  h3tt_throw_aff_map(&N);
    hr3_pmap_t P = hr3_pmap_inv_compose(&M, &N);
    for (int32_t k = 0; k < 5; k++)
      { hr3_point_t p = hr3_point_throw();
        hr3_point_t q = hr3_pmap_point(&p, &Minv);
        hr3_point_t r = hr3_pmap_point(&q, &N);
        h3tt_check_pmap_point("p", &p, FALSE, &P, FALSE, &r, "hr3_pmap_inv_comp failed");
      }
   }

void test_hr3_pmap_translation(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr3_pmap_translation ---\n"); }
    r3_t r; r3_throw_cube(&r);
    hr3_pmap_t M = hr3_pmap_translation(&r);
    for (int32_t k = 0; k < 5; k++)
      { r3_t p; r3_throw_cube(&p);
        r3_t q; r3_add(&r, &p, &q);
        h3tt_check_pmap_r3_point("p", &p, &M, FALSE, &q, "hr3_pmap_translation failed");
      }
  }

void test_hr3_pmap_u_v_rotation(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_pmap_u_v_rotation ---\n"); }
    r3_t u; r3_throw_dir(&u);
    r3_t v; r3_throw_dir(&v);
    hr3_pmap_t M = hr3_pmap_u_v_rotation(&u, &v);
    
    /* Check that it is a rotation: */
    for (int32_t d = 0; d < 2; d++)
      { r4x4_t *Q = (d == 0 ? &(M.dir) : &(M.inv));
        affirm(Q->c[0][0] > 0, "hr3_pmap_u_v_rotation failed - not affine (1)");
        for (int32_t i = 1; i < NH; i++)
          { affirm(Q->c[i][0] == 0, "hr3_pmap_u_v_rotation failed - not affine (2)");
            affirm(Q->c[0][i] == 0, "hr3_pmap_u_v_rotation failed - not linear");
            for (int32_t k = 1; k <= i; k++)
              { double dot = Q->c[i][0]*Q->c[k][0] + Q->c[i][1]*Q->c[k][1] + Q->c[i][2]*Q->c[k][2];
                affirm(fabs(dot - (i == k ? 1 : 0)) > 1.0e-11, "hr3_pmap_u_v_rotation failed - not orthonormal");
              }
          }
        affirm(r4x4_det(Q) > 0, "hr3_pmap_u_v_rotation failed - not orient preserving");
      }

    /* Check that it maps {u} to {v}: */
    h3tt_check_pmap_r3_point("u", &u, &M, FALSE, &v, "hr3_pmap_inv_comp failed");
    h3tt_check_pmap_r3_point("v", &v, &M, TRUE, &u, "hr3_pmap_inv_comp failed");
  }

void test_hr3_pmap_aff_from_mat_and_disp(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr3_pmap_aff_from_mat_disp ---\n"); }
    
    r3_t disp; r3_throw_cube(&disp);
    r3x3_t mat;
    for (int32_t i = 0; i < 2; i++)
      { r3_t s; r3_throw_cube(&s);
        mat.c[i][0] = s.c[0];
        mat.c[i][1] = s.c[1];
      }
    hr3_pmap_t A = hr3_pmap_aff_from_mat_and_disp(&mat, &disp);
    for (int32_t k = 0; k < 5; k++)
      { /* Should take the origin to {disp}* */
        r3_t o = (r3_t){{ 0.0, 0.0 }};
        h3tt_check_pmap_r3_point("o", &o, &A, FALSE, &disp, "hr3_pmap_aff_from_mat_disp failed");
        /* Test with a few other points: */
        for (int32_t kp = 0; kp < 3; kp++)
          { r3_t p; r3_throw_cube(&p);
            r3_t q; r3x3_map_row(&p, &mat, &q);
            r3_add(&disp, &q, &q);
            h3tt_check_pmap_r3_point("p", &p, &A, FALSE, &q, "hr3_pmap_aff_from_mat_disp failed");
          }
      }
  }

void test_hr3_pmap_scaling(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_pmap_scaling ---\n"); }

    r3_t scale; 
    do { r3_throw_cube(&scale); } while ((scale.c[0] == 0) || (scale.c[1] == 0));
    
    hr3_pmap_t A = hr3_pmap_scaling(&scale);

    /* Check whether the map works: */
    r3_t oc = (r3_t){{ 0.0, 0.0 }};
    r3_t pc = (r3_t){{ 1.0, 0.0 }};
    r3_t qc = (r3_t){{ 0.0, 1.0 }};
    r3_t uc = (r3_t){{ 1.0, 1.0 }};
    
    r3_t pcM = (r3_t){{ scale.c[0], 0.0 }};
    r3_t qcM = (r3_t){{ 0.0, scale.c[1] }};
    r3_t ucM = (r3_t){{ scale.c[0], scale.c[1] }};

    h3tt_check_pmap_r3_point("o", &oc, &A, FALSE, &oc, "hr3_pmap_scaling failed");
    h3tt_check_pmap_r3_point("p", &pc, &A, FALSE, &pcM, "hr3_pmap_scaling failed");
    h3tt_check_pmap_r3_point("q", &qc, &A, FALSE, &qcM, "hr3_pmap_scaling failed");
    h3tt_check_pmap_r3_point("u", &uc, &A, FALSE, &ucM, "hr3_pmap_scaling failed");
  }

void test_hr3_pmap_aff_from_four_points(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr3_pmap_aff_from_four_points ---\n"); }
    r3_t o; r3_throw_cube(&o);
    r3_t p; r3_throw_cube(&p);
    r3_t q; r3_throw_cube(&q);
    r3_t r; r3_throw_cube(&r);
    hr3_pmap_t M = hr3_pmap_aff_from_four_points(&o, &p, &q, &r);
    /* Check whether the map works: */
    r3_t oo = (r3_t){{ 0.0, 0.0, 0.0 }};
    r3_t pp = (r3_t){{ 1.0, 0.0, 0.0 }};
    r3_t qq = (r3_t){{ 0.0, 1.0, 0.0 }};
    r3_t rr = (r3_t){{ 0.0, 0.0, 1.0 }};

    h3tt_check_pmap_r3_point("o", &oo, &M, FALSE, &o, "hr3_pmap_aff_from_four_points failed");
    h3tt_check_pmap_r3_point("p", &pp, &M, FALSE, &p, "hr3_pmap_aff_from_four_points failed");
    h3tt_check_pmap_r3_point("q", &qq, &M, FALSE, &q, "hr3_pmap_aff_from_four_points failed");
    h3tt_check_pmap_r3_point("r", &rr, &M, FALSE, &r, "hr3_pmap_aff_from_four_points failed");
  }

void test_hr3_pmap_point(bool_t inv, bool_t verbose)
  { char *tag = (inv ? "_inv" : "");
    if (verbose) { fprintf(stderr, "--- hr3_pmap%s_point ---\n", tag); }
    hr3_point_t p = hr3_point_throw();
    hr3_pmap_t M; h3tt_throw_pmap(&M);

    /* Compute expected images of {p} by {M.dir} and {M.inv}: */
    hr3_point_t q_exp;
    r4x4_map_row(&(p.c), (inv ? &(M.inv) : &(M.dir)), &(q_exp.c));

    /* Compare with library func results: */
    h3tt_check_pmap_point("p", &p, FALSE, &M, inv, &q_exp, txtcat3("hr3_pmap", tag,  "_point failed"));
  }
    
void test_hr3_pmap_plane(bool_t inv, bool_t verbose)
  { char *tag = (inv ? "_inv" : "");
    if (verbose) { fprintf(stderr, "--- hr3_pmap%s_plane ---\n", tag); }
    hr3_plane_t A = hr3_plane_throw();
    hr3_pmap_t M; h3tt_throw_pmap(&M);
    
    /* Compute expected images of {p} by {M.dir} and {M.inv}: */
    hr3_plane_t B_exp;
    r4x4_map_col((inv ? &(M.dir) : &(M.inv)), &(A.f), &(B_exp.f));

    /* Compare with library func results: */
    h3tt_check_pmap_plane("A", &A, FALSE, &M, inv, &B_exp, txtcat3("hr3_pmap", tag, "_plane failed"));
    
    /* Pick a point on the line: */
    hr3_point_t p;
    (void)r4_pick_ortho (&(A.f), &(p.c));
    hr3_point_t q = (inv ? hr3_pmap_inv_point(&p, &M) : hr3_pmap_point(&p, &M));
    double dot = r4_dot(&(B_exp.f), &(q.c));
    affirm(fabs(dot) < 1.0e-11, txtcat3(txtcat3("hr3_pmap", tag, "_plane"), " incompat with ", txtcat3("hr3_pmap", tag, "_point")));
  }

void test_hr3_pmap_is_identity(bool_t verbose)
  { fprintf(stderr, "!! {hr3_pmap_is_identity} NOT TESTED\n"); }

void test_hr3_pmap_identity(bool_t verbose)
  { fprintf(stderr, "!! {hr3_pmap_identity} NOT TESTED\n"); }

void test_hr3_pmap_print(bool_t verbose)
  { fprintf(stderr, "!! {hr3_pmap_print} NOT TESTED\n"); }

void test_hr3_pmap_gen_print(bool_t verbose)
  { fprintf(stderr, "!! {hr3_pmap_gen_print} NOT TESTED\n"); }

void test_hr3_pmap_is_affine(bool_t verbose)
  { fprintf(stderr, "!! {test_hr3_pmap_is_affine} NOT TESTED\n"); }

