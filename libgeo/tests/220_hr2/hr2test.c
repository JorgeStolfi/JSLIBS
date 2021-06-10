/* hr2test --- test program for hr2.h  */
/* Last edited on 2021-06-09 19:54:25 by jstolfi */

#include <hr2.h>

#include <r2.h>
#include <r2x2.h>
#include <r3.h>
#include <r3x3.h>
#include <rn.h>

#include <affirm.h>
#include <jsrandom.h>
#include <flt.h>

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */
 
void do_check_eq(double x, double y, char *msg, char *file, int32_t lnum, const char *func);
  /* If {x} and {y} differ, prints them, prints {msg}, and stops. */

#define check_eq(x,y,msg) \
  do_check_eq((x), (y), (msg), __FILE__, __LINE__, __FUNCTION__)
  
void do_check_eps(double x, double y, double eps, char *msg, char *file, int32_t lnum, const char *func);
  /* If {x} and {y} differ by more than {eps}, prints them, prints {msg}, and stops. */

#define check_eps(x, y, eps, msg) \
  do_check_eps((x), (y), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void do_check_hr2_eps
  ( hr2_point_t *a, 
    hr2_point_t *x, 
    hr2_point_t *y, 
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If points {x} and {y} differ by more than {eps} 
    (apart from homogeneous scaling), prints {a} (if not NULL), 
    prints {x}, {y}, prints {msg}, and stops. */

#define check_hr2_eps(a, x, y, eps, msg)                                 \
  do_check_hr2_eps((a), (x), (y), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void do_check_eq(double x, double y, char *msg, char *file, int32_t lnum, const char *func)
  { if (x != y)
      { fprintf(stderr, " ** %+20.16e %+20.16e differ\n", x, y);
        programerror(msg, file, lnum, func);
      }
  }

void do_check_eps(double x, double y, double eps, char *msg, char *file, int32_t lnum, const char *func)
  { double diff = fabs(x - y);
    if (diff > eps)
      { fprintf(stderr, " ** %+20.16e %+20.16e", x, y);
        fprintf(stderr, " diff = %+20.16e  max = %+20.16e\n", diff, eps);
        programerror(msg, file, lnum, func);
      }
  }

void do_check_hr2_eps
  ( hr2_point_t *a, 
    hr2_point_t *x, 
    hr2_point_t *y, 
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  )
  {
    double diff = hr2_pt_pt_diff(x, y);
    if (diff > eps)
      { fprintf(stderr, " ** ");
        if (a != NULL)
          { rn_gen_print(stderr, NH, a->c.c, "%+2.0f", "[ ", " ", " ]");
            fprintf(stderr, " -> "); 
          }
        r3_t xx;
        r3_dir(&(x->c), &xx);
        rn_gen_print(stderr, NH, xx.c, "%+11.8f", "[ ", " ", " ]");
        fprintf(stderr, " should be "); 
        r3_t yy;
        r3_dir(&(y->c), &yy);
        rn_gen_print(stderr, NH, yy.c, "%+11.8f", "[ ", " ", " ]");
        fprintf(stderr, " sdiff = %20.16e  max = %+20.16e\n", diff, eps);
        programerror(msg, file, lnum, func);
      }
  }

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_hr2(int32_t verbose);
void test_hr2_pmap(int32_t verbose);
void throw_pmap(hr2_pmap_t *m);
double frac (double x);

void check_pmap
  ( char *name, 
    double w, 
    double x, 
    double y, 
    hr2_pmap_t *M, 
    hr2_point_t *q, 
    bool_t flip,
    char *msg
  );
  /* If {flip} is FALSE, maps the point {[w,x,y]} by {M} and compares it
    with {q}. If {flip} is true, tries reversing the signs of {w,x,y} in
    all combinations, and proceeds as above. In either case, returns
    silently if some attempt produced a match (modulo rounding errors);
    aborts with error {msg} if all attempts failed. */

int32_t main (int32_t argc, char **argv)
  {
    int32_t i;
    srand(1993);
    srandom(1993);

    for (i = 0; i < 100; i++) test_hr2(i < 3);
    for (i = 0; i < 100; i++) test_hr2_pmap(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

double frac (double x)
  { int32_t i = (int32_t)x;
    double f = x - (double)i;
    return (f);
  }

void test_hr2(int32_t verbose)
  {
    hr2_point_t p, q, r;
    hr2_line_t L, M, N;
    r2_t pc;
    int32_t i;

    if (verbose)
      { fprintf(stderr,
          "sizeof(hr2_point_t) = %lud  %d*sizeof(double) = %lud\n",
          sizeof(hr2_point_t), NH, NH*sizeof(double)
        );
      }

    if (verbose) { fprintf(stderr, "--- hr2_from_r2 ---\n"); }
    r2_throw_cube(&pc);
    p = hr2_from_r2(&pc);
    affirm(p.c.c[0] == 1.0, "hr2_from_r2 error(1)");
    
    if (verbose) { fprintf(stderr, "--- r2_from_hr2 ---\n"); }
    r3_throw_cube(&(p.c));
    pc = r2_from_hr2(&p);
    q = hr2_from_r2(&pc);
    { double tol = fabs(p.c.c[0])*1.0e-12;
      for (i = 1; i <= NC; i++)
        { double di = q.c.c[i]*p.c.c[0] - p.c.c[i];
          affirm(fabs(di) < tol, "r2_from_hr2 error(1)");
        }
    }
    
    if (verbose) { fprintf(stderr, "--- hr2_pt_pt_diff ---\n"); }
    { /* Check zero distance: */
      r3_throw_cube(&(p.c));
      double dpp = hr2_pt_pt_diff(&p, &p); 
      check_eq(dpp, 0.0, "hr2_pt_pt_diff(p,p) error(1)");
      
      /* Check symmetry: */
      r3_throw_cube(&(p.c));
      r3_throw_cube(&(q.c));
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
      check_eps(dob, dex, 1.0e-8, "hr2_pt_pt_diff error(3)");
      /* Check invariance under rotations: */
      for (i = 0; i < NH; i++)
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
          check_eps(drt, dex, 1.0e-8, "hr2_pt_pt_diff error(4)");
        }
    }
    
    if (verbose) { fprintf(stderr, "--- hr2_side ---\n"); }
    r3_throw_cube(&(p.c));
    r3_throw_cube(&(L.f));
    { double dd = r3_dot(&(p.c), &(L.f));
      sign_t sgn = hr2_side(&p, &L);
      affirm(((dd == 0.0) && (sgn == 0)) || (dd*sgn > 0), "hr2_side error(1)");
    }

    if (verbose) { fprintf(stderr, "--- hr2_orient ---\n"); }
    r3_throw_cube(&(p.c));
    r3_throw_cube(&(q.c));
    r3_throw_cube(&(r.c));
    { double dd = r3_det(&(p.c), &(q.c), &(r.c));
      sign_t sgn1 = hr2_orient(&p, &q, &r);
      affirm(((dd == 0.0) && (sgn1 == 0)) || (dd*sgn1 > 0), "hr2_orient error(1)");
      sign_t sgn2 = hr2_orient(&p, &p, &q);
      affirm(sgn2 == 0, "hr2_orient error(2)");
    }

    if (verbose) { fprintf(stderr, "--- hr2_join ---\n"); }
    r3_throw_cube(&(p.c));
    r3_throw_cube(&(q.c));
    r3_throw_cube(&(r.c)); /* Random test point */
    L = hr2_join(&p, &q);
    { double tp = r3_norm(&(L.f))*r3_norm(&(p.c))*1.0e-12;
      double dp = r3_dot(&(p.c), &(L.f));
      affirm(fabs(dp) < tp, "hr2_join error(1)");
      
      double tq = r3_norm(&(L.f))*r3_norm(&(q.c))*1.0e-12;
      double dq = r3_dot(&(q.c), &(L.f));
      affirm(fabs(dq) < tq, "hr2_join error(2)");
      
      sign_t sgnL = hr2_side(&r, &L);
      sign_t sgnO = hr2_orient(&p, &q, &r);
      affirm(sgnL == sgnO, "hr2_join error(3)");
    }

    if (verbose) { fprintf(stderr, "--- hr2_meet ---\n"); }
    r3_throw_cube(&(L.f)); p.c = L.f;
    r3_throw_cube(&(M.f)); q.c = M.f;
    r = hr2_meet(&L, &M);
    N = hr2_join(&p, &q);
    for (i = 0; i < NH; i++)
      { check_eq(r.c.c[i], N.f.c[i], "hr2_meet error(1)"); }
      
    if (verbose) { fprintf(stderr, "--- hr2_point_point_dir ---\n"); }
    r3_throw_cube(&(p.c)); p.c.c[0] = fabs(p.c.c[0]);
    r3_throw_cube(&(q.c)); q.c.c[0] = fabs(q.c.c[0]);
    { r2_t pc = r2_from_hr2(&p);
      r2_t qc = r2_from_hr2(&q);
      r2_t upq = hr2_point_point_dir(&p, &q);
      r2_t vpq;
      r2_sub(&qc, &pc, &vpq);
      r2_dir(&vpq, &vpq);
      for (i = 0; i < NC; i++)
        { check_eps(upq.c[i], vpq.c[i], 1.0e-12, "hr2_point_point_dir error"); }
    }

    if (verbose) { fprintf(stderr, "--- hr2_line_dir ---\n"); }
    r3_throw_cube(&(p.c)); p.c.c[0] = fabs(p.c.c[0]) + 0.00001;
    r3_throw_cube(&(q.c)); q.c.c[0] = fabs(q.c.c[0]) + 0.00001;
    L = hr2_join(&p, &q);
    { r2_t dL = hr2_line_dir(&L);
      r2_t eL = hr2_point_point_dir(&p, &q);;
      double tol = 1.0e-12;
      for (i = 0; i < NC; i++)
        { check_eps(dL.c[i], eL.c[i], tol, "hr2_line_dir error"); }
    }
      
    if (verbose) { fprintf(stderr, "--- hr2_line_normal ---\n"); }
    r3_throw_cube(&(L.f));
    { r2_t nL = hr2_line_normal(&L);
      r2_t mL = (r2_t){{L.f.c[1], L.f.c[2]}};
      r2_dir(&mL, &mL);
      double tol = 1.0e-12;
      for (i = 0; i < NC; i++)
        { check_eps(nL.c[i], mL.c[i], tol, "hr2_line_normal error"); }
    }

    //  if (verbose)
    //    { 
    //      fprintf(stderr, "!! hr2_ NOT TESTED\n");
    //    }
  }

void test_hr2_pmap(int32_t verbose)
  {
    hr2_pmap_t A/* , B, C */;
    hr2_point_t p, q, r, u;
    /* int32_t i, j, k; */

    /* Size: */
    if (verbose)
      { fprintf(stderr,
          "sizeof(hr2_pmap_t) = %lud  2*%d*%d*sizeof(double) = %lud\n",
          sizeof(hr2_pmap_t), NH, NH, 2*NH*NH*sizeof(double)
        );
      }

    if (verbose) { fprintf(stderr, "--- hr2_pmap_from_mat_disp ---\n"); }
    r2_t disp;
    r2_throw_cube(&disp);
    r2x2_t mat;
    for (int32_t i = 0; i < 2; i++)
      { r2_t s; r2_throw_cube(&s);
        mat.c[i][0] = s.c[0];
        mat.c[i][1] = s.c[1];
      }
    A = hr2_pmap_from_mat_and_disp(&mat, &disp);
    for (int32_t k = 0; k < 5; k++)
      { r2_t p2; r2_throw_cube(&p2);
        r2_t q2; r2x2_map_row(&p2, &mat, &q2);
        r2_add(&disp, &q2, &q2);
        q = hr2_from_r2(&q2);
        double s = dabrandom(0.1, 10.0);
        r3_scale(s, &(q.c), &(q.c));
        check_pmap("p", 1.0, p2.c[0], p2.c[1], &A, &q, FALSE, "hr2_pmap_from_mat_disp failed");
      }
        
    if (verbose) { fprintf(stderr, "--- hr2_pmap_from_points ---\n"); }
    r3_throw_cube(&(p.c));
    r3_throw_cube(&(q.c));
    r3_throw_cube(&(r.c));
    r3_throw_cube(&(u.c));
    A = hr2_pmap_from_points(&p, &q, &r, &u);
    /* Check whether the map works: */
    
    check_pmap("p", 1.0, 0.0, 0.0, &A, &p, FALSE, "hr2_pmap_from_points failed");
    check_pmap("q", 0.0, 1.0, 0.0, &A, &q, FALSE, "hr2_pmap_from_points failed");
    check_pmap("r", 0.0, 0.0, 1.0, &A, &r, FALSE, "hr2_pmap_from_points failed");
    check_pmap("u", 1.0, 1.0, 1.0, &A, &u, TRUE,  "hr2_pmap_from_points failed");
    
    /* TO BE COMPLETED !!! */
    
    // if (verbose) { fprintf(stderr, "--- r2x2_map_row, r2x2_map_col ---\n"); }
    // throw_matrix(&A);
    // r2_throw_cube(&a);
    // r2x2_map_row(&a, &A, &b);
    // r2x2_map_col(&A, &a, &c);
    // r2_zero(&bb);
    // r2_zero(&cc);
    // for (i = 0; i < N; i++)
    //   { for (j = 0; j < N; j++)
    //       { bb.c[j] += a.c[i] * A.c[i][j];
    //         cc.c[i] += A.c[i][j] * a.c[j];
    //       }
    //   }
    // r = r2_dist(&b, &bb);
    // affirm(r < 0.000000001 * r2_norm(&bb), "r2_map_row error");
    // s = r2_dist(&c, &cc);
    // affirm(s < 0.000000001 * r2_norm(&cc), "r2_map_col error");
    // 
    // if (verbose) { fprintf(stderr, "--- r2x2_mul ---\n"); }
    // throw_matrix(&A);
    // throw_matrix(&B);
    // r2x2_mul(&A, &B, &C);
    // for (i = 0; i < N; i++)
    //   { for (j = 0; j < N; j++)
    //       { double sum = 0.0;
    //         for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[k][j]; }
    //         check_eps(C.c[i][j], sum, 0.000000001 * fabs(sum),
    //           "r2x2_mul error"
    //         );
    //       }
    //   }
    // 
    // if (verbose) { fprintf(stderr, "--- r2x2_det ---\n"); }
    // throw_matrix(&A);
    // for (i = 0; i < N; i++)
    //   { int32_t k = (i + 1) % N;
    //     for (j = 0; j < N; j++)
    //       { /* Check for linearity */
    //         r = drandom();
    //         A.c[i][j] = r;
    //         rr = r2x2_det(&A);
    // 
    //         s = drandom();
    //         A.c[i][j] = s;
    //         ss = r2x2_det(&A);
    // 
    //         t = drandom();
    //         A.c[i][j] = r*(1-t) + s*t;
    //         tt = r2x2_det(&A);
    //         mag = fabs(rr) + fabs(ss) + fabs(tt);
    //         check_eps(tt, rr*(1.0 - t) + ss*t, 000000001 * mag,
    //           "r2x2_det error(1)"
    //         );
    //       }
    // 
    //     /* Row swap test: */
    //     r = r2x2_det(&A);
    //     for (j = 0; j < N; j++)
    //       { double t = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = t; }
    //     rr = r2x2_det(&A);
    //     mag = fabs(r) + fabs(rr);
    //     check_eps(r, -rr, 000000001 * mag, "r2x2_det error(2)");
    // 
    //     /* Col swap test: */
    //     r = r2x2_det(&A);
    //     for (j = 0; j < N; j++)
    //       { double t = A.c[j][i]; A.c[j][i] = A.c[j][k]; A.c[j][k] = t; }
    //     rr = r2x2_det(&A);
    //     mag = fabs(r) + fabs(rr);
    //     check_eps(r, -rr, 000000001 * mag, "r2x2_det error(3)");
    //   }
    // 
    // if (verbose) { fprintf(stderr, "--- r2x2_inv ---\n"); }
    // throw_matrix(&A);
    // r2x2_inv(&A, &B);
    // r2x2_mul(&A, &B, &C);
    // for (i = 0; i < N; i++)
    //   { for (j = 0; j < N; j++)
    //       { double val = (i == j ? 1.0 : 0.0);
    //         affirm((C.c[i][j] - val) < 000000001, "r2x2_inv error");
    //       }
    //   }
    // 
    // if (verbose) { fprintf(stderr, "--- r2x2_print ---\n"); }
    // if (verbose)
    //   { throw_matrix (&A);
    //     fprintf(stderr, "A = ");
    //     r2x2_print(stderr, &A);
    //     fputc('\n', stderr);
    //   }

    if (verbose)
      { 
        fprintf(stderr, "!! hr2_pmap_is_identity NOT TESTED\n");
        fprintf(stderr, "!! hr2_pmap_point NOT TESTED\n");
        fprintf(stderr, "!! hr2_pmap_line NOT TESTED\n");
        fprintf(stderr, "!! hr2_pmap_translation NOT TESTED\n");
        fprintf(stderr, "!! hr2_pmap_rotation NOT TESTED\n");
        fprintf(stderr, "!! hr2_pmap_comp NOT TESTED\n");
        fprintf(stderr, "!! hr2_pmap_inv NOT TESTED\n");
        fprintf(stderr, "!! hr2_pmap_from_point_pairs NOT TESTED\n");
      }
  }  

void throw_pmap(hr2_pmap_t *m)
  {
    int32_t i, j;
    r3_t a;
    for (i = 0; i < NH; i++)
      { r3_throw_cube(&a);
        for (j = 0; j < NH; j++) { m->dir.c[i][j] = a.c[j]; }
      }
    r3x3_inv(&(m->dir), &(m->inv));
  }

void check_pmap
  ( char *name, 
    double w, 
    double x, 
    double y, 
    hr2_pmap_t *M, 
    hr2_point_t *q, 
    bool_t flip,
    char *msg
  )
  { 
    hr2_point_t p, pM; 
    if (! flip)
      { /* {w,x,y} must map to {p} as they are: */
        p = (hr2_point_t){{{w, x, y}}};
        pM = hr2_pmap_point(&p, M); 
      }
    else
      { /* Try flipping the signs of {w,x,y} in all combinations: */
        double dmin = +INF;
        int32_t sw, sx, sy;
        for (sy = -1; sy <= +1; sy += 2)
          for (sx = -1; sx <= +1; sx += 2)
            for (sw = -1; sw <= +1; sw += 2)
              { hr2_point_t u = (hr2_point_t){{{sw*w, sx*x, sy*y}}};
                hr2_point_t uM = hr2_pmap_point(&u, M); 
                double d = hr2_pt_pt_diff(&uM, q);
                if (d < dmin) { dmin = d; p = u; pM = uM; }
              }
      }
    check_hr2_eps(&p, &pM, q, 0.00000001, msg);

  }
