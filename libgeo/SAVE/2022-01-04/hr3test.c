/* hr3test --- test program for hr3.h  */
/* Last edited on 2022-01-04 08:42:35 by stolfi */

#include <hr3.h>

#include <r3.h>
#include <r3x3.h>
#include <r4.h>
#include <r4x4.h>
#include <r6.h>

#include <affirm.h>
#include <jsrandom.h>
#include <flt.h>

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define NH 4
  /* Number of homogeneous coordinates in a point. */

#define NC 3
  /* Number of Cartesian coordinates in a point. */
 
#define NP 6
  /* Number of Pluecker coordinates in a line. */
 
#define check_eq(x, y, msg) \
  if ((x) != (y)) \
    { fprintf(stderr, " ** %20.16e %20.16e\n", (x), (y)); \
      affirm(0, msg); \
    } \
  else { }

#define check_eps(x, y, eps, msg) \
  if (fabs((x) - (y)) > (eps)) \
    { fprintf(stderr, " ** %20.16e %20.16e %20.16e\n", (x), (y), (eps)); \
      affirm(0, msg); \
    } \
  else { }

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_hr3(int32_t verbose);
void test_hr3_pmap(int32_t verbose);
void throw_pmap(hr3_pmap_t *m);
double frac (double x);

int32_t main (int32_t argc, char **argv)
  {
    int32_t i;
    srand(1993);
    srandom(1993);

    for (i = 0; i < 100; i++) test_hr3(i < 3);
    for (i = 0; i < 100; i++) test_hr3_pmap(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

double frac (double x)
  { int32_t i = (int32_t)x;
    double f = x - (double)i;
    return (f);
  }

void test_hr3(int32_t verbose)
  {
    hr3_point_t p, q, r, s;
    hr3_line_t F, G;
    hr3_plane_t L, M, N, O;
    r3_t pc;
    int32_t i;

    if (verbose)
      { fprintf(stderr,
          "sizeof(hr3_point_t) = %lu  %d*sizeof(double) = %lu\n",
          sizeof(hr3_point_t), NH, NH*sizeof(double)
        );
      }

    if (verbose) { fprintf(stderr, "--- hr3_from_r3 ---\n"); }
    r3_throw_cube(&pc);
    p = hr3_from_r3(&pc);
    affirm(p.c.c[0] == 1.0, "hr3_from_r3 error(1)");
    
    if (verbose) { fprintf(stderr, "--- r3_from_hr3 ---\n"); }
    r4_throw_cube(&(p.c));
    pc = r3_from_hr3(&p);
    q = hr3_from_r3(&pc);
    { double tol = fabs(p.c.c[0])*1.0e-12;
      for (i = 1; i < NC; i++)
        { double di = q.c.c[i]*p.c.c[0] - p.c.c[i];
          affirm(fabs(di) < tol, "r3_from_hr3 error(1)");
        }
    }
    
    if (verbose) { fprintf(stderr, "--- hr3_side ---\n"); }
    r4_throw_cube(&(p.c));
    r4_throw_cube(&(L.f));
    { double dd = r4_dot(&(p.c), &(L.f));
      sign_t sgn = hr3_side(&p, &L);
      affirm(((dd == 0.0) && (sgn == 0)) || (dd*sgn > 0), "hr3_side error(1)");
    }

    if (verbose) { fprintf(stderr, "--- hr3_orient ---\n"); }
    p = (hr3_point_t){{{1.0, 0.0, 0.0, 0.0}}};
    q = (hr3_point_t){{{0.0, 1.0, 0.0, 0.0}}};
    r = (hr3_point_t){{{0.0, 0.0, 1.0, 0.0}}};
    s = (hr3_point_t){{{0.0, 0.0, 0.0, 1.0}}};
    { sign_t sgn0 = hr3_orient(&p, &q, &r, &s);
      affirm(sgn0 == +1, "hr3_orient error(0.0)");
      
      sign_t sgn1 = hr3_orient(&q, &p, &r, &s);
      affirm(sgn1 == -1, "hr3_orient error(0.1)");
      
      /* Apply some positive transformation: */
      for (i = 0; i < NH; i++)
        { double t = s.c.c[i] + 0.500*p.c.c[i];
          p.c.c[i] = p.c.c[i] + 0.250*q.c.c[i];
          q.c.c[i] = q.c.c[i] + 0.125*r.c.c[i];
          r.c.c[i] = r.c.c[i] + 0.750*s.c.c[i];
          s.c.c[i] = t;
        }
        
      sign_t sgn2 = hr3_orient(&p, &q, &r, &s);
      affirm(sgn2 == +1, "hr3_orient error(0.2)");
    }
    /* Check with random points: */
    r4_throw_cube(&(p.c));
    r4_throw_cube(&(q.c));
    r4_throw_cube(&(r.c));
    r4_throw_cube(&(s.c));
    { double dd = r4_det(&(p.c), &(q.c), &(r.c), &(s.c));
      sign_t sgn1 = hr3_orient(&p, &q, &r, &s);
      affirm(((dd == 0.0) && (sgn1 == 0)) || (dd*sgn1 > 0), "hr3_orient error(1)");
      sign_t sgn2 = hr3_orient(&p, &p, &q, &r);
      affirm(sgn2 == 0, "hr3_orient error(2)");
    }

    if (verbose) { fprintf(stderr, "--- hr3_line_from_two_points ---\n"); }
    p = (hr3_point_t){{{1.0, 0.0, 0.0, 0.0}}};
    q = (hr3_point_t){{{0.0, 1.0, 0.0, 0.0}}};
    F = hr3_line_from_two_points(&p, &q);
    { double tol = 1.0e-12;
      hr3_line_t G = (hr3_line_t){{{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}}};
      for (i = 0; i < NP; i++)
        { check_eps(F.k.c[i], G.k.c[i], tol, "hr3_line_from_two_points error(1.0)"); }
    }
    /* Check with random points: */
    r4_throw_cube(&(p.c));
    r4_throw_cube(&(q.c));
    F = hr3_line_from_two_points(&p, &q);
    { double tol = r6_norm(&(F.k))*r4_norm(&(p.c))*1.0e-12;
      hr3_plane_t L = hr3_plane_from_line_and_point(&F, &p);
      for (i = 0; i < NH; i++)
        { check_eps(L.f.c[i], 0.0, tol, "hr3_line_from_two_points error(1.1)"); }
    }
    {
      double tol = r6_norm(&(F.k))*r4_norm(&(q.c))*1.0e-12;
      hr3_plane_t L = hr3_plane_from_line_and_point(&F, &q);
      for (i = 0; i < NH; i++)
        { check_eps(L.f.c[i], 0.0, tol, "hr3_line_from_two_points error(1.2)"); }
    }
      
    if (verbose) { fprintf(stderr, "--- hr3_line_dir ---\n"); }
    r4_throw_cube(&(p.c)); p.c.c[0] = fabs(p.c.c[0]) + 0.001;
    r4_throw_cube(&(q.c)); q.c.c[0] = fabs(q.c.c[0]) + 0.001;
    F = hr3_line_from_two_points(&p, &q);
    { r3_t upq = hr3_line_dir(&F);
      r3_t vpq = hr3_point_point_dir(&p, &q);
      for (i = 0; i < NC; i++)
        { check_eps(upq.c[i], vpq.c[i], 1.0e-12, "hr3_line_dir error"); }
    }

    if (verbose) { fprintf(stderr, "--- hr3_plane_from_three_points ---\n"); }
    p = (hr3_point_t){{{1.0, 0.0, 0.0, 0.0}}};
    q = (hr3_point_t){{{0.0, 1.0, 0.0, 0.0}}};
    r = (hr3_point_t){{{0.0, 0.0, 1.0, 0.0}}};
    { L = hr3_plane_from_three_points(&p, &q, &r);
      M = (hr3_plane_t){{{0.0, 0.0, 0.0, 1.0}}};
      for (i = 0; i < NH; i++)
        { check_eq(L.f.c[i], M.f.c[i], "hr3_plane_from_three_points error(1.1)"); }
    }
    { L = hr3_plane_from_three_points(&q, &p, &r);
      M = (hr3_plane_t){{{0.0, 0.0, 0.0, -1.0}}};
      for (i = 0; i < NH; i++)
        { check_eq(L.f.c[i], M.f.c[i], "hr3_plane_from_three_points error(1.2)"); }
    }
    /* Check with random points: */
    r4_throw_cube(&(p.c));
    r4_throw_cube(&(q.c));
    r4_throw_cube(&(r.c));
    r4_throw_cube(&(s.c)); /* Random test point */
    L = hr3_plane_from_three_points(&p, &q, &r);
    { double tp = r4_norm(&(L.f))*r4_norm(&(p.c))*1.0e-12;
      double dp = r4_dot(&(p.c), &(L.f));
      affirm(fabs(dp) < tp, "hr3_plane_from_three_points error(1)");
      
      double tq = r4_norm(&(L.f))*r4_norm(&(q.c))*1.0e-12;
      double dq = r4_dot(&(q.c), &(L.f));
      affirm(fabs(dq) < tq, "hr3_plane_from_three_points error(2)");
      
      double tr = r4_norm(&(L.f))*r4_norm(&(r.c))*1.0e-12;
      double dr = r4_dot(&(r.c), &(L.f));
      affirm(fabs(dr) < tr, "hr3_plane_from_three_points error(3)");
      
      sign_t sgnL = hr3_side(&s, &L);
      sign_t sgnO = hr3_orient(&p, &q, &r, &s);
      affirm(sgnL == sgnO, "hr3_plane_from_three_points error(4)");
    }
      
    if (verbose) { fprintf(stderr, "--- hr3_point_point_dir ---\n"); }
    r4_throw_cube(&(p.c)); p.c.c[0] = fabs(p.c.c[0]);
    r4_throw_cube(&(q.c)); q.c.c[0] = fabs(q.c.c[0]);
    { r3_t pc = r3_from_hr3(&p);
      r3_t qc = r3_from_hr3(&q);
      r3_t upq = hr3_point_point_dir(&p, &q);
      r3_t vpq;
      r3_sub(&qc, &pc, &vpq);
      r3_dir(&vpq, &vpq);
      for (i = 0; i < NC; i++)
        { check_eps(upq.c[i], vpq.c[i], 1.0e-12, "hr3_point_point_dir error"); }
    }

    if (verbose) { fprintf(stderr, "--- hr3_plane_normal ---\n"); }
    r4_throw_cube(&(L.f));
    { r3_t nL = hr3_plane_normal(&L);
      r3_t mL = (r3_t){{L.f.c[1], L.f.c[2], L.f.c[3]}};
      r3_dir(&mL, &mL);
      double tol = 1.0e-12;
      for (i = 0; i < NC; i++)
        { check_eps(nL.c[i], mL.c[i], tol, "hr3_plane_normal error"); }
    }
    
    if (verbose) { fprintf(stderr, "--- hr3_point_from_three_planes ---\n"); }
    L = (hr3_plane_t){{{0.0, 1.0, 0.0, 0.0}}};
    M = (hr3_plane_t){{{0.0, 0.0, 1.0, 0.0}}};
    N = (hr3_plane_t){{{0.0, 0.0, 0.0, 1.0}}};
    s = hr3_point_from_three_planes(&L, &M, &N);
    { hr3_point_t t = (hr3_point_t){{{1.0, 0.0, 0.0, 0.0}}};
      for (i = 0; i < NH; i++)
        { check_eq(s.c.c[i], t.c.c[i], "hr3_point_from_three_planes error(0)"); }
    }
    /* Check with random planes: */
    r4_throw_cube(&(L.f)); p.c = L.f;
    r4_throw_cube(&(M.f)); q.c = M.f;
    r4_throw_cube(&(N.f)); r.c = N.f;
    s = hr3_point_from_three_planes(&L, &M, &N);
    O = hr3_plane_from_three_points(&r, &q, &p);
    for (i = 0; i < NH; i++)
      { check_eq(s.c.c[i], O.f.c[i], "hr3_point_from_three_planes error(1)"); }
      
    if (verbose) { fprintf(stderr, "--- hr3_plane_from_line_and_point ---\n"); }
    p = (hr3_point_t){{{1.0, 0.0, 0.0, 0.0}}};
    q = (hr3_point_t){{{0.0, 1.0, 0.0, 0.0}}};
    r = (hr3_point_t){{{0.0, 0.0, 1.0, 0.0}}};
    F = hr3_line_from_two_points(&p, &q);
    L = hr3_plane_from_line_and_point(&F, &r);
    M = hr3_plane_from_three_points(&p, &q, &r);
    { double sL = r4_norm(&(L.f));
      double sM = r4_norm(&(M.f));
      double tol = sL*sM*1.0e-12;
      for (i = 0; i < NH; i++)
        { check_eps(sM*L.f.c[i], sL*M.f.c[i], tol, "hr3_plane_from_line_and_point error(1)"); }
    }
   
    if (verbose) { fprintf(stderr, "--- hr3_line_from_two_planes ---\n"); }
    L = (hr3_plane_t){{{0.0, 1.0, 0.0, 0.0}}};
    M = (hr3_plane_t){{{0.0, 0.0, 1.0, 0.0}}};
    F = hr3_line_from_two_planes(&L, &M);
    p = (hr3_point_t){{{1.0, 0.0, 0.0, 0.0}}};
    q = (hr3_point_t){{{0.0, 0.0, 0.0, 1.0}}};
    G = hr3_line_from_two_points(&p, &q);
    { double sF = r6_norm(&(F.k));
      double sG = r6_norm(&(G.k));
      double tol = sF*sG*1.0e-12;
      for (i = 0; i < NP; i++)
        { check_eps(sG*F.k.c[i], sF*G.k.c[i], tol, "hr3_line_from_two_planes error(1)"); }
    }
    /* Check with random points: */
    r4_throw_cube(&(p.c));
    r4_throw_cube(&(q.c));
    r4_throw_cube(&(r.c));
    r4_throw_cube(&(s.c));
    L = hr3_plane_from_three_points(&p, &q, &r);
    M = hr3_plane_from_three_points(&q, &r, &s);
    F = hr3_line_from_two_planes(&L, &M);
    G = hr3_line_from_two_points(&q, &r);
    { sign_t sgn = hr3_orient(&p, &q, &r, &s);
      double sF = r6_norm(&(F.k));
      double sG = r6_norm(&(G.k));
      double tol = sF*sG*1.0e-12;
      for (i = 0; i < NP; i++)
        { check_eps(sG*F.k.c[i], sgn*sF*G.k.c[i], tol, "hr3_line_from_two_planes error(2)"); }
    }

    if (verbose) { fprintf(stderr, "--- hr3_point_from_line_and_plane ---\n"); }
    p = (hr3_point_t){{{1.0, 0.0, 0.0, 0.0}}};
    q = (hr3_point_t){{{0.0, 1.0, 0.0, 0.0}}};
    F = hr3_line_from_two_points(&p, &q);
    L = (hr3_plane_t){{{1.0, 0.0, 0.0, 0.0}}};
    s = hr3_point_from_line_and_plane(&F, &L);
    { double tol = 1.0e-12;
      for (i = 0; i < NH; i++)
        { check_eps(s.c.c[i], -q.c.c[i], tol, "hr3_point_from_line_and_plane error(0)"); }
    }
    /* Check with random planes: */
    r4_throw_cube(&(p.c));
    r4_throw_cube(&(q.c));
    r4_throw_cube(&(r.c));
    r4_throw_cube(&(s.c));
    F = hr3_line_from_two_points(&p, &q);
    L = hr3_plane_from_three_points(&q, &r, &s);
    { sign_t sgn = hr3_orient(&p, &q, &r, &s);
      hr3_point_t t = hr3_point_from_line_and_plane(&F, &L);
      double st = r4_norm(&(t.c));
      double sq = r4_norm(&(q.c));
      double tol = st*sq*1.0e-12;
      for (i = 0; i < NH; i++)
        { check_eps(sq*t.c.c[i], sgn*st*q.c.c[i], tol, "hr3_point_from_line_and_plane error(1)"); }
    }
  }

void test_hr3_pmap(int32_t verbose)
  {
    // hr3_pmap_t A, B, C;
    // hr3_point_t p, q, r;
    // int32_t i, j, k;

    /* Size: */
    if (verbose)
      { fprintf(stderr,
          "sizeof(hr3_pmap_t) = %lu  2*%d*%d*sizeof(double) = %lu\n",
          sizeof(hr3_pmap_t), NH, NH, 2*NH*NH*sizeof(double)
        );
      }

    /* TO BE COMPLETED !!! */
    
    // if (verbose) { fprintf(stderr, "--- r3x3_map_row, r3x3_map_col ---\n"); }
    // throw_matrix(&A);
    // r3_throw_cube(&a);
    // r3x3_map_row(&a, &A, &b);
    // r3x3_map_col(&A, &a, &c);
    // r3_zero(&bb);
    // r3_zero(&cc);
    // for (i = 0; i < N; i++)
    //   { for (j = 0; j < N; j++)
    //       { bb.c[j] += a.c[i] * A.c[i][j];
    //         cc.c[i] += A.c[i][j] * a.c[j];
    //       }
    //   }
    // r = r3_dist(&b, &bb);
    // affirm(r < 0.000000001 * r3_norm(&bb), "r3_map_row error");
    // s = r3_dist(&c, &cc);
    // affirm(s < 0.000000001 * r3_norm(&cc), "r3_map_col error");
    // 
    // if (verbose) { fprintf(stderr, "--- r3x3_mul ---\n"); }
    // throw_matrix(&A);
    // throw_matrix(&B);
    // r3x3_mul(&A, &B, &C);
    // for (i = 0; i < N; i++)
    //   { for (j = 0; j < N; j++)
    //       { double sum = 0.0;
    //         for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[k][j]; }
    //         check_eps(C.c[i][j], sum, 0.000000001 * fabs(sum),
    //           "r3x3_mul error"
    //         );
    //       }
    //   }
    // 
    // if (verbose) { fprintf(stderr, "--- r3x3_det ---\n"); }
    // throw_matrix(&A);
    // for (i = 0; i < N; i++)
    //   { int32_t k = (i + 1) % N;
    //     for (j = 0; j < N; j++)
    //       { /* Check for linearity */
    //         r = drandom();
    //         A.c[i][j] = r;
    //         rr = r3x3_det(&A);
    // 
    //         s = drandom();
    //         A.c[i][j] = s;
    //         ss = r3x3_det(&A);
    // 
    //         t = drandom();
    //         A.c[i][j] = r*(1-t) + s*t;
    //         tt = r3x3_det(&A);
    //         mag = fabs(rr) + fabs(ss) + fabs(tt);
    //         check_eps(tt, rr*(1.0 - t) + ss*t, 000000001 * mag,
    //           "r3x3_det error(1)"
    //         );
    //       }
    // 
    //     /* Row swap test: */
    //     r = r3x3_det(&A);
    //     for (j = 0; j < N; j++)
    //       { double t = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = t; }
    //     rr = r3x3_det(&A);
    //     mag = fabs(r) + fabs(rr);
    //     check_eps(r, -rr, 000000001 * mag, "r3x3_det error(2)");
    // 
    //     /* Col swap test: */
    //     r = r3x3_det(&A);
    //     for (j = 0; j < N; j++)
    //       { double t = A.c[j][i]; A.c[j][i] = A.c[j][k]; A.c[j][k] = t; }
    //     rr = r3x3_det(&A);
    //     mag = fabs(r) + fabs(rr);
    //     check_eps(r, -rr, 000000001 * mag, "r3x3_det error(3)");
    //   }
    // 
    // if (verbose) { fprintf(stderr, "--- r3x3_inv ---\n"); }
    // throw_matrix(&A);
    // r3x3_inv(&A, &B);
    // r3x3_mul(&A, &B, &C);
    // for (i = 0; i < N; i++)
    //   { for (j = 0; j < N; j++)
    //       { double val = (i == j ? 1.0 : 0.0);
    //         affirm((C.c[i][j] - val) < 000000001, "r3x3_inv error");
    //       }
    //   }
    // 
    // if (verbose) { fprintf(stderr, "--- r3x3_print ---\n"); }
    // if (verbose)
    //   { throw_matrix (&A);
    //     fprintf(stderr, "A = ");
    //     r3x3_print(stderr, &A);
    //     fputc('\n', stderr);
    //   }
  }  

void throw_pmap(hr3_pmap_t *m)
  {
    int32_t i, j;
    r4_t a;
    for (i = 0; i < NH; i++)
      { r4_throw_cube(&a);
        for (j = 0; j < NH; j++) { m->dir.c[i][j] = a.c[j]; }
      }
    r4x4_inv(&(m->dir), &(m->inv));
  }
