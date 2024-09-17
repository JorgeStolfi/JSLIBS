/* hi2test --- test program for hi2.h  */
/* Last edited on 2024-09-04 20:50:47 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <i2.h>
/* #include <i2x2.h> */
#include <i3.h>
/* #include <i3x3.h> */

#include <affirm.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <jswsize.h>

#include <hi2.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */
 
void do_check_eq(int64_t x, int64_t y, char *msg, char *file, int32_t lnum, const char *func);
  /* If {x} and {y} differ, prints them, prints {msg}, and stops. */

#define check_eq(x,y,msg) \
  do_check_eq((x), (y), (msg), __FILE__, __LINE__, __FUNCTION__)
 
void do_check_i2_coll(i2_t *u, i2_t *v, char *msg, char *file, int32_t lnum, const char *func);
  /* If {u} and {v} are not collinear with same direction (or both null),
    prints them, prints {msg} and stops. */
 
#define check_i2_coll(u,v,msg) \
  do_check_i2_coll((u), (v), (msg), __FILE__, __LINE__, __FUNCTION__)
 
void do_check_i3_coll(i3_t *u, i3_t *v, char *msg, char *file, int32_t lnum, const char *func);
  /* If {u} and {v} are not collinear with same direction (or both null),
    prints them, prints {msg} and stops. */
 
#define check_i3_coll(u,v,msg) \
  do_check_i3_coll((u), (v), (msg), __FILE__, __LINE__, __FUNCTION__)
 
void do_check_eq(int64_t x, int64_t y, char *msg, char *file, int32_t lnum, const char *func)
  { if (x != y)
      { fprintf(stderr, (" ** %+" int64_d_fmt " %+" int64_d_fmt " differ\n"), x, y);
        programerror(msg, file, lnum, func);
      }
  }

void do_check_i2_coll(i2_t *u, i2_t *v, char *msg, char *file, int32_t lnum, const char *func)
  {
    int32_t ug = (int32_t)gcd(abs(u->c[0]), abs(u->c[1])); 
    int32_t vg = (int32_t)gcd(abs(v->c[0]), abs(v->c[1])); 
    demand((ug == 0) == (vg == 0), msg);
    int32_t i;
    for (i = 0; i < NC; i++) { do_check_eq(u->c[i]/ug, v->c[i]/vg, msg, file, lnum, func); }
  }

void do_check_i3_coll(i3_t *u, i3_t *v, char *msg, char *file, int32_t lnum, const char *func)
  {
    int32_t ug = (int32_t)gcd(gcd(abs(u->c[0]), abs(u->c[1])), abs(u->c[2])); 
    int32_t vg = (int32_t)gcd(gcd(abs(v->c[0]), abs(v->c[1])), abs(u->c[2]));  
    demand((ug == 0) == (vg == 0), msg);
    int32_t i;
    for (i = 0; i < NH; i++) { do_check_eq(u->c[i]/ug, v->c[i]/vg, msg, file, lnum, func); }
  }

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_hi2(int32_t verbose);
/* void test_hi2_pmap(int32_t verbose); */
/* void throw_pmap(hi2_pmap_t *m); */

// void check_pmap
//   ( char *name, 
//     double w, 
//     double x, 
//     double y, 
//     hi2_pmap_t *M, 
//     hi2_point_t *q, 
//     bool_t twirl,
//     char *msg
//   );
//   /* If {twirl} is FALSE, maps the point {[w,x,y]} by {M} and compares it
//     with {q}. If {twirl} is true, tries reversing the signs of {w,x,y} in
//     all combinations, and proceeds as above. In either case, returns
//     silently if some attempt produced a match (modulo rounding errors);
//     aborts with error {msg} if all attempts failed. */

int32_t main (int32_t argc, char **argv)
  {
    int32_t i;
    srand(1993);
    srandom(1993);

    for (i = 0; i < 100; i++) test_hi2(i < 3);
    /* for (i = 0; i < 100; i++) test_hi2_pmap(i < 3); */
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_hi2(int32_t verbose)
  {
    hi2_point_t p, q, r;
    hi2_line_t L, M, N;
    int32_t trad = 4615;
    i2_t pc;
    int32_t i;

    if (verbose)
      { fprintf(stderr,
          "sizeof(hi2_point_t) = %lu  %d*sizeof(int32_t) = %lu\n",
          sizeof(hi2_point_t), NH, NH*sizeof(int32_t)
        );
      }

    if (verbose) { fprintf(stderr, "--- hi2_from_i2 ---\n"); }
    i2_throw_cube(trad, &pc);
    p = hi2_from_i2(&pc);
    affirm(p.c.c[0] == 1, "hi2_from_i2 error(w)");
    affirm(p.c.c[1] == pc.c[0], "hi2_from_i2 error(x)");
    affirm(p.c.c[2] == pc.c[1], "hi2_from_i2 error(y)");
    
    if (verbose) { fprintf(stderr, "--- hi2_side ---\n"); }
    i3_throw_cube(trad, &(p.c));
    i3_throw_cube(trad, &(L.f));
    { int64_t dd = i3_dot(&(p.c), &(L.f));
      sign_t sgn = hi2_side(&p, &L);
      affirm(((dd == 0) && (sgn == 0)) || (dd*sgn > 0), "hi2_side error(1)");
    }

    if (verbose) { fprintf(stderr, "--- hi2_orient ---\n"); }
    i3_throw_cube(trad, &(p.c));
    i3_throw_cube(trad, &(q.c));
    i3_throw_cube(trad, &(r.c));
    { int64_t dd = i3_det(&(p.c), &(q.c), &(r.c));
      sign_t sgn1 = hi2_orient(&p, &q, &r);
      affirm(((dd == 0) && (sgn1 == 0)) || (dd*sgn1 > 0), "hi2_orient error(1)");
      sign_t sgn2 = hi2_orient(&p, &p, &q);
      affirm(sgn2 == 0, "hi2_orient error(2)");
    }

    if (verbose) { fprintf(stderr, "--- hi2_point_point_dir ---\n"); }
    i3_throw_cube(trad, &(p.c)); p.c.c[0] = abs(p.c.c[0]);
    i3_throw_cube(trad, &(q.c)); q.c.c[0] = abs(q.c.c[0]);
    { i2_t upq = hi2_point_point_dir(&p, &q);
      
      /* Move {p,q} to the hither or infinity, count flips in {spq}: */
      sign_t spq = +1;
      if (p.c.c[0] < 0) { spq = -spq; for (i=0; i < NH; i++) { p.c.c[i] = -p.c.c[i]; }}
      if (q.c.c[0] < 0) { spq = -spq; for (i=0; i < NH; i++) { q.c.c[i] = -q.c.c[i]; }}

      /* Adjust {p,q} to the same weight: */
      int32_t pw = p.c.c[0]; assert(pw >= 0);
      int32_t qw = q.c.c[0]; assert(qw >= 0);
      for (i=0; i < NH; i++) { p.c.c[i] *= qw; q.c.c[i] *= pw; }
      
      /* Now get the direction vector: */
      i2_t vpq = (i2_t){{ spq*(q.c.c[1] - p.c.c[1]), spq*(q.c.c[2] - p.c.c[2]) }};
      
      /* Compare the two vectors: */
      check_i2_coll(&upq, &vpq, "hi2_point_point_dir error");
    }

    if (verbose) { fprintf(stderr, "--- hi2_join ---\n"); }
    i3_throw_cube(trad, &(p.c));
    i3_throw_cube(trad, &(q.c));
    i3_throw_cube(trad, &(r.c)); /* Random test point */
    L = hi2_join(&p, &q);
    { affirm(i3_dot(&(p.c), &(L.f)) == 0, "hi2_join error(1)");
      affirm(i3_dot(&(q.c), &(L.f)) == 0, "hi2_join error(2)");
      sign_t sgnL = hi2_side(&r, &L);
      sign_t sgnO = hi2_orient(&p, &q, &r);
      affirm(sgnL == sgnO, "hi2_join error(3)");
    }

    if (verbose) { fprintf(stderr, "--- hi2_meet ---\n"); }
    i3_throw_cube(trad, &(L.f)); p.c = L.f;
    i3_throw_cube(trad, &(M.f)); q.c = M.f;
    r = hi2_meet(&L, &M);
    /* Check by duality: */
    N = hi2_join(&p, &q);
    check_i3_coll(&(r.c), &(N.f), "hi2_meet error(1)");
      
    if (verbose) { fprintf(stderr, "--- hi2_line_dir ---\n"); }
    i3_throw_cube(trad, &(p.c)); p.c.c[0] = abs(p.c.c[0]);
    i3_throw_cube(trad, &(q.c)); q.c.c[0] = abs(q.c.c[0]);
    L = hi2_join(&p, &q);
    { i2_t dL = hi2_line_dir(&L);
      i2_t eL = hi2_point_point_dir(&p, &q);;
      check_i2_coll(&dL, &eL, "hi2_line_dir error");
    }
      
    if (verbose) { fprintf(stderr, "--- hi2_line_normal ---\n"); }
    i3_throw_cube(trad, &(L.f));
    { i2_t nL = hi2_line_normal(&L);
      i2_t mL = (i2_t){{L.f.c[1], L.f.c[2]}};
      check_i2_coll(&nL, &mL, "hi2_line_normal error");
    }

    if (verbose) { fprintf(stderr, "--- hi2_dist_sqr ---\n"); }
    i3_throw_cube(trad, &(p.c)); p.c.c[0] = abs(p.c.c[0]);
    i3_throw_cube(trad, &(q.c)); q.c.c[0] = abs(q.c.c[0]);
    { urat64_t dpq = hi2_dist_sqr(&p, &q);
      
      /* Move {p,q} to the hither or infinity, count flips in {spq}: */
      if (p.c.c[0] < 0) { for (i=0; i < NH; i++) { p.c.c[i] = -p.c.c[i]; }}
      if (q.c.c[0] < 0) { for (i=0; i < NH; i++) { q.c.c[i] = -q.c.c[i]; }}

      /* Adjust {p,q} to the same weight: */
      int32_t pw = p.c.c[0]; assert(pw >= 0);
      int32_t qw = q.c.c[0]; assert(qw >= 0);
      for (i=0; i < NH; i++) { p.c.c[i] *= qw; q.c.c[i] *= pw; }
      assert(p.c.c[0] == q.c.c[0]);
      
      /* Now get the distance squared: */
      int64_t nd = i3_dist_sqr(&(p.c), &(q.c));
      int64_t dd = p.c.c[0]*(int64_t)q.c.c[0];
      
      urat64_t epq = (urat64_t){ nd, dd };
      
      /* Compare the two fractions: */
      demand(urat64_compare(&dpq, &epq) == 0, "hi2_dist_sqr error");
    }

    if (verbose) { fprintf(stderr, "--- hi2_in_circle ---\n"); }
    /* Test on special circle: */
    /* Beware of overflow coords must be in {-M..+M} with {M <= ???}: */
    {
      int32_t cx = int32_abrandom(-100,+100);
      int32_t cy = int32_abrandom(-100,+100);
      int32_t cm = int32_abrandom(1,10);
      
      hi2_point_t a = (hi2_point_t){{{ 1, -2*cm + cx, -1*cm + cy }}}; /* Rectangle corner 1. */
      hi2_point_t b = (hi2_point_t){{{ 1, +2*cm + cx, -1*cm + cy }}}; /* Rectangle corner 2. */
      hi2_point_t c = (hi2_point_t){{{ 1, -2*cm + cx, +1*cm + cy }}}; /* Rectangle corner 3. */
      hi2_point_t d = (hi2_point_t){{{ 1, +2*cm + cx, +1*cm + cy }}}; /* Rectangle corner 4. */
      hi2_point_t e = (hi2_point_t){{{ 1, 00*cm + cx, 00*cm + cy }}}; /* Center of rectangle. */

      hi2_point_t f; /* Probe point. */
      
      if (verbose)
        { /* Plot result: */
          int32_t iw, ix, iy;
          iw = 10;
          for (iy = -30; iy <= 30; iy++)
            { fprintf(stderr, "%4d/%-4d ", iy, iw);
              for (ix = -30; ix <= +30; ix++)
                { f = (hi2_point_t){{{ iw, ix*cm + iw*cx, iy*cm + iw*cy }}};
                  sign_t Sabcf = hi2_in_circle(&a, &b, &c, &f);
                  fprintf(stderr, "%s", (Sabcf == 0 ? "O" : (Sabcf < 0 ? "-" : "+")));
                }
              fprintf(stderr, "\n");
            }
        }
      
      /* Test by distance: */
      for (i = 0; i < 10; i ++)
        { if (i == 0)
            { f = a; }
          else if (i == 1)
            { f = b; }
          else if (i == 2)
            { f = c; }
          else if (i == 3)
            { f = d; }
          else if (i == 4)
            { f = e; }
          else
            { i3_throw_cube(trad, &(f.c)); }
          /* Any three CCW corners of the rectangle must give the same circle: */
          sign_t Sabcf = hi2_in_circle(&a, &b, &c, &f);
          sign_t Sabdf = hi2_in_circle(&a, &b, &d, &f);
          check_eq(Sabcf, Sabdf, "hi2_in_circle error(1)");
          sign_t Sbdaf = hi2_in_circle(&b, &d, &a, &f);
          check_eq(Sabcf, Sbdaf, "hi2_in_circle error(2)");
          
          /* Any three CW corners of the rectangle must give the inverted circle: */
          sign_t Sbacf = hi2_in_circle(&b, &a, &c, &f);
          check_eq(Sabcf, -Sbacf, "hi2_in_circle error(3)");
          
          /* Point {f} is inside circle {a,b,c} iff {a} is inside {b,c,f}: */
          sign_t Sbfca = hi2_in_circle(&b, &f, &c, &a);
          check_eq(Sabcf, Sbfca, "hi2_in_circle error(4)");

          /* Point {f} is inside circle {a,b,c} iff {b} is outside {a,c,f}: */
          sign_t Sacfb = hi2_in_circle(&b, &c, &f, &a);
          check_eq(Sabcf, -Sacfb, "hi2_in_circle error(5)");

          /* Check with distances: */
          urat64_t dae = hi2_dist_sqr(&a, &e);
          urat64_t dfe = hi2_dist_sqr(&f, &e);
          sign_t T = urat64_compare(&dfe, &dae);
          check_eq(Sabcf, -T, "hi2_in_circle error(6)");
        }
    }
  }

// void test_hi2_pmap(int32_t verbose)
//   {
//     hi2_pmap_t A/* , B, C */;
//     hi2_point_t p, q, r, u;
//     /* int32_t i, j, k; */
// 
//     /* Size: */
//     if (verbose)
//       { fprintf(stderr,
//           "sizeof(hi2_pmap_t) = %d  2*%d*%d*sizeof(double) = %d\n",
//           sizeof(hi2_pmap_t), NH, NH, 2*NH*NH*sizeof(double)
//         );
//       }
// 
//     if (verbose) { fprintf(stderr, "--- hi2_pmap_from_points ---\n"); }
//     i3_throw_cube(trad, &(p.c));
//     i3_throw_cube(trad, &(q.c));
//     i3_throw_cube(trad, &(r.c));
//     i3_throw_cube(trad, &(u.c));
//     A = hi2_pmap_from_points(&p, &q, &r, &u);
//     /* Check whether the map works: */
//     
//     check_pmap("p", 1.0, 0.0, 0.0, &A, &p, FALSE, "hi2_pmap_from_points failed");
//     check_pmap("q", 0.0, 1.0, 0.0, &A, &q, FALSE, "hi2_pmap_from_points failed");
//     check_pmap("r", 0.0, 0.0, 1.0, &A, &r, FALSE, "hi2_pmap_from_points failed");
//     check_pmap("u", 1.0, 1.0, 1.0, &A, &u, TRUE,  "hi2_pmap_from_points failed");
//     
//     /* TO BE COMPLETED !!! */
//     
//     // if (verbose) { fprintf(stderr, "--- i2x2_map_row, i2x2_map_col ---\n"); }
//     // throw_matrix(&A);
//     // i2_throw_cube(trad, &a);
//     // i2x2_map_row(&a, &A, &b);
//     // i2x2_map_col(&A, &a, &c);
//     // i2_zero(&bb);
//     // i2_zero(&cc);
//     // for (i = 0; i < N; i++)
//     //   { for (j = 0; j < N; j++)
//     //       { bb.c[j] += a.c[i] * A.c[i][j];
//     //         cc.c[i] += A.c[i][j] * a.c[j];
//     //       }
//     //   }
//     // r = i2_dist(&b, &bb);
//     // affirm(r < 0.000000001 * i2_norm(&bb), "i2_map_row error");
//     // s = i2_dist(&c, &cc);
//     // affirm(s < 0.000000001 * i2_norm(&cc), "i2_map_col error");
//     // 
//     // if (verbose) { fprintf(stderr, "--- i2x2_mul ---\n"); }
//     // throw_matrix(&A);
//     // throw_matrix(&B);
//     // i2x2_mul(&A, &B, &C);
//     // for (i = 0; i < N; i++)
//     //   { for (j = 0; j < N; j++)
//     //       { double sum = 0.0;
//     //         for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[k][j]; }
//     //         check_eps(C.c[i][j], sum, 0.000000001 * fabs(sum),
//     //           "i2x2_mul error"
//     //         );
//     //       }
//     //   }
//     // 
//     // if (verbose) { fprintf(stderr, "--- i2x2_det ---\n"); }
//     // throw_matrix(&A);
//     // for (i = 0; i < N; i++)
//     //   { int32_t k = (i + 1) % N;
//     //     for (j = 0; j < N; j++)
//     //       { /* Check for linearity */
//     //         r = drandom();
//     //         A.c[i][j] = r;
//     //         rr = i2x2_det(&A);
//     // 
//     //         s = drandom();
//     //         A.c[i][j] = s;
//     //         ss = i2x2_det(&A);
//     // 
//     //         t = drandom();
//     //         A.c[i][j] = r*(1-t) + s*t;
//     //         tt = i2x2_det(&A);
//     //         mag = fabs(rr) + fabs(ss) + fabs(tt);
//     //         check_eps(tt, rr*(1.0 - t) + ss*t, 000000001 * mag,
//     //           "i2x2_det error(1)"
//     //         );
//     //       }
//     // 
//     //     /* Row swap test: */
//     //     r = i2x2_det(&A);
//     //     for (j = 0; j < N; j++)
//     //       { double t = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = t; }
//     //     rr = i2x2_det(&A);
//     //     mag = fabs(r) + fabs(rr);
//     //     check_eps(r, -rr, 000000001 * mag, "i2x2_det error(2)");
//     // 
//     //     /* Col swap test: */
//     //     r = i2x2_det(&A);
//     //     for (j = 0; j < N; j++)
//     //       { double t = A.c[j][i]; A.c[j][i] = A.c[j][k]; A.c[j][k] = t; }
//     //     rr = i2x2_det(&A);
//     //     mag = fabs(r) + fabs(rr);
//     //     check_eps(r, -rr, 000000001 * mag, "i2x2_det error(3)");
//     //   }
//     // 
//     // if (verbose) { fprintf(stderr, "--- i2x2_inv ---\n"); }
//     // throw_matrix(&A);
//     // i2x2_inv(&A, &B);
//     // i2x2_mul(&A, &B, &C);
//     // for (i = 0; i < N; i++)
//     //   { for (j = 0; j < N; j++)
//     //       { double val = (i == j ? 1.0 : 0.0);
//     //         affirm((C.c[i][j] - val) < 000000001, "i2x2_inv error");
//     //       }
//     //   }
//     // 
//     // if (verbose) { fprintf(stderr, "--- i2x2_print ---\n"); }
//     // if (verbose)
//     //   { throw_matrix (&A);
//     //     fprintf(stderr, "A = ");
//     //     i2x2_print(stderr, &A);
//     //     fputc('\n', stderr);
//     //   }
//   }  
// 
// void throw_pmap(hi2_pmap_t *m)
//   {
//     int32_t i, j;
//     i3_t a;
//     for (i = 0; i < NH; i++)
//       { i3_throw_cube(trad, &a);
//         for (j = 0; j < NH; j++) { m->dir.c[i][j] = a.c[j]; }
//       }
//     i3x3_inv(&(m->dir), &(m->inv));
//   }
// 
// void check_pmap
//   ( char *name, 
//     double w, 
//     double x, 
//     double y, 
//     hi2_pmap_t *M, 
//     hi2_point_t *q, 
//     bool_t twirl,
//     char *msg
//   )
//   { 
//     hi2_point_t p, pM; 
//     if (! twirl)
//       { /* {w,x,y} must map to {p} as they are: */
//         p = (hi2_point_t){{{w, x, y}}};
//         pM = hi2_map_point(&p, M); 
//       }
//     else
//       { /* Try flipping the signs of {w,x,y} in all combinations: */
//         double dmin = +INF;
//         int32_t sw, sx, sy;
//         for (sy = -1; sy <= +1; sy += 2)
//           for (sx = -1; sx <= +1; sx += 2)
//             for (sw = -1; sw <= +1; sw += 2)
//               { hi2_point_t u = (hi2_point_t){{{sw*w, sx*x, sy*y}}};
//                 hi2_point_t uM = hi2_map_point(&u, M); 
//                 double d = hi2_pt_pt_diff(&uM, q);
//                 if (d < dmin) { dmin = d; p = u; pM = uM; }
//               }
//       }
//     check_hi2_eps(&p, &pM, q, 0.00000001, msg);
//   }
