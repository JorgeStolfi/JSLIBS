/* hr2test --- test program for hr2.h  */
/* Last edited on 2022-03-01 12:13:48 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <flt.h>
#include <jsrandom.h>
#include <affirm.h>

#include <rn.h>
#include <r3x3.h>
#include <r3.h>
#include <r2x2.h>
#include <r2.h>

#include <hr2.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

/* For the module: */

hr2_point_t hr2_point_throw_cube(void);
  /* Returns a random point in the square {[-1.0 _ +1.0]^{×2}},
    with unit weight coordinate. */

hr2_point_t hr2_point_throw_cube(void)
  { r2_t pc; r2_throw_cube(&pc);
    hr2_point_t p = hr2_from_r2(&pc);
    return p;
  }


/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_hr2(bool_t verbose);
void test_hr2_pmap(bool_t verbose);
void test_hr2_pmap_aff(bool_t verbose);

void throw_aff_map(hr2_pmap_t *M);

void test_hr2_to_from_r2(bool_t verbose);
void test_hr2_pt_pt_diff(bool_t verbose);
void test_hr2_side(bool_t verbose);
void test_hr2_orient(bool_t verbose);
void test_hr2_join(bool_t verbose);
void test_hr2_meet(bool_t verbose);
void test_hr2_point_point_dir(bool_t verbose);
void test_hr2_line_dir(bool_t verbose);
void test_hr2_line_normal(bool_t verbose);
void test_hr2_pmap_scaling(bool_t verbose);
void test_hr2_pmap_congruence_from_point_and_dir(bool_t verbose, bool_t flip);
void test_hr2_pmap_similarity_from_two_points(bool_t verbose, bool_t flip);
void test_hr2_pmap_from_four_points(bool_t verbose);

void test_hr2_pmap_gen_print(bool_t verbose);
void test_hr2_pmap_identity(bool_t verbose);
void test_hr2_pmap_inv(bool_t verbose);
void test_hr2_pmap_compose(bool_t verbose);
void test_hr2_pmap_inv_comp(bool_t verbose);
void test_hr2_pmap_is_affine(bool_t verbose);
void test_hr2_pmap_is_identity(bool_t verbose);
void test_hr2_pmap_line(bool_t verbose);
void test_hr2_pmap_mismatch_sqr(bool_t verbose);
void test_hr2_pmap_point(bool_t verbose);
void test_hr2_pmap_print(bool_t verbose);
void test_hr2_pmap_r2_point(bool_t verbose);
void test_hr2_pmap_aff_from_mat_and_disp(bool_t verbose);
void test_hr2_pmap_aff_from_point_pairs(bool_t verbose);
void test_hr2_pmap_aff_from_points(bool_t verbose);
void test_hr2_pmap_aff_mismatch_sqr(bool_t verbose);
void test_hr2_pmap_deform_sqr(bool_t verbose);
void test_hr2_pmap_rotation(bool_t verbose);
void test_hr2_pmap_rotation_and_scaling(bool_t verbose);
void test_hr2_pmap_scaling(bool_t verbose);
void test_hr2_pmap_translation(bool_t verbose);


    /* TEST: bool_t hr2_pmap_is_affine(r3x3_t *A); */
    /* TEST: bool_t hr2_pmap_is_identity(hr2_pmap_t *M); */
    /* TEST: double hr2_pmap_deform_sqr(r2_t h[], r3x3_t *A); */
    /* TEST: double hr2_pmap_mismatch_sqr(hr2_pmap_t *M, int32_t np, r2_t p1[], r2_t p2[]); */
    /* TEST: double r2_aff_map_mismatch_sqr( r3x3_t *A, r3x3_t *B); */
    /* TEST: hr2_line_t hr2_pmap_line(hr2_line_t *L, hr2_pmap_t *M); */
    /* TEST: hr2_pmap_t hr2_pmap_aff_from_mat_and_disp(r2x2_t *E, r2_t *d); */
    /* TEST: hr2_pmap_t hr2_pmap_aff_from_point_pairs(int32_t np, r2_t p1[], r2_t p2[], double w[]); */
    /* TEST: hr2_pmap_t hr2_pmap_identity(void); */
    /* TEST: hr2_pmap_t hr2_pmap_inv(hr2_pmap_t *M); */
    /* TEST: hr2_pmap_t hr2_pmap_compose(hr2_pmap_t *M, hr2_pmap_t *N); */
    /* TEST: hr2_pmap_t hr2_pmap_inv_comp(hr2_pmap_t *M, hr2_pmap_t *N); */
    /* TEST: hr2_pmap_t hr2_pmap_rotation(double ang); */
    /* TEST: hr2_pmap_t hr2_pmap_rotation_and_scaling(double ang, double scale); */
    /* TEST: hr2_pmap_t hr2_pmap_scaling(r2_t *scale); */
    /* TEST: hr2_pmap_t hr2_pmap_translation(r2_t *vec); */
    /* TEST: hr2_point_t hr2_pmap_point(hr2_point_t *p, hr2_pmap_t *M); */
    /* TEST: r2_t hr2_pmap_r2_point(r2_t *p, r3x3_t *A); */
    /* TEST: typedef enum  */
    /* TEST: typedef struct hr2_line_t { r3_t f; } hr2_line_t;   /\* {f.c[0..2]} are the line's coefficients. *\/ */
    /* TEST: typedef struct hr2_pmap_t { r3x3_t dir; r3x3_t inv; } hr2_pmap_t; */
    /* TEST: typedef struct hr2_point_t { r3_t c; } hr2_point_t; /\* {c.c[0..2]} are the points's coordinates. *\/ */
    /* TEST: typedef void hr2_pmap_opt_report_proc_t (hr2_pmap_t *M, double F); */
    /* TEST: void hr2_pmap_gen_print */
    /* TEST: void hr2_pmap_print (FILE *wr, hr2_pmap_t *M); */

void throw_pmap(hr2_pmap_t *M);
  /* Fills {M} with a random projective map. */

void do_check_eq(double x, double y, char *msg, char *file, int32_t lnum, const char *func);
  /* If {x} and {y} differ, prints them, prints {msg}, and stops. */

#define check_eq(x,y,msg) \
  do_check_eq((x), (y), (msg), __FILE__, __LINE__, __FUNCTION__)
  
void do_check_eps(double x, double y, double eps, char *msg, char *file, int32_t lnum, const char *func);
  /* If {x} and {y} differ by more than {eps}, prints them, prints {msg}, and stops. */

#define check_eps(x, y, eps, msg) \
  do_check_eps((x), (y), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void do_check_hr2_eps
  ( char *name, 
    hr2_point_t *a, 
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

#define check_hr2_eps(name, a, x, y, eps, msg)                                 \
  do_check_hr2_eps((name), (a), (x), (y), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void show_hr2_point(char *tag, hr2_point_t *x);
  /* Prints the pont {x} on a single line, prefixed by {tag}.
    If finite, also prints the Cartesian coordinates. */

/* IMPLEMENTATIONS */

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
  ( char *name,
    hr2_point_t *a, 
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
        if (name != NULL) { fprintf(stderr, "(%s) ", name); }
        if (a != NULL) { rn_gen_print(stderr, NH, a->c.c, "%+2.0f", "[ ", " ", " ]");  }
        fprintf(stderr, "\n");

        show_hr2_point("was", x);
        show_hr2_point("should be", y);
        fprintf(stderr, " diff = %20.16e  max = %+20.16e\n", diff, eps);

        programerror(msg, file, lnum, func);
      }
  }
        
void show_hr2_point(char *tag, hr2_point_t *x)
  { 
    /* Normalize the {hr2_point_t} to unit {\RR^3} norm: */
    hr2_point_t xx;
    r3_dir(&(x->c), &(xx.c));
    fprintf(stderr, "  %10s", tag); 
    rn_gen_print(stderr, NH, xx.c.c, "%+11.8f", "[ ", " ", " ]");
    if (xx.c.c[0] != 0.0) 
      { /* Print the Cartesian coordinates: */
        r2_t p = r2_from_hr2(&xx); r2_gen_print(stderr, &p, "%+11.8f", " = ( ", " ", " )");
      }
    fprintf(stderr, "\n");
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
  );
  /* If {flip} is FALSE, maps the point {[w,x,y]} by {M} and compares it
    with {q}. If {flip} is true, tries reversing the signs of {w,x,y} in
    all combinations, and proceeds as above. In either case, returns
    silently if some attempt produced a match (modulo rounding errors);
    aborts with error {msg} if all attempts failed. */

void check_pmap_point
  ( char *name, 
    hr2_point_t *p, 
    hr2_pmap_t *M, 
    hr2_point_t *q, 
    bool_t flip,
    char *msg
  );
  /* Same as {check_pmap} using the homogeneous coordinates of {p} as {w,x,y}. */


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
    
    test_hr2_pmap_point(verbose);
    test_hr2_pmap_r2_point(verbose);
    test_hr2_pmap_line(verbose);

    test_hr2_pmap_is_identity(verbose);
    test_hr2_pmap_identity(verbose);
    test_hr2_pmap_inv(verbose);
    test_hr2_pmap_compose(verbose);
    test_hr2_pmap_inv_comp(verbose);
    test_hr2_pmap_from_four_points(verbose);
    test_hr2_pmap_mismatch_sqr(verbose);
    test_hr2_pmap_deform_sqr(verbose);

    test_hr2_pmap_print(verbose);
    test_hr2_pmap_gen_print(verbose);
 
    /* TO BE COMPLETED !!! */
    
    /* ??? bool_t hr2_pmap_is_identity(hr2_pmap_t *M); */
    /* ??? hr2_point_t hr2_pmap_point(hr2_point_t *p, hr2_pmap_t *M); */
    /* ??? hr2_line_t hr2_pmap_line(hr2_line_t *L, hr2_pmap_t *M); */
    /* ??? hr2_pmap_t hr2_pmap_translation(r2_t *vec); */
    /* ??? hr2_pmap_t hr2_pmap_rotation(double ang); */

    /* ??? hr2_pmap_t hr2_pmap_comp(hr2_pmap_t *M, hr2_pmap_t *n); */
    /* ??? hr2_pmap_t hr2_pmap_inv(hr2_pmap_t *M); */
    /* ??? hr2_pmap_t hr2_pmap_similarity_from_two_points(hr2_point_t *p, hr2_point_t *q); */

    /* ??? typedef void hr2_pmap_opt_report_proc_t (hr2_pmap_t *M, double F); */
    /* ??? typedef enum  */
    /* ??? double hr2_pmap_mismatch_sqr(hr2_pmap_t *M, int32_t np, r2_t p1[], r2_t p2[]); */
    /* ??? double hr2_pmap_deform_sqr(r2_t h[], r3x3_t *M); */

    // if (verbose) { fprintf(stderr, "--- r2x2_map_row, r2x2_map_col ---\n"); }
    // throw_matrix(&A);
    // r2_throw_cube(&a);
    // r2x2_map_row(&a, &A, &b);
    // r2x2_map_col(&A, &a, &c);
    // r2_zero(&bb);
    // r2_zero(&cc);
    // for (i = 0; i < N; i++)
    //   { for (int32_t j = 0; j < N; j++)
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
    //   { for (int32_t j = 0; j < N; j++)
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
    //     for (int32_t j = 0; j < N; j++)
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
    //     for (int32_t j = 0; j < N; j++)
    //       { double t = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = t; }
    //     rr = r2x2_det(&A);
    //     mag = fabs(r) + fabs(rr);
    //     check_eps(r, -rr, 000000001 * mag, "r2x2_det error(2)");
    // 
    //     /* Col swap test: */
    //     r = r2x2_det(&A);
    //     for (int32_t j = 0; j < N; j++)
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
    //   { for (int32_t j = 0; j < N; j++)
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

void test_hr2_to_from_r2(bool_t verbose)
  { 
    /* TEST: hr2_point_t hr2_from_r2(r2_t *c); */
    /* TEST: r2_t r2_from_hr2(hr2_point_t *p); */

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
    /* TEST: double hr2_pt_pt_diff(hr2_point_t *p, hr2_point_t *q); */

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
      check_eps(dob, dex, 1.0e-8, "hr2_pt_pt_diff error(3)");
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
          check_eps(drt, dex, 1.0e-8, "hr2_pt_pt_diff error(4)");
        }
    }
  }

void test_hr2_side(bool_t verbose)
  { 
    /* TEST: sign_t hr2_side(hr2_point_t *p, hr2_line_t *L);  */

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
    /* TEST: sign_t hr2_orient(hr2_point_t *p, hr2_point_t *q, hr2_point_t *r); */

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
    /* TEST: hr2_line_t hr2_join(hr2_point_t *p, hr2_point_t *q); */

    if (verbose) { fprintf(stderr, "--- hr2_join ---\n"); }
    hr2_point_t p; r3_throw_cube(&(p.c));
    hr2_point_t q; r3_throw_cube(&(q.c));
    hr2_point_t r; r3_throw_cube(&(r.c)); /* Random test point */
    hr2_line_t L = hr2_join(&p, &q);
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
  }
    
void test_hr2_meet(bool_t verbose)
  { 
    /* TEST: hr2_point_t hr2_meet(hr2_line_t *K, hr2_line_t *L); */

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
    /* TEST: r2_t hr2_point_point_dir(hr2_point_t *frm, hr2_point_t *tto); */

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
        { check_eps(upq.c[i], vpq.c[i], 1.0e-12, "hr2_point_point_dir error"); }
    }
  }

void test_hr2_line_dir(bool_t verbose)
  { 
    /* TEST: r2_t hr2_line_dir(hr2_line_t *L); */

    if (verbose) { fprintf(stderr, "--- hr2_line_dir ---\n"); }
    hr2_point_t p; r3_throw_cube(&(p.c)); p.c.c[0] = fabs(p.c.c[0]) + 0.00001;
    hr2_point_t q; r3_throw_cube(&(q.c)); q.c.c[0] = fabs(q.c.c[0]) + 0.00001;
    hr2_line_t L = hr2_join(&p, &q);
    { r2_t dL = hr2_line_dir(&L);
      r2_t eL = hr2_point_point_dir(&p, &q);;
      double tol = 1.0e-12;
      for (int32_t i = 0; i < NC; i++)
        { check_eps(dL.c[i], eL.c[i], tol, "hr2_line_dir error"); }
    }
  }
    
void test_hr2_line_normal(bool_t verbose)
  { 
    /* TEST: r2_t hr2_line_normal(hr2_line_t *L); */

    if (verbose) { fprintf(stderr, "--- hr2_line_normal ---\n"); }
    hr2_line_t L; r3_throw_cube(&(L.f));
    { r2_t nL = hr2_line_normal(&L);
      r2_t mL = (r2_t){{ L.f.c[1], L.f.c[2] }};
      double mLmag = r2_dir(&mL, &mL);
      assert(mLmag != 0);
      double tol = 1.0e-12;
      for (int32_t i = 0; i < NC; i++)
        { check_eps(nL.c[i], mL.c[i], tol, "hr2_line_normal error"); }
    }
  }

void test_hr2_pmap_aff_from_mat_and_disp(bool_t verbose)
  {
    /* TEST: hr2_pmap_t hr2_pmap_aff_from_mat_and_disp(r2x2_t *mat, r2_t *disp); */

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
      { r2_t p2; r2_throw_cube(&p2);
        r2_t q2; r2x2_map_row(&p2, &mat, &q2);
        r2_add(&disp, &q2, &q2);
        hr2_point_t q = hr2_from_r2(&q2);
        double s = dabrandom(0.1, 10.0);
        r3_scale(s, &(q.c), &(q.c));
        check_pmap("p", 1.0, p2.c[0], p2.c[1], &A, &q, FALSE, "hr2_pmap_aff_from_mat_disp failed");
      }
  }
    
void test_hr2_pmap_scaling(bool_t verbose)
  { 
    /* TEST: hr2_pmap_t hr2_pmap_scaling(r2_t *scale) */
        
    if (verbose) { fprintf(stderr, "--- hr2_pmap_scaling ---\n"); }
    r2_t scale; 
    do { r2_throw_cube(&scale); } while ((scale.c[0] == 0) || (scale.c[1] == 0));
    
    hr2_pmap_t A = hr2_pmap_scaling(&scale);

    /* Check whether the map works: */
    r2_t oc = (r2_t){{ 0.0, 0.0 }};
    r2_t pc = (r2_t){{ scale.c[0], 0.0 }};
    r2_t qc = (r2_t){{ 0.0, scale.c[1] }};
    r2_t uc = (r2_t){{ scale.c[0], scale.c[1] }};

    hr2_point_t o = hr2_from_r2(&oc);
    check_pmap("o", 1.0, 0.0, 0.0, &A, &o, FALSE, "hr2_pmap_scaling failed");
    hr2_point_t p = hr2_from_r2(&pc);
    check_pmap("p", 1.0, 1.0, 0.0, &A, &p, FALSE, "hr2_pmap_scaling failed");
    hr2_point_t q = hr2_from_r2(&qc);
    check_pmap("q", 1.0, 0.0, 1.0, &A, &q, FALSE, "hr2_pmap_scaling failed");
    hr2_point_t u = hr2_from_r2(&uc);
    check_pmap("u", 1.0, 1.0, 1.0, &A, &u, FALSE, "hr2_pmap_scaling failed");
  }

void test_hr2_pmap_congruence_from_point_and_dir(bool_t verbose, bool_t flip)
  { 
    bool_t debug = FALSE;

    /* TEST: hr2_pmap_t hr2_pmap_congruence_from_point_and_dir(r2_t *p, r2_t *u, bool_t flip); */
        
    if (verbose) { fprintf(stderr, "--- hr2_pmap_congruence_from_point_and_dir ---\n"); }
    r2_t pc; r2_throw_cube(&pc);
    r2_t u; r2_throw_dir(&u);

    /* Determine the vector {v} that is going to be the image of vector {(0,1)}: */
    r2_t v = (r2_t){{ -u.c[1], +u.c[0] }}; 
    if (flip) { r2_neg(&v, &v); }
    
    if (debug) 
      { fprintf(stderr, "  flip = %c\n", "FT"[flip]);
        r2_gen_print(stderr, &pc, "%12.8f", "  pc = ( ", " ", " )\n");
        r2_gen_print(stderr, &u,  "%12.8f", "  u =  ( ", " ", " )\n");
        r2_gen_print(stderr, &v,  "%12.8f", "  v =  ( ", " ", " )\n");
      }
    
    hr2_pmap_t A = hr2_pmap_congruence_from_point_and_dir(&pc, &u, flip);

    /* Check whether the map works: */
    hr2_point_t p = hr2_from_r2(&pc);
    check_pmap("p", 1.0, 0.0, 0.0, &A, &p, FALSE, "hr2_pmap_congruence_from_point_and_dir failed");
    r2_t qc; r2_add(&pc, &u, &qc);
    hr2_point_t q = hr2_from_r2(&qc);
    check_pmap("q", 1.0, 1.0, 0.0, &A, &q, FALSE, "hr2_pmap_congruence_from_point_and_dir failed");
    r2_t rc; r2_add(&pc, &v, &rc);
    hr2_point_t r = hr2_from_r2(&rc);
    check_pmap("r", 1.0, 0.0, 1.0, &A, &r, FALSE, "hr2_pmap_congruence_from_point_and_dir failed");
  }
    
void test_hr2_pmap_similarity_from_two_points(bool_t verbose, bool_t flip)
  {
    bool_t debug = FALSE;

    /* TEST: hr2_pmap_t hr2_pmap_similarity_from_two_points(r2_t *p, r2_t *q, bool_t flip); */
        
    if (verbose) { fprintf(stderr, "--- hr2_pmap_similarity_from_two_points ---\n"); }
    r2_t pc; r2_throw_cube(&pc);
    r2_t qc; r2_throw_cube(&qc);

    /* Determine the vector {v} that is going to be the image of vector {(0,1)}: */
    r2_t u; r2_sub(&qc, &pc, &u); 
    r2_t v = (r2_t){{ -u.c[1], +u.c[0] }};
    if (flip) { r2_neg(&v, &v); }
    
    if (debug) 
      { fprintf(stderr, "  flip = %c\n", "FT"[flip]);
        r2_gen_print(stderr, &pc, "%12.8f", "  pc = ( ", " ", " )\n");
        r2_gen_print(stderr, &qc, "%12.8f", "  qc = ( ", " ", " )\n");
        r2_gen_print(stderr, &u,  "%12.8f", "  u =  ( ", " ", " )\n");
        r2_gen_print(stderr, &v,  "%12.8f", "  v =  ( ", " ", " )\n");
      }
    
    hr2_pmap_t A = hr2_pmap_similarity_from_two_points(&pc, &qc, flip);

    /* Check whether the map works: */
    hr2_point_t p = hr2_from_r2(&pc);
    check_pmap("p", 1.0, 0.0, 0.0, &A, &p, FALSE, "hr2_pmap_similarity_from_two_points failed");
    hr2_point_t q = hr2_from_r2(&qc);
    check_pmap("q", 1.0, 1.0, 0.0, &A, &q, FALSE, "hr2_pmap_similarity_from_two_points failed");
    r2_t rc; r2_add(&pc, &v, &rc);
    hr2_point_t r = hr2_from_r2(&rc);
    check_pmap("r", 1.0, 0.0, 1.0, &A, &r, FALSE, "hr2_pmap_similarity_from_two_points failed");
  }

void test_hr2_pmap_from_four_points(bool_t verbose)
  { 
    /* TEST: hr2_pmap_t hr2_pmap_from_four_points(hr2_point_t *p, hr2_point_t *q, hr2_point_t *r, hr2_point_t *u); */

    if (verbose) { fprintf(stderr, "--- hr2_pmap_from_four_points ---\n"); }
    hr2_point_t p; r3_throw_cube(&(p.c));
    hr2_point_t q; r3_throw_cube(&(q.c));
    hr2_point_t r; r3_throw_cube(&(r.c));
    hr2_point_t u; r3_throw_cube(&(u.c));
    hr2_pmap_t A = hr2_pmap_from_four_points(&p, &q, &r, &u);
    /* Check whether the map works: */
    
    check_pmap("p", 1.0, 0.0, 0.0, &A, &p, FALSE, "hr2_pmap_from_four_points failed");
    check_pmap("q", 0.0, 1.0, 0.0, &A, &q, FALSE, "hr2_pmap_from_four_points failed");
    check_pmap("r", 0.0, 0.0, 1.0, &A, &r, FALSE, "hr2_pmap_from_four_points failed");
    check_pmap("u", 1.0, 1.0, 1.0, &A, &u, TRUE,  "hr2_pmap_from_four_points failed");
  }

void throw_pmap(hr2_pmap_t *M)
  {
    r3_t a;
    for (int32_t i = 0; i < NH; i++)
      { r3_throw_cube(&a);
        for (int32_t j = 0; j < NH; j++) { M->dir.c[i][j] = a.c[j]; }
      }
    r3x3_inv(&(M->dir), &(M->inv));
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
    check_hr2_eps(name, &p, &pM, q, 0.00000001, msg);
  }

/* ###################################################################### */
/* INTERNAL PROTOTYPES */

void check_num_eps(char *name, double x, double y, double eps, char *msg);
  /* If {x} and {y} differ by more than {eps}, prints {name}, {x}, {y}, and {msg}, and stops. */

void check_r2_eps
  ( char *name, 
    r2_t *a, 
    r2_t *x, 
    r2_t *y, 
    double eps, 
    char *msg
  );
  /* If points {x} and {y} differ by more than {eps} 
    (apart from homogeneous scaling), prints {name}, {a} (if not NULL), 
    {x}, {y}, and {msg}, and stops. */

void check_pmap_r2_point
  ( char *name, 
    r2_t *p, 
    hr2_pmap_t *M, 
    r2_t *q,
    char *msg
  );
  /* Maps the point {p} by {A} and compares it
    with {q}.  Aborts with error {msg} if failed. */
  
void print_pmap(char *name, hr2_pmap_t *M);
  /* Prints {M} to stderr. */

void test_hr2_pmap_aff(bool_t verbose)
  {
    /* Size: */
    if (verbose)
      { fprintf(stderr,
          "sizeof(hr2_pmap_t) = %lu  %d*%d*sizeof(double) = %lu\n",
          sizeof(hr2_pmap_t), NC+1, NC, (NC+1)*NC*sizeof(double)
        );
      }

    test_hr2_pmap_is_affine(verbose);
    test_hr2_pmap_aff_from_mat_and_disp(verbose);

    test_hr2_pmap_translation(verbose);
    test_hr2_pmap_rotation(verbose);
    test_hr2_pmap_scaling(verbose);
    test_hr2_pmap_rotation_and_scaling(verbose);
    test_hr2_pmap_aff_mismatch_sqr(verbose);
    test_hr2_pmap_congruence_from_point_and_dir(verbose, FALSE);
    test_hr2_pmap_congruence_from_point_and_dir(verbose, TRUE);
    test_hr2_pmap_similarity_from_two_points(verbose, FALSE);
    test_hr2_pmap_similarity_from_two_points(verbose, TRUE);
    test_hr2_pmap_aff_from_points(verbose);
    test_hr2_pmap_aff_from_point_pairs(verbose);

    /* TO BE COMPLETED !!! */
    
    if (verbose)
      { 
        /* fprintf(stderr, "!! r2_aff_map_from_XXX NOT TESTED\n"); */
      }
  }
    
void test_hr2_pmap_r2_point(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_r2_point ---\n"); }
    hr2_pmap_t M; throw_aff_map(&M);
    r2_t pc; r2_throw_cube(&pc);
    hr2_point_t p = hr2_from_r2(&pc);
    hr2_point_t q = hr2_pmap_point(&p, &M);
    r2_t qc = (r2_t){{ q.c.c[1]/q.c.c[0],  q.c.c[2]/q.c.c[0] }};
    check_pmap_r2_point("p", &pc, &M, &qc, "hr2_pmap_r2_point failed");
  }

void test_hr2_pmap_inv(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_inv ---\n"); }
    hr2_pmap_t M; throw_aff_map(&M);
    hr2_pmap_t N = hr2_pmap_inv(&M);
    for (int32_t k = 0; k < 5; k++)
      { hr2_point_t p = hr2_point_throw_cube();
        hr2_point_t q = hr2_pmap_point(&p, &M);
        check_pmap_point("p", &q, &N, &p, FALSE, "hr2_pmap_inv failed");
      }
 }

void test_hr2_pmap_compose(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_compose ---\n"); }
    hr2_pmap_t M; throw_aff_map(&M);
    hr2_pmap_t N; throw_aff_map(&N);
    hr2_pmap_t P = hr2_pmap_compose(&M, &N);
    for (int32_t k = 0; k < 5; k++)
      { hr2_point_t p = hr2_point_throw_cube();
        hr2_point_t q = hr2_pmap_point(&p, &M);
        hr2_point_t r = hr2_pmap_point(&q, &N);
        check_pmap_point("p", &p, &P, &r, FALSE, "hr2_pmap_compose failed");
      }
   }

void test_hr2_pmap_inv_comp(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_inv_comp ---\n"); }
    hr2_pmap_t M; throw_aff_map(&M);
    hr2_pmap_t Minv = hr2_pmap_inv(&M);
    hr2_pmap_t N; throw_aff_map(&N);
    hr2_pmap_t P = hr2_pmap_inv_comp(&M, &N);
    for (int32_t k = 0; k < 5; k++)
      { hr2_point_t p = hr2_point_throw_cube();
        hr2_point_t q = hr2_pmap_point(&p, &Minv);
        hr2_point_t r = hr2_pmap_point(&q, &N);
        check_pmap_point("p", &p, &P, &r, FALSE, "hr2_pmap_inv_comp failed");
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
        check_pmap_r2_point("p", &p, &M, &q, "hr2_pmap_translation failed");
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
        check_pmap_r2_point("p", &p, &M, &q, "hr2_pmap_rotation_and_scaling failed");
      }
  }

void test_hr2_pmap_aff_mismatch_sqr(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_aff_mismatch_sqr ---\n"); }
    
    bool_t debug = FALSE;

    hr2_pmap_t M, N;
    throw_aff_map(&M);
    throw_aff_map(&N);
    if (debug)
      { print_pmap("M", &M);
        print_pmap("N", &N);
      }
    double mis2 = hr2_pmap_aff_mismatch_sqr(&M, &N);
    /* The integrand |(M - N)(u(t))|^2 is a sinusoidal function with freqs 0,1,2. */
    /* Can be integrated numerically with 5 or more samples. */
    int32_t nang = 3;
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
    check_num_eps("mis", mis2, cis2, 0.0000001, "hr2_pmap_aff_mismatch_sqr failed");
  }

void test_hr2_pmap_aff_from_points(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr2_pmap_aff_from_points ---\n"); }
    r2_t o; r2_throw_cube(&o);
    r2_t p; r2_throw_cube(&p);
    r2_t q; r2_throw_cube(&q);
    hr2_pmap_t M = hr2_pmap_aff_from_points (&o, &p, &q);
    /* Check whether the map works: */
    r2_t oo = (r2_t){{ 0.0, 0.0 }};
    check_pmap_r2_point("o", &oo, &M, &o, "hr2_pmap_aff_from_points failed");
    r2_t pp = (r2_t){{ 1.0, 0.0 }};
    check_pmap_r2_point("p", &pp, &M, &p, "hr2_pmap_aff_from_points failed");
    r2_t qq = (r2_t){{ 0.0, 1.0 }};
    check_pmap_r2_point("q", &qq, &M, &q, "hr2_pmap_aff_from_points failed");
  }

void test_hr2_pmap_aff_from_point_pairs(bool_t verbose)
  {
    /* Decide how many data point pairs to use: */
    
    int32_t np;
    if (drandom() < 0.250)
      { np = int32_abrandom(0,3); }
    else
      { np = int32_abrandom(4,20); }
     
    if (verbose) { fprintf(stderr, "--- hr2_pmap_aff_from_point_pairs  np = %d ---\n", np); }
     
    bool_t exact = (np <= 3);
    verbose = verbose | exact;

    bool_t debug = verbose;

    if (debug) { fprintf(stderr, "np = %d\n", np); }
    
    /* Pick a set of original points: */
    r2_t p[np];
    if ((np >= 1) && (np <= 3)) { p[0] = (r2_t){{ 5.0, 5.0 }}; }
    if ((np >= 2) && (np <= 3)) { p[1] = (r2_t){{ 7.0, 5.0 }}; }
    if (np == 3) { p[2] = (r2_t){{ 5.0, 7.0 }}; }
    if (np >= 4)
      { 
        for (int32_t ip = 0; ip < np; ip++) { r2_throw_cube(&(p[ip])); }
      }
    
    /* Pick a set of weights: */
    double w[np];
    for (int32_t ip = 0; ip < np; ip++) { w[ip] = dabrandom(0.1, 1.0); }

    /* Choose the expected map {M}: */
    hr2_pmap_t M; 
    if (np == 0)
      { /* Identity: */
        M = hr2_pmap_identity();
      }
    else if (np == 1)
      { /* M translation: */
        r2_t v; r2_throw_cube(&v);
        M = hr2_pmap_translation(&v);
      }
    else if (np == 2)
      { /* M rotation, scaling, and translation: */
        double ang = dabrandom(0.0, 2.0)*M_PI;
        double scale = dabrandom(0.5, 2.0);
        hr2_pmap_t R = hr2_pmap_rotation_and_scaling(ang, scale);
        r2_t v; r2_throw_cube(&v);
        hr2_pmap_t T = hr2_pmap_translation(&v);
        M = hr2_pmap_compose(&R, &T);
      }
    else
      { throw_aff_map(&M); }
    
    /* Map the points through {M}, and add perturbation: */
    r2_t q[np];
    for (int32_t ip = 0; ip < np; ip++) 
      { q[ip] = hr2_pmap_r2_point(&(p[ip]), &M);
        if (! exact)
          { r2_t d; r2_throw_cube(&d);
            r2_mix(1.0, &(q[ip]), 0.1, &d, &(q[ip]));
          }
      }
    
    /* Compute the map: */
    hr2_pmap_t N = hr2_pmap_aff_from_point_pairs(np, p, q, w);
    
    if (debug)
      { print_pmap("M", &M);
        print_pmap("N", &N);
      }
    /* Compare the maps with {hr2_pmap_aff_mismatch_sqr}: */
    double mis2 = hr2_pmap_aff_mismatch_sqr(&M, &N);
    double mis = sqrt(mis2); 
    if (exact)
      { check_num_eps("mis", mis, 0.0, 0.0000001, "hr2_pmap_aff_from_point_pairs failed"); }
    else 
      { /* Compare the maps by testing with given points: */
        assert (np > 0);
        double sum_d2A = 0.0;
        double sum_d2B = 0.0;
        for (int32_t ip = 0; ip < np; ip++)
          { r2_t pi = p[ip]; /* M source point. */
            r2_t qi = q[ip]; /* Given destination point. */
            r2_t piA = hr2_pmap_r2_point(&pi, &M); /* Where the map {M} sent it. */
            double d2A = r2_dist_sqr(&piA, &qi);
            sum_d2A += d2A;
            r2_t piB = hr2_pmap_r2_point(&pi, &N); /* Where the solution {N} sent it. */
            double d2B = r2_dist_sqr(&piB, &qi);
            sum_d2B += d2B;
          }
        double dA = sqrt(sum_d2A/np); /* RMS error of original data points rel to {M}. */
        double dB = sqrt(sum_d2B/np); /* RMS error of original data points rel to {N}. */
        if (verbose) { fprintf(stderr, "rms error M = %12.7f N = %12.7f\n", dA, dB); }
        double bad = fmax(0.0, dB - dA); /* Positive if {N} is worse than {M}, 0 if better. */
        check_num_eps("bad", bad, 0.0, np*0.0000001, "hr2_pmap_aff_from_point_pairs failed");
      }
  }

void test_hr2_pmap_point(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_point} NOT TESTED\n"); }

void test_hr2_pmap_line(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_line} NOT TESTED\n"); }

void test_hr2_pmap_is_identity(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_is_identity} NOT TESTED\n"); }

void test_hr2_pmap_identity(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_identity} NOT TESTED\n"); }

void test_hr2_pmap_mismatch_sqr(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_mismatch_sqr} NOT TESTED\n"); }

void test_hr2_pmap_deform_sqr(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_deform_sqr} NOT TESTED\n"); }

void test_hr2_pmap_print(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_print} NOT TESTED\n"); }

void test_hr2_pmap_gen_print(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_gen_print} NOT TESTED\n"); }

void test_hr2_pmap_is_affine(bool_t verbose)
  { fprintf(stderr, "!! {test_hr2_pmap_is_affine} NOT TESTED\n"); }

void test_hr2_pmap_rotation(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_rotation} NOT TESTED\n"); }

void throw_aff_map(hr2_pmap_t *M)
  {
    for (int32_t i = 0; i < NH; i++)
      { r2_t p;
        r2_throw_cube(&p);
        M->dir.c[i][0] = (i == 0 ? 1.0 : 0.0);
        M->dir.c[i][1] = p.c[0];
        M->dir.c[i][2] = p.c[1];
      }
    r3x3_inv(&(M->dir), &(M->inv));
  }
  
void print_pmap(char *name, hr2_pmap_t *M)
  { fprintf(stderr, "%s = ", name);
    hr2_pmap_gen_print(stderr, M, "%+10.6f", NULL,  NULL,NULL,NULL, NULL,NULL,NULL, "\n");
  }

void check_pmap_point
  ( char *name, 
    hr2_point_t *p, 
    hr2_pmap_t *M, 
    hr2_point_t *q, 
    bool_t flip,
    char *msg
  )
  { check_pmap(name, p->c.c[0], p->c.c[1], p->c.c[2], M, q, flip, msg); }

void check_pmap_r2_point
  ( char *name, 
    r2_t *p, 
    hr2_pmap_t *M, 
    r2_t *q,
    char *msg
  )
  { 
    r2_t pA = hr2_pmap_r2_point(p, M); 
    check_r2_eps(name, p, &pA, q, 0.00000001, msg);
  }

void check_num_eps(char *name, double x, double y, double eps, char *msg)
  { double diff = fabs(x - y);
    if (diff > eps)
      { fprintf(stderr, " ** %s: %+20.16e %+20.16e", name, x, y);
        fprintf(stderr, " diff = %+20.16e  max = %+20.16e - %s\n", diff, eps, msg);
        exit(1);
      }
  }

void check_r2_eps
  ( char *name, 
    r2_t *a, 
    r2_t *x, 
    r2_t *y, 
    double eps, 
    char *msg
  )
  {
    double diff = r2_dist(x, y);
    if (diff > eps)
      { fprintf(stderr, " ** %s: ", name);
        if (a != NULL)
          { rn_gen_print(stderr, NC, a->c, "%+2.0f", "[ ", " ", " ]");
            fprintf(stderr, " -> "); 
          }
        rn_gen_print(stderr, NC, x->c, "%+11.8f", "[ ", " ", " ]");
        fprintf(stderr, " should be "); 
        rn_gen_print(stderr, NC, y->c, "%+11.8f", "[ ", " ", " ]");
        fprintf(stderr, " sdiff = %20.16e  max = %+20.16e - %s\n", diff, eps, msg);
        exit(1);
      }
  }

