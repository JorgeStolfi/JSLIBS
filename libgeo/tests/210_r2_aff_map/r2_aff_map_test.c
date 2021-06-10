/* hr2test --- test program for hr2.h  */
/* Last edited on 2021-06-09 19:51:49 by jstolfi */

#define _GNU_SOURCE_
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <flt.h>
#include <jsrandom.h>
#include <affirm.h>

#include <r2.h>
#include <r2x2.h>
#include <rn.h>

#include <r2_aff_map.h>

#define NC 2
 
/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);
void test_r2_aff_map(int32_t i);
void throw_aff_map(r2_aff_map_t *A);

void test_r2_aff_map_apply(int32_t i, bool_t verbose);
void test_r2_aff_map_invert(int32_t i, bool_t verbose);
void test_r2_aff_map_compose(int32_t i, bool_t verbose);
void test_r2_aff_map_translation(int32_t i, bool_t verbose);
void test_r2_aff_map_rot_scale(int32_t i, bool_t verbose);
void test_r2_aff_map_disp_sqr(int32_t i, bool_t verbose);
void test_r2_aff_map_mismatch_sqr(int32_t i, bool_t verbose);
void test_r2_aff_map_from_points(int32_t i, bool_t verbose);
void test_r2_aff_map_from_point_pairs(int32_t i, bool_t verbose);

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

void check_aff_map
  ( char *name, 
    r2_t *p, 
    r2_aff_map_t *A, 
    r2_t *q,
    char *msg
  );
  /* Maps the point {p} by {A} and compares it
    with {q}.  Aborts with error {msg} if failed. */
  
void print_aff_map(char *name, r2_aff_map_t *A);
  /* Prints {A} to stderr. */

int32_t main (int32_t argc, char **argv)
  {
    int32_t i;
    srand(1993);

    for (i = 0; i < 100; i++) test_r2_aff_map(i);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_r2_aff_map(int32_t i)
  {
    bool_t verbose = (i < 6);

    /* Size: */
    if (verbose)
      { fprintf(stderr,
          "sizeof(r2_aff_map_t) = %lud  %d*%d*sizeof(double) = %lu\n",
          sizeof(r2_aff_map_t), NC+1, NC, (NC+1)*NC*sizeof(double)
        );
      }

    test_r2_aff_map_apply(i, verbose);
    test_r2_aff_map_invert(i, verbose);
    test_r2_aff_map_compose(i, verbose);
    test_r2_aff_map_translation(i, verbose);
    test_r2_aff_map_rot_scale(i, verbose);
    test_r2_aff_map_disp_sqr(i, verbose);
    test_r2_aff_map_mismatch_sqr(i, verbose);
    test_r2_aff_map_from_points(i, verbose);
    test_r2_aff_map_from_point_pairs(i, verbose);

    /* TO BE COMPLETED !!! */
    
    if (verbose)
      { 
        /* fprintf(stderr, "!! r2_aff_map_from_XXX NOT TESTED\n"); */
      }
  }
    
void test_r2_aff_map_apply(int32_t i, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_aff_map_apply ---\n"); }
    r2_aff_map_t A; throw_aff_map(&A);
    r2_t p; r2_throw_cube(&p);
    r2_t q; r2x2_map_row(&p, &(A.mat), &q);
    r2_add(&(A.disp), &q, &q);
    check_aff_map("p", &p, &A, &q, "r2_aff_map_apply failed");
  }

void test_r2_aff_map_invert(int32_t i, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_aff_map_invert ---\n"); }
    r2_aff_map_t A; throw_aff_map(&A);
    r2_aff_map_t B; r2_aff_map_invert(&A, &B);
    for (int32_t k = 0; k < 5; k++)
      { r2_t p; r2_throw_cube(&p);
        r2_t q; r2_aff_map_apply(&p, &A, &q);
        check_aff_map("p", &q, &B, &p, "r2_aff_map_invert failed");
      }
 }

void test_r2_aff_map_compose(int32_t i, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_aff_map_compose ---\n"); }
    r2_aff_map_t A; throw_aff_map(&A);
    r2_aff_map_t B; throw_aff_map(&B);
    r2_aff_map_t C; r2_aff_map_compose(&A, &B, &C);
    for (int32_t k = 0; k < 5; k++)
      { r2_t p; r2_throw_cube(&p);
        r2_t q; r2_aff_map_apply(&p, &A, &q);
        r2_t r; r2_aff_map_apply(&q, &B, &r);
        check_aff_map("p", &p, &C, &r, "r2_aff_map_compose failed");
      }
   }

void test_r2_aff_map_translation(int32_t i, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_aff_map_translation ---\n"); }
    r2_t r; r2_throw_cube(&r);
    r2_aff_map_t A = r2_aff_map_translation(&r);
    for (int32_t k = 0; k < 5; k++)
      { r2_t p; r2_throw_cube(&p);
        r2_t q; r2_add(&r, &p, &q);
        check_aff_map("p", &p, &A, &q, "r2_aff_map_translation failed");
      }
  }

void test_r2_aff_map_rot_scale(int32_t i, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_aff_map_rot_scale ---\n"); }
    double ang = dabrandom(-1.58, +1.58);
    double scale = dabrandom(0.5, 2.0);
    r2_aff_map_t A = r2_aff_map_rot_scale(ang, scale);
    double ca = cos(ang);
    double sa = sin(ang);
    for (int32_t k = 0; k < 5; k++)
      { r2_t p; r2_throw_cube(&p);
        r2_t q; 
        q.c[0] = + ca*scale*p.c[0] - sa*scale*p.c[1];
        q.c[1] = + sa*scale*p.c[0] + ca*scale*p.c[1];
        check_aff_map("p", &p, &A, &q, "r2_aff_map_rot_scale failed");
      }
  }

void test_r2_aff_map_disp_sqr(int32_t i, bool_t verbose)
  {
    bool_t debug = (i == 0);
    if (verbose) { fprintf(stderr, "--- r2_aff_map_disp_sqr ---\n"); }
    
    r2_aff_map_t A, B, R;
    if (i == 0)
      { A = (r2_aff_map_t){ .mat = (r2x2_t){{ { 1, 2 }, { 3, 4 } }}, .disp = (r2_t){{ 5, 6 }} };
        B = (r2_aff_map_t){ .mat = (r2x2_t){{ { 4, 5 }, { 6, 7 } }}, .disp = (r2_t){{ 8, 9 }} };
        R = (r2_aff_map_t){ .mat = (r2x2_t){{ { 3, 3 }, { 3, 3 } }}, .disp = (r2_t){{ 3, 3 }} };
      }
    else
      {
        throw_aff_map(&A);
        throw_aff_map(&B);
        throw_aff_map(&R);
      }
    int32_t iz = int32_abrandom(0, 1);
    int32_t jz = int32_abrandom(0, 1);
    int32_t kz = int32_abrandom(0, 1);
    R.mat.c[iz][jz] = 0;
    R.disp.c[kz] = 0;
    if (debug)
      { print_aff_map("A", &A);
        print_aff_map("B", &B);
        print_aff_map("R", &R);
      }
    double dabs2, drel2;
    r2_aff_map_disp_sqr(&A, &B, &R, &dabs2, &drel2);
    for (int32_t k = 0; k < 5; k++)
      { double cabs2 = 0, crel2 = 0;
        for (int32_t j = 0; j < NC;  j++)
          { for (int32_t i = 0; i < NC;  i++)
              { double rij = R.mat.c[i][j];
                if (rij != 0.0)
                  { double d = A.mat.c[i][j] - B.mat.c[i][j];
                    cabs2 += d*d;
                    crel2 += (d/rij)*(d/rij);
                  }
              }
            double rj = R.disp.c[j];
            if (rj !=0)
              { double d = A.disp.c[j] - B.disp.c[j];
                cabs2 += d*d;
                crel2 += (d/rj)*(d/rj);
              }
          }
        check_num_eps("drel2", drel2, crel2, 0.0000001, "r2_aff_map_disp_sqr failed");
        check_num_eps("dabs2", dabs2, cabs2, 0.0000001, "r2_aff_map_disp_sqr failed");
      }
  }

void test_r2_aff_map_mismatch_sqr(int32_t i, bool_t verbose)
  {
    bool_t debug = (i < 3);
    if (verbose) { fprintf(stderr, "--- r2_aff_map_mismatch_sqr ---\n"); }
    
    r2_aff_map_t A, B;
    if (i == 0)
      { A = (r2_aff_map_t){ .mat = (r2x2_t){{ { 10, 0 }, { 0, 10 } }}, .disp = (r2_t){{ 5, 6 }} };
        B = (r2_aff_map_t){ .mat = (r2x2_t){{ { 20, 0 }, { 0, 20 } }}, .disp = (r2_t){{ 5, 6 }} };
      }
    else if (i == 1)
      { A = (r2_aff_map_t){ .mat = (r2x2_t){{ { +10,   0 }, {   0, +10 } }}, .disp = (r2_t){{ 5, 6 }} };
        B = (r2_aff_map_t){ .mat = (r2x2_t){{ {   0, +10 }, { -10,   0 } }}, .disp = (r2_t){{ 5, 6 }} };
      }
    else if (i == 2)
      { A = (r2_aff_map_t){ .mat = (r2x2_t){{ { 10, 0 }, { 0,  0 } }}, .disp = (r2_t){{  5,  6 }} };
        B = (r2_aff_map_t){ .mat = (r2x2_t){{ { 10, 0 }, { 0, 10 } }}, .disp = (r2_t){{  5,  6 }} };
      }
    else
      {
        throw_aff_map(&A);
        throw_aff_map(&B);
      }
    if (debug)
      { print_aff_map("A", &A);
        print_aff_map("B", &B);
      }
    double mis2 = r2_aff_map_mismatch_sqr(&A, &B);
    /* The integrand |(A - B)(u(t))|^2 is a sinusoidal function with freqs 0,1,2. */
    /* Can be integrated numerically with 5 or more samples. */
    int32_t nang = 3;
    double sum_d2 = 0; 
    for (int32_t i = 0; i < nang; i++)
      { double ang = 2*M_PI*((double)i)/((double)nang);
        r2_t p = (r2_t){{ cos(ang), sin(ang) }};
        r2_t q; r2_aff_map_apply(&p, &A, &q);
        r2_t r; r2_aff_map_apply(&p, &B, &r);
        double d2 = r2_dist_sqr(&q, &r);
        sum_d2 += d2;
      }
    double cis2 = sum_d2/nang;
    check_num_eps("mis", mis2, cis2, 0.0000001, "r2_aff_map_mismatch_sqr failed");
  }

void test_r2_aff_map_from_points(int32_t i, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- r2_aff_map_from_points ---\n"); }
    r2_t o; r2_throw_cube(&o);
    r2_t p; r2_throw_cube(&p);
    r2_t q; r2_throw_cube(&q);
    r2_aff_map_t A = r2_aff_map_from_points (&o, &p, &q);
    /* Check whether the map works: */
    r2_t oo = (r2_t){{ 0.0, 0.0 }};
    check_aff_map("o", &oo, &A, &o, "r2_aff_map_from_points failed");
    r2_t pp = (r2_t){{ 1.0, 0.0 }};
    check_aff_map("p", &pp, &A, &p, "r2_aff_map_from_points failed");
    r2_t qq = (r2_t){{ 0.0, 1.0 }};
    check_aff_map("q", &qq, &A, &q, "r2_aff_map_from_points failed");
  }

void test_r2_aff_map_from_point_pairs(int32_t i, bool_t verbose)
  {
    /* Decide how many data point pairs to use: */
    int32_t np = (i <= 3 ? i : int32_abrandom(4,20));
     
     
    /* bool_t debug = verbose; */
    bool_t debug = TRUE;

    bool_t exact = (np <= 3);
    verbose = verbose | exact | debug;

    if (verbose) { fprintf(stderr, "--- r2_aff_map_from_point_pairs ---\n"); }
    
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

    /* Pick a goalmap: */
    r2_aff_map_t A; 
    if (np == 0)
      { /* Identity: */
        r2x2_ident(&(A.mat));
        A.disp = (r2_t){{ 0.0, 0.0 }};
      }
    else if (np == 1)
      { /* A translation: */
        r2_t v; r2_throw_cube(&v);
        A = r2_aff_map_translation(&v);
      }
    else if (np == 2)
      { /* A rotation, scaling, and translation: */
        double ang = dabrandom(0.0, 2.0)*M_PI;
        double scale = dabrandom(0.5, 2.0);
        r2_aff_map_t R = r2_aff_map_rot_scale(ang, scale);
        r2_t v; r2_throw_cube(&v);
        r2_aff_map_t T = r2_aff_map_translation(&v);
        r2_aff_map_compose(&R, &T, &A);
      }
    else
      { throw_aff_map(&A); }
    
    /* Map the points through {A}, and add perturbation: */
    r2_t q[np];
    for (int32_t ip = 0; ip < np; ip++) 
      { r2_aff_map_apply(&(p[ip]), &A, &(q[ip]));
        if (! exact)
          { r2_t d; r2_throw_cube(&d);
            r2_mix(1.0, &(q[ip]), 0.1, &d, &(q[ip]));
          }
      }
    
    /* Compute the map: */
    r2_aff_map_t B = r2_aff_map_from_point_pairs(np, p, q, w);
    
    if (debug)
      { print_aff_map("A", &A);
        print_aff_map("B", &B);
      }

    /* Compare the maps with {r2_aff_map_mismatch_sqr}: */
    double mis2 = r2_aff_map_mismatch_sqr(&A, &B);
    double mis = sqrt(mis2); 
    if (exact)
      { check_num_eps("mis", mis, 0.0, 0.0000001, "r2_aff_map_from_point_pairs failed"); }
    else 
      { /* Compare the maps by testing with given points: */
        assert (np > 0);
        double sum_d2A = 0.0;
        double sum_d2B = 0.0;
        for (int32_t ip = 0; ip < np; ip++)
          { r2_t pi = p[ip]; /* A source point. */
            r2_t qi = q[ip]; /* Given destination point. */
            r2_t piA; r2_aff_map_apply(&pi, &A, &piA); /* Where the map {A} sent it. */
            double d2A = r2_dist_sqr(&piA, &qi);
            sum_d2A += d2A;
            r2_t piB; r2_aff_map_apply(&pi, &B, &piB); /* Where the solution {B} sent it. */
            double d2B = r2_dist_sqr(&piB, &qi);
            sum_d2B += d2B;
          }
        double dA = sqrt(sum_d2A/np); /* RMS error of original data points rel to {A}. */
        double dB = sqrt(sum_d2B/np); /* RMS error of original data points rel to {B}. */
        if (verbose) { fprintf(stderr, "rms error A = %12.7f B = %12.7f\n", dA, dB); }
        double bad = fmax(0.0, dB - dA); /* Positive if {B} is worse than {A}, 0 if better. */
        check_num_eps("bad", bad, 0.0, np*0.0000001, "r2_aff_map_from_point_pairs failed");
      }
  }

void throw_aff_map(r2_aff_map_t *A)
  {
    for (int32_t i = 0; i < NC; i++)
      { r2_t p;
        r2_throw_cube(&p);
        A->mat.c[i][0] = p.c[0];
        A->mat.c[i][1] = p.c[1];
      }
    r2_throw_cube(&(A->disp));
  }
  
void print_aff_map(char *name, r2_aff_map_t *A)
  { fprintf(stderr, "%s = ", name);
    r2_aff_map_gen_print(stderr, A, "%+10.6f", "%+10.6f", "[ "," "," ]","[ "," "," ]"," + ");
    fprintf(stderr, "\n");
  }

void check_aff_map
  ( char *name, 
    r2_t *p, 
    r2_aff_map_t *A, 
    r2_t *q,
    char *msg
  )
  { 
    r2_t pA; 
    r2_aff_map_apply(p, A, &pA); 
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

