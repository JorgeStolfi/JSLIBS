#define PROG_NAME "hr2_pmap_enc_dec_test"
#define PROG_DESC "test of {hr2_pmap*_encode.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-09-17 19:42:48 by stolfi */ 
/* Created on 2024-09-07 by J. Stolfi, UNICAMP */

#define hr2_pmap_enc_dec_test_COPYRIGHT \
  "Copyright � 2024  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <ix.h>
#include <r2.h>
#include <hr2.h>
#include <hr2_pmap.h>
#include <rn.h>
#include <r3x3.h>
#include <bool.h>
#include <jsrandom.h>
#include <affirm.h>
#include <hr2_test_tools.h>
#include <hr2_pmap_test_tools.h>

#include <hr2_pmap_translation_encode.h>
#include <hr2_pmap_congruence_encode.h>
#include <hr2_pmap_similarity_encode.h>
#include <hr2_pmap_affine_encode.h>
#include <hr2_pmap_generic_encode.h>
#include <hr2_pmap_encode.h>
#include <hr2_pmap_throw_type.h>

void test_hr2_pmap_encode_num_parameters(bool_t verbose);

void test_hr2_pmap_translation_encode__hr2_pmap_translation_decode(bool_t verbose);
void test_hr2_pmap_congruence_encode__hr2_pmap_congruence_decode(bool_t verbose);
void test_hr2_pmap_similarity_encode__hr2_pmap_similarity_decode(bool_t verbose);
void test_hr2_pmap_affine_encode__hr2_pmap_affine_decode(bool_t verbose);
void test_hr2_pmap_generic_encode__hr2_pmap_generic_decode(bool_t verbose);
void test_hr2_pmap_throw_type__hr2_pmap_is_type(bool_t verbose);
void test_hr2_pmap_encode__hr2_pmap_decode(bool_t verbose);

int32_t main(int32_t argn, char **argv);

void hpedt_print_map(char *name, hr2_pmap_t *M);
  /* Prints the map {M} (direct and inverse) to {stderr},
    prefixed by a line"  {name} = ".  */

r3x3_t hpedt_throw_r3x3(double mag);
  /* Returns a random 3x3 matrix with entries in {[-mag _ +mag]}. */

void hpedt_print_encoding(int32_t ny, double y[]);
  /* Prints {y[0..ny-1]} on one line, preceded by "encoding". */

void hpedt_check_same_map(hr2_pmap_t *M, hr2_pmap_t *N);
  /* Compares the original map {M} with the encoded and decoded map {N}. 
    The matrices should be the same apart from homogeneous scaling. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srand(4615*417);
    
    test_hr2_pmap_encode_num_parameters(TRUE);
  
    for (int32_t i = 0; i < 20; i++)
      { bool_t verbose = TRUE;
        test_hr2_pmap_translation_encode__hr2_pmap_translation_decode(verbose);
        test_hr2_pmap_congruence_encode__hr2_pmap_congruence_decode(verbose);
        test_hr2_pmap_similarity_encode__hr2_pmap_similarity_decode(verbose);
        test_hr2_pmap_affine_encode__hr2_pmap_affine_decode(verbose);
        test_hr2_pmap_generic_encode__hr2_pmap_generic_decode(verbose);
        test_hr2_pmap_throw_type__hr2_pmap_is_type(verbose);
        test_hr2_pmap_encode__hr2_pmap_decode(verbose);
      }
    return 0;
  }

void test_hr2_pmap_encode_num_parameters(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "> --- %s ---\n", __FUNCTION__); }

    demand(0 == hr2_pmap_encode_num_parameters(hr2_pmap_type_IDENTITY), "IDENTITY");
    demand(2 == hr2_pmap_encode_num_parameters(hr2_pmap_type_TRANSLATION), "TRANSLATION");
    demand(3 == hr2_pmap_encode_num_parameters(hr2_pmap_type_CONGRUENCE), "CONGRUENCE");
    demand(4 == hr2_pmap_encode_num_parameters(hr2_pmap_type_SIMILARITY), "SIMILARITY");
    demand(6 == hr2_pmap_encode_num_parameters(hr2_pmap_type_AFFINE), "AFFINE");
    demand(9 == hr2_pmap_encode_num_parameters(hr2_pmap_type_GENERIC), "GENERIC");

    if (verbose) { fprintf(stderr, "< --- %s ---\n", __FUNCTION__); }
  }

void test_hr2_pmap_translation_encode__hr2_pmap_translation_decode(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "> --- %s ---\n", __FUNCTION__); }
   
    /* Create a translation {M}: */
    double xmax = 10.0, ymax = 20.0;
    r2_t disp = (r2_t){{ dabrandom(-xmax,+xmax), dabrandom(-ymax,+ymax) }};
    hr2_pmap_t M = hr2_pmap_translation(&disp);
    r3x3_normalize(&(M.inv));  /* Assumes {N.dir} is normalized. */
    if (verbose) { hpedt_print_map("original", &M); }

    double eqtol = 0.0; /* Tolerance for type matching. */
    demand(hr2_pmap_is_type(&M, hr2_pmap_type_TRANSLATION, +1, eqtol), "wrong map type or sign");

    /* Test encoding: */
    int32_t ny = 2;
    double y[ny];
    hr2_pmap_translation_encode(&M, y);
    if (verbose) { hpedt_print_encoding(ny, y); }
    demand(fabs(y[0] - disp.c[0]) <= eqtol, "{hr2_pmap_translation_encode} failed (0)");
    demand(fabs(y[1] - disp.c[1]) <= eqtol, "{hr2_pmap_translation_encode} failed (1)");

    /* Test decoding: */
    hr2_pmap_t N = hr2_pmap_identity();
    hr2_pmap_translation_decode(y, &N);
    r3x3_normalize(&(N.inv));  /* Assumes {N.dir} is normalized. */
    if (verbose) { hpedt_print_map("decoded", &N); }

    demand(hr2_pmap_is_type(&N, hr2_pmap_type_TRANSLATION, +1, eqtol), "wrong map type or sign");

    hpedt_check_same_map(&M, &N);

    if (verbose) { fprintf(stderr, "< --- %s ---\n", __FUNCTION__); }
  }

void test_hr2_pmap_congruence_encode__hr2_pmap_congruence_decode(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "> --- %s ---\n", __FUNCTION__); }
   
    double eqtol = 1.0e-13; /* Tolerance for type matching. */

    /* Choose displacement {disp} and rotation angle {ang}: */
    double xmax = 10.0, ymax = 20.0;
    double angmax = M_PI; /* Better be less than {pi/2}. */
    r2_t disp = (r2_t){{ dabrandom(-xmax,+xmax), dabrandom(-ymax,+ymax) }};
    double ang = dabrandom(-angmax, +angmax);
    r2_t u = (r2_t){{ cos(ang), sin(ang) }};
    
    /* Generate a congruence {M0} with positive handedness: */
    hr2_pmap_t M0 = hr2_pmap_congruence_from_point_and_dir(&disp, &u, +1);
    if (verbose) { hpedt_print_map("reference", &M0); }
    demand(hr2_pmap_is_type(&M0, hr2_pmap_type_CONGRUENCE, +1, eqtol), "wrong map type or sign");
    demand(r3x3_det(&(M0.dir)) > 0, "{hr2_pmap_congruence_from_point_and_dir} wrong handedness (1)");
    
    /* Get the encoding {y0[0..2]} of {M0}: */
    int32_t ny = 3;
    double y0[ny];
    hr2_pmap_congruence_encode(&M0, y0);
    if (verbose) { hpedt_print_encoding(ny, y0); }
    demand(fabs(y0[0] - disp.c[0]) <= eqtol, "{hr2_pmap_congruence_encode} failed (0)");
    demand(fabs(y0[1] - disp.c[1]) <= eqtol, "{hr2_pmap_congruence_encode} failed (1)");
    demand(fabs(y0[2] - ang) <= eqtol, "{hr2_pmap_congruence_encode} failed (2)");
 
    for (sign_t sgn = +1; sgn >= -1; sgn -= 2)
      {
        if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
        if (verbose) { fprintf(stderr, "  testing with sgn = %+d\n", sgn); }
        
        hr2_pmap_t M = M0;
        hr2_pmap_set_sign(&M, sgn);
        if (verbose) { hpedt_print_map("hand-made", &M); }
        demand(sgn*r3x3_det(&(M.dir)) > 0, "{hr2_pmap_congruence_from_point_and_dir} wrong handedness (2)");

        /* Test encoding: */
        double y[ny];
        hr2_pmap_congruence_encode(&M, y);
        if (verbose) { hpedt_print_encoding(ny, y); }
        for (int32_t jy = 0; jy < ny; jy++)
          { demand(fabs(y[jy] - y0[jy]) <= eqtol, "{hr2_pmap_congruence_encode} failed (3)"); }

        /* Test decoding: */
        hr2_pmap_t N = hr2_pmap_identity();
        hr2_pmap_congruence_decode(y, &N);
        if (verbose) { hpedt_print_map("decoded", &N); }
        hpedt_check_same_map(&M0, &N);

        if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
      }

    if (verbose) { fprintf(stderr, "< --- %s ---\n", __FUNCTION__); }
  }

void test_hr2_pmap_similarity_encode__hr2_pmap_similarity_decode(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "> --- %s ---\n", __FUNCTION__); }
   
    double eqtol = 1.0e-13; /* Tolerance for type matching. */

    /* Choose two points {p,q}: */
    double xmax = 10.0, ymax = 20.0;
    r2_t p = (r2_t){{ dabrandom(-xmax,+xmax), dabrandom(-ymax,+ymax) }};
    r2_t u = (r2_t){{ 1.0 + dabrandom(-xmax,+xmax), dabrandom(-ymax,+ymax) }};
    double ang = atan2(u.c[1], u.c[0]);
    r2_t q; r2_add(&p, &u, &q);
     
    /* Generate a similarity {M0} with positive handedness: */
    hr2_pmap_t M0 = hr2_pmap_similarity_from_two_points(&p, &q, +1);
    if (verbose) { hpedt_print_map("reference", &M0); }
    demand(hr2_pmap_is_type(&M0, hr2_pmap_type_SIMILARITY, +1, eqtol), "wrong map type or sign");

    /* Get the encoding {y0[0..3]} of {M0}: */
    int32_t ny = 4;
    double y0[ny];
    hr2_pmap_similarity_encode(&M0, y0);
    if (verbose) { hpedt_print_encoding(ny, y0); }
    demand(fabs(y0[0] - p.c[0]) <= eqtol, "{hr2_pmap_similarity_encode} failed (0)");
    demand(fabs(y0[1] - p.c[1]) <= eqtol, "{hr2_pmap_similarity_encode} failed (1)");
    demand(fabs(y0[2] - ang) <= eqtol, "{hr2_pmap_similarity_encode} failed (2)");
    /* Too complicated to test {y0[3]}. */

    for (sign_t sgn = +1; sgn >= -1; sgn -= 2)
      {
        if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
        if (verbose) { fprintf(stderr, "  testing with sgn = %+d\n", sgn); }
        
        hr2_pmap_t M = M0;
        hr2_pmap_set_sign(&M, sgn);
        if (verbose) { hpedt_print_map("hand-made", &M); }
        demand(sgn*r3x3_det(&(M.dir)) > 0, "{hr2_pmap_similarity_from_two_points} wrong handedness (2)");
        
        /* Test encoding: */
        double y[ny];
        hr2_pmap_similarity_encode(&M, y);
        if (verbose) { hpedt_print_encoding(ny, y); }
        for (int32_t jy = 0; jy < ny; jy++)
          { demand(fabs(y[jy] - y0[jy]) <= eqtol, "{hr2_pmap_similarity_encode} failed (3)"); }

        /* Test decoding: */
        hr2_pmap_t N = hr2_pmap_identity();
        hr2_pmap_similarity_decode(y, &N);
        if (verbose) { hpedt_print_map("decoded", &N); }
        hpedt_check_same_map(&M0, &N);

        if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
      }
      
    if (verbose) { fprintf(stderr, "< --- %s ---\n", __FUNCTION__); }
  }

void test_hr2_pmap_affine_encode__hr2_pmap_affine_decode(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "> --- %s ---\n", __FUNCTION__); }
   
    bool_t debug = TRUE;
    
    double eqtol = 1.0e-13; /* Tolerance for type matching. */
   
    /* Choose a positive affine frame {o,u,v}: */
    double xmax = 10.0, ymax = 20.0;
    r2_t o = (r2_t){{ dabrandom(-xmax,+xmax), dabrandom(-ymax,+ymax) }};
    /* Get vectors {u,v} with {+1} handedness: */
    r2_t u = (r2_t){{ 1.0 + dabrandom(-xmax,+xmax), dabrandom(-ymax,+ymax) }};
    r2_t v = (r2_t){{ dabrandom(-xmax,+xmax), 1.0 + dabrandom(-ymax,+ymax) }};
    double det = r2_det(&u, &v);
    if (det < 0) { r2_neg(&v, &v); }
    if (verbose) { r2_gen_print(stderr, &u, "%+18.10f", "    u = [ ", " ", " ]\n"); }
    if (verbose) { r2_gen_print(stderr, &v, "%+18.10f", "    v = [ ", " ", " ]\n"); }
    
    /* Generate an affine map {M0} with positive handedness: */
    r2_t p; r2_add(&o, &u, &p);
    r2_t q; r2_add(&o, &v, &q);
    hr2_pmap_t M0 = hr2_pmap_aff_from_three_points(&o, &p, &q);
    if (verbose) { hpedt_print_map("reference", &M0); }
    demand(hr2_pmap_is_type(&M0, hr2_pmap_type_AFFINE, +1, eqtol), "wrong map type or sign");
    demand(r3x3_det(&(M0.dir)) > 0, "{hr2_pmap_aff_from_three_points} wrong handedness (1)");

    /* Get the encoding {y0[0..3]} of {M0}: */
    int32_t ny = 6;
    double y0[ny];
    hr2_pmap_affine_encode(&M0, y0);
    if (verbose) { hpedt_print_encoding(ny, y0); }
    demand(fabs(y0[0] - o.c[0]) <= eqtol, "{hr2_pmap_affine_encode} failed (0)");
    demand(fabs(y0[1] - o.c[1]) <= eqtol, "{hr2_pmap_affine_encode} failed (1)");
    /* Too complicated to test {y0[2..5]}. */
  
    for (sign_t sgn = +1; sgn >= -1; sgn -= 2)
      {
        if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
        if (verbose) { fprintf(stderr, "  testing with sgn = %+d\n", sgn); }
        
        hr2_pmap_t M = hr2_pmap_aff_from_three_points(&o, &p, &q);
        hr2_pmap_set_sign(&M, sgn);
        if (verbose) { hpedt_print_map("hand-made", &M); }
        demand(sgn*r3x3_det(&(M.dir)) > 0, "{hr2_pmap_aff_from_three_points} wrong handedness (2)");
        
        /* Test encoding: */
        double y[ny];
        hr2_pmap_affine_encode(&M, y);
        if (verbose) { hpedt_print_encoding(ny, y); }
        for (int32_t jy = 0; jy < ny; jy++)
          { if (debug) { fprintf(stderr, "    diff[%d] = %24.16e\n", jy, y[jy] - y0[jy]); }
            demand(fabs(y[jy] - y0[jy]) <= eqtol, "{hr2_pmap_affine_encode} failed (3)");
          }

        /* Test decoding: */
        hr2_pmap_t N = hr2_pmap_identity();
        hr2_pmap_affine_decode(y, &N);
        if (verbose) { hpedt_print_map("decoded", &N); }
        hpedt_check_same_map(&M0, &N);

        if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
      }

    if (verbose) { fprintf(stderr, "< --- %s ---\n", __FUNCTION__); }
  }

void test_hr2_pmap_generic_encode__hr2_pmap_generic_decode(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "> --- %s ---\n", __FUNCTION__); }
   
    double eqtol = 1.0e-13; /* Tolerance for type matching. */
   
    /* Generate a projective map {M0} with positive handedness: */
    hr2_pmap_t M0;
    M0.dir = hpedt_throw_r3x3(10.0);
    double det = r3x3_det(&(M0.dir));
    if (det < 0) { for (int32_t j = 0; j < 3; j++) { M0.dir.c[2][j] = - M0.dir.c[2][j]; } }
    r3x3_inv(&(M0.dir), &(M0.inv));
    if (verbose) { hpedt_print_map("reference", &M0); }
    demand(hr2_pmap_is_type(&M0, hr2_pmap_type_GENERIC, +1, eqtol), "wrong map type or sign");
    
    /* Get the encoding {y0[0..8]} of {M0}: */
    int32_t ny = 9;
    double y0[ny];
    hr2_pmap_generic_encode(&M0, y0);
    if (verbose) { hpedt_print_encoding(ny, y0); }
    int32_t kv = 0;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++) 
          { double Aij = M0.dir.c[i][j];
            demand(fabs(y0[kv] - Aij) <= eqtol, "{hr2_pmap_generic_encode} failed (1)");
            kv++;
          }
      }
 
    for (sign_t sgn = +1; sgn >= -1; sgn -= 2)
      {
        if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
        if (verbose) { fprintf(stderr, "  testing with sgn = %+d\n", sgn); }
        
        hr2_pmap_t M = M0;
        hr2_pmap_set_sign(&M, sgn);
        if (verbose) { hpedt_print_map("hand-made", &M); }
        demand(sgn*r3x3_det(&(M.dir)) > 0, "generic map: wrong handedness (2)");

        /* Test encoding: */
        double y[ny]; /* Encoded elements. */
        hr2_pmap_generic_encode(&M, y);
        if (verbose) { hpedt_print_encoding(ny, y); }
        for (int32_t jy = 0; jy < ny; jy++)
          { demand(fabs(y[jy] - y0[jy]) <= eqtol, "{hr2_pmap_generic_encode} failed (3)"); }

        /* Test decoding: */
        hr2_pmap_t N = hr2_pmap_identity();
        hr2_pmap_generic_decode(y, &N);
        if (verbose) { hpedt_print_map("decoded", &N); }
        hpedt_check_same_map(&M0, &N);

        if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
      }

    if (verbose) { fprintf(stderr, "< --- %s ---\n", __FUNCTION__); }
  }

void test_hr2_pmap_encode__hr2_pmap_decode(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "> --- %s ---\n", __FUNCTION__); }
   
    double eqtol = 1.0e-13; /* Tolerance for type matching. */
   
    for (hr2_pmap_type_t type = hr2_pmap_type_FIRST; type <= hr2_pmap_type_LAST; type++)
      { 
        char *xtype = hr2_pmap_type_to_string(type);
        if (verbose) { fprintf(stderr, "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"); }
        
        /* Set {M0} to a random map of type {type} with {+1} handedness: */
        hr2_pmap_t M0 = hr2_pmap_throw_type(type, +1);
        if (verbose) { hpedt_print_map("reference", &M0); }
        
        /* Get the encoding {y0[0..ny-1]} of {M0}: */
        int32_t ny = hr2_pmap_encode_num_parameters(type);
        double y0[ny]; /* Encoded elements. */
        hr2_pmap_encode(&M0, type, ny, y0);
        if (verbose) { hpedt_print_encoding(ny, y0); }

        for (int32_t sgn = +1; sgn >= -1; sgn -= 2)
          { 
            if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
            if (verbose) { fprintf(stderr, "  testing with type = %s sgn = %+d\n", xtype, sgn); }

            /* Make from {M0} a map {M} with handedness {sgn}: */
            hr2_pmap_t M = M0;
            hr2_pmap_set_sign(&M, sgn);
            if (verbose) { hpedt_print_map("hand-made", &M); }

            /* Test encoding: */
            double y[ny]; /* Encoded elements. */
            hr2_pmap_encode(&M, type, ny, y);
            if (verbose) { hpedt_print_encoding(ny, y); }
            for (int32_t ky = 0; ky < ny; ky++)
              { demand(fabs(y[ky] - y0[ky]) <= eqtol, "{hr2_pmap_encode} failed"); }

            /* Test decoding: */
            hr2_pmap_t N = hr2_pmap_identity();
            hr2_pmap_decode(ny, y, type, sgn, &N);
            if (verbose) { hpedt_print_map("decoded", &N); }
            hpedt_check_same_map(&M, &N);
            
            if (verbose) { fprintf(stderr, "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); }
          }
          
        if (verbose) { fprintf(stderr, "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"); }
      }
      
    if (verbose) { fprintf(stderr, "< --- %s ---\n", __FUNCTION__); }
  }

r3x3_t hpedt_throw_r3x3(double mag)
  { r3x3_t A;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { A.c[i][j] = mag*(2*drandom() - 1.0); }
      }
    return A;
  }

void hpedt_print_map(char *name, hr2_pmap_t *M)
  {
    fprintf(stderr, "  %s = \n", name);
    hr2_pmap_gen_print(stderr, M, "%12.7f", "    ", "[ ","  "," ]\n    ", "[ "," "," ]", "\n");
  }

void hpedt_print_encoding(int32_t ny, double y[])
  { 
    fprintf(stderr, "  encoding parameters:\n");
    fprintf(stderr, "    ");
    for (int32_t kv = 0; kv < ny; kv++)
      { fprintf(stderr, " %+12.8f", y[kv]); }
    fprintf(stderr, "\n");
  }

void hpedt_check_same_map(hr2_pmap_t *M, hr2_pmap_t *N)
  { for (int32_t dir = 0; dir <= 1; dir++)
      { r3x3_t *A = (dir == 0 ? &(M->dir) : &(M->inv));
        double Am = r3x3_norm(A);
        demand(Am > 1.0e-100, "Bad {M}");
        
        r3x3_t *B = (dir == 0 ? &(N->dir) : &(N->inv));
        double Bm = r3x3_norm(B);
        demand(Am > 1.0e-100, "Bad {N}");
        
        double tol = 1.0e-12;
        for (int32_t i = 0; i < 3; i++)
          for (int32_t j = 0; j < 3; j++) 
            { double Aij = A->c[i][j]/Am;
              double Bij = B->c[i][j]/Bm;
              double err = fabs(Aij - Bij);
              if (err > tol)
                { fprintf(stderr, "    A[%d][%d] = %.16f  B[%d][%d] = %.16f ", i, j, Aij, i, j, Bij);
                  fprintf(stderr, "  err = %24.16f\n", err);
                  demand(FALSE, "** maps differ -- encode/decode failed\n");
                }
            }
      }
  }

void test_hr2_pmap_throw_type__hr2_pmap_is_type(bool_t verbose)
  { if (verbose) { fprintf(stderr, "--- hr2_pmap_is_type ---\n"); }
    
    double tol = 1.0e-14;
    hr2_pmap_t M; 
    for (hr2_pmap_type_t type = hr2_pmap_type_FIRST; type <= hr2_pmap_type_LAST; type++)
      { 
        for (sign_t sgn = +1; sgn >= -1; sgn -= 2)
          { M = hr2_pmap_throw_type(type, sgn);
            demand(hr2_pmap_is_type(&M, type, sgn, tol), "{hr2_pmap_is_type} failed (1)");

            auto bool_t ok(hr2_pmap_t *A);

            /* Check whether perturbing {M} at specific elements falsifies the test: */
            double eps = 50*tol;
            r3x3_t P; /* Perturbation matrix */
            switch(type)
              { case hr2_pmap_type_IDENTITY:
                  P = (r3x3_t){{ { eps, eps, eps }, { eps, eps, eps }, { eps, eps, eps } }}; break;
                case hr2_pmap_type_TRANSLATION:
                  P = (r3x3_t){{ { 0.0, 0.0, 0.0 }, { eps, eps, eps }, { eps, eps, eps } }}; break;
                case hr2_pmap_type_CONGRUENCE:
                  P = (r3x3_t){{ { eps, 0.0, 0.0 }, { eps, eps, eps }, { eps, eps, eps } }}; break;
                case hr2_pmap_type_SIMILARITY:
                  P = (r3x3_t){{ { 0.0, 0.0, 0.0 }, { eps, eps, eps }, { eps, eps, eps } }}; break;
                case hr2_pmap_type_AFFINE:
                  P = (r3x3_t){{ { 0.0, 0.0, 0.0 }, { eps, 0.0, 0.0 }, { eps, 0.0, 0.0 } }}; break;
                case hr2_pmap_type_GENERIC:
                  P = (r3x3_t){{ { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } }}; break;
                case hr2_pmap_type_NONE:
                  demand(FALSE, "invalid projective map type");
                default:
                  demand(FALSE, "unimplemented map type");
              }
            hr2_test_perturbed_pmap(&M, &ok, &P, "hr2_pmap_is_type");
            return;

            bool_t ok(hr2_pmap_t *A)
              { return hr2_pmap_is_type(A, type, sgn, tol); }
          }
      }
  }
