#define PROG_NAME "hr2_pmap_opt_test"
#define PROG_DESC "test of {hr2_pmap_opt.h} and related modules"
#define PROG_VERS "1.0"

/* Last edited on 2023-10-21 16:00:44 by stolfi */ 
/* Created on 2020-07-11 by J. Stolfi, UNICAMP */
/* Based on {test_align.c} by J. Stolfi, UNICAMP */

#define hr2_pmap_opt_test_COPYRIGHT \
  "Copyright © 2020  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <wt_table.h>
#include <ix.h>
#include <r2.h>
#include <hr2.h>
#include <rn.h>
#include <r3x3.h>
#include <bool.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>
#include <assert.h>

#include <hr2_pmap_opt.h>
#include <hr2_pmap_opt_translation.h>
#include <hr2_pmap_opt_congruence.h>
#include <hr2_pmap_opt_similarity.h>
#include <hr2_pmap_opt_affine.h>
#include <hr2_pmap_opt_generic.h>

typedef struct hpot_options_t 
  { char *outPrefix;      /* Prefix for output files. */
    hr2_pmat_type_t type; /* Map type. */
    hr2_pmap_t optimum;   /* Expected optimum. */
    r3x3_t deviation;     /* Max elementwise adjustment in direct array of map. */
    int32_t ns;           /* Number of evaluations for plotting. */
  } hpot_options_t;
  /* Command line parameters.
    
    The matrix {deviation} is the search region radius: it determines the
    entries of the projective map that are to be adjusted, and the
    amount of adjustment. Namely, entry {A.dir[i][j]} can be ajusted by
    {±deviation[i][j]}. In particular, element {A.dir[i][j]} will be fixed at
    {optimum.dir[i][j]} if {deviation[i][j]} is zero. The elements of {deviation} must be
    non-negative. Because of the homoegneity property, at least one of
    the elements of {deviation} must be zero. */

void hpot_test_enc_dec_translation(void);
void hpot_test_enc_dec_congruence(void);
void hpot_test_enc_dec_similarity(void);
void hpot_test_enc_dec_affine(void);
void hpot_test_enc_dec_generic(void);
  /* Tests {hr2_pmap_opt_translation_encode} and
    {hr2_pmap_opt_translation_encode}. */

void hpot_normalize_matrix(r3x3_t *A);
  /* Tries to scale all elements of {A} so that element {A[0][0]} becomes 1.
    If that element is too small, scales so that the sum of squares of elements is 1.
    if all elements are too small, leaves the matrix unchanged. */

void hpot_print_encoding(int32_t nv, double rad[], double y[]);
  /* Prints {y[0..nv-1]} on one line, then {rad[0..nv-1]} on the next line,
    all preceded by a line saying "encoding and scaling factors". */

void hpot_check_decoded_map(hr2_pmap_t *M, hr2_pmap_t *N);
  /* Compares the original map {M} with the encoded and decoded map {N}. 
    Assumes that the direct and inverse matrices are normalized
    as per {hpot_normalize_matrix}. */

void hpot_one(hpot_options_t *o, char *method, bool_t verbose);
  /* Tests affine map adjustment algorithm {method} ("quad" etc). */

void hpot_choose_initial_guess(hr2_pmap_t *Aopt, r3x3_t *R, hr2_pmap_t *Aini);
  /* Stores into {*Aini} the initail guess for the affine map
    to be adjusted, chosen at random in the box {Aopt±(R/2)}.  
    At least one of them must be zero. */

r3x3_t hpot_throw_elem_var_matrix(void);
  /* Generate a matrix {R} with non-negative elements,
    at most 8 of them nonzero. */

void hpot_get_var_elems
  ( r3x3_t *M, 
    r3x3_t *R, 
    double *Mp[],
    double Re[],
    int32_t *nvP
  );
  /* Identifies the elements of {*M} that are to be adjusted 
    (that is, whose corresponding element of {*R} is nonzero.
    Returns the number {nv} of such elements (a number in {0..8}) in {*nvP},
    the addresses of those elemens of {*M} in {Mp[0..nv-1]},
    and the values of the corresponding elements of {*R}
    in {Re[0..nv-1]}.  These vectors should have at leat 8 elements.
    The procedure fails if all 9 elements of {*R} are nonzero.*/

void hpot_print_map(char *name, hr2_pmap_t *M);
  /* Prints the map {M} (direct and inverse) to {stderr},
    prefixed by a line"  {name} = ". */

void hpot_show_disp
  ( hr2_pmap_t *A, 
    char *Bname, 
    hr2_pmap_t *B,
    r3x3_t *R
  );
  /* Shows the discrepancy beweeen maps {A} and {B}.
    
    The arrays are compared with {hr2_pmap_diff_sqr(A,B)} and with
    {hpot_diff_rel_sqr(D,R)} where {D=A.dir-B.dir}. If {print} is true
    also prints {D}. If the maps {A} and {B} are affine, they are also
    compared with {hr2_pmap_aff_discr_sqr(A,B)}. */

void hpot_debug_map
  ( char *label,
    bool_t print, 
    hr2_pmap_t *A, 
    hr2_pmap_t *Aini, 
    hr2_pmap_t *Aopt, 
    r3x3_t *R,
    double F2A
  );
  /* Prints the affine map {*A} and the corresponding raw goal function value {F2A}.
    If {Aini} is not {NULL}, prints also the difference between 
    {*A} and the corresponding {*Aini}.  Ditto for the difference to {Aopt},
    if {Aopt} is not {NULL}. */

void hpot_plot_goal
  ( char *outPrefix,
    char *method,
    hr2_pmap_opt_func_t *F2, 
    hr2_pmap_t *A,
    r3x3_t *U,
    r3x3_t *V,
    int32_t ns
  );
  /* Writes a file "{outPrefix}_{method}.dat" with a random 2D slice of the goal
    function {F2(B)}, where {B} ranges in the neighborhood of {A}, in
    the plane definde by the displacements {±U} and {±V}.
    
    Specifically, generates the affine map as {B = A + u*U + v*V[k]},
    where {u} is {tu/ns} for {tu} in {-ns..+ns}, and {v} is {tv/ns} for
    {tv} in {-ns..+ns}. Plots {F2(B)} as a function of {u} and {v}. */

void hpot_choose_plot_directions
  ( r3x3_t *R, 
    r3x3_t *U, 
    r3x3_t *V
  );
  /* Chooses two random matrices {U,V} that are orthogonal and of equal
    length when viewed as vectors of {\RR^9}. They have
    non-zero elements only where {R} itself is non-zero. Their norms are
    adjusted so that the rectangle {{a*U + b*V : a,b \in [-1_+1]}} lies
    within the box {±R}. */
  
double hpot_mismatch_sqr_1(hr2_pmap_t *A, hr2_pmap_t *B);
  /* A squared mismatch function between {*A} and {*B},
    defined as {hr2_pmap_diff_sqr(A,B)}. */
  
double hpot_mismatch_sqr_2(hr2_pmap_t *A, hr2_pmap_t *B);
  /* A squared mismatch function between {*A} and {*B},
    defined a positive definite quadratic funtion of 
    the six elements of {A^{-1} B - I} where {I} is the 
    identity.  This metric is used instead of the simple squared
    Euclidean distance (sum of squares of element differences)
    because the latter would make elementwise optimzation too easy. */
  
double hpot_diff_rel_sqr(r3x3_t *D, r3x3_t *R);
  /* Computes the sum of squared elements of {D}, each divided by the
    corresponding element of {R}.  The elements that are zero in {R} 
    must be zero in {D} too, and are ignored. */

hpot_options_t *hpot_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */
  
void hpot_parse_next_r3x3(argparser_t *pp, r3x3_t *R);
  /* Parses a 3x3 matrix from the command line, as 9 numbers 
    {R[0][0] R[0][1] R[0][2]  R[1][0] ... R[2][2]}. */

int32_t main(int32_t argn, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srand(4615*417);
    bool_t verbose = TRUE;
    
    hpot_options_t *o = hpot_parse_options(argc, argv);
    
    hpot_test_enc_dec_translation();
    hpot_test_enc_dec_congruence();
    hpot_test_enc_dec_similarity();
    hpot_test_enc_dec_affine();
    hpot_test_enc_dec_generic();

    hpot_one(o, "quad", verbose);

    return 0;
  }
  
void hpot_test_enc_dec_translation(void)
  { 
    fprintf(stderr, "> --- %s -----------------------------------------------------\n", __FUNCTION__);
   
    /* Choose ranges and the test matrix: */
    double xmax = 10.0, ymax = 20.0;
    r2_t disp = (r2_t){{ dabrandom(-xmax,+xmax), dabrandom(-ymax,+ymax) }};
    hr2_pmap_t M = hr2_pmap_translation(&disp);
    hpot_normalize_matrix(&(M.inv));  /* Assumes {M.dir} is normalized. */
    hpot_print_map("original", &M);

    /* Test encoding: */
    int32_t nv = 2;
    double rad[nv], y[nv];
    rad[0] = xmax;
    rad[1] = ymax;
    hr2_pmap_opt_translation_encode(&M, rad, y);
    hpot_print_encoding(nv, rad, y);
    demand(fabs(y[0] - disp.c[0]/rad[0]) < 1.0e-14, "encode failed (0)");
    demand(fabs(y[1] - disp.c[1]/rad[1]) < 1.0e-14, "encode failed (1)");

    /* Test decoding: */
    hr2_pmap_t N = hr2_pmap_identity();
    hr2_pmap_opt_translation_decode(y, rad, &N);
    hpot_normalize_matrix(&(N.inv));  /* Assumes {N.dir} is normalized. */
    hpot_print_map("decoded", &N);
    hpot_check_decoded_map(&M, &N);
    
    fprintf(stderr, "< --- %s -----------------------------------------------------\n", __FUNCTION__);
  }

void hpot_test_enc_dec_congruence(void)
  { 
    fprintf(stderr, "> --- %s -----------------------------------------------------\n", __FUNCTION__);
   
    /* Choose ranges and the test matrix: */
    double xmax = 10.0, ymax = 20.0;
    double angmax = 0.5; /* Better be less than {pi/2}. */
    r2_t disp = (r2_t){{ dabrandom(-xmax,+xmax), dabrandom(-ymax,+ymax) }};
    double ang = dabrandom(-angmax, +angmax);
    r2_t u = (r2_t){{ cos(ang), sin(ang) }};
    hr2_pmap_t M = hr2_pmap_congruence_from_point_and_dir(&disp, &u, FALSE);
    hpot_normalize_matrix(&(M.inv));
    hpot_print_map("original", &M);

    /* Test encoding: */
    int32_t nv = 3;
    double rad[nv], y[nv];
    rad[0] = xmax;
    rad[1] = ymax;
    rad[2] = angmax;
    hr2_pmap_opt_congruence_encode(&M, rad, y);
    hpot_print_encoding(nv, rad, y);
    demand(fabs(y[0] - disp.c[0]/rad[0]) < 1.0e-14, "encode failed (0)");
    demand(fabs(y[1] - disp.c[1]/rad[1]) < 1.0e-14, "encode failed (1)");
    demand(fabs(y[2] - ang/rad[2]) < 1.0e-14, "encode failed (2)");

    /* Test decoding: */
    hr2_pmap_t N = hr2_pmap_identity();
    hr2_pmap_opt_congruence_decode(y, rad, &N);
    hpot_normalize_matrix(&(N.inv));  /* Assumes {N.dir} is normalized. */
    hpot_print_map("decoded", &N);
    hpot_check_decoded_map(&M, &N);

     fprintf(stderr, "< --- %s -----------------------------------------------------\n", __FUNCTION__);
  }

void hpot_test_enc_dec_similarity(void)
  { 
    fprintf(stderr, "> --- %s -----------------------------------------------------\n", __FUNCTION__);
   
    /* Choose ranges and the test matrix: */
    double xmax = 10.0, ymax = 20.0;
    r2_t p = (r2_t){{ dabrandom(-xmax,+xmax), dabrandom(-ymax,+ymax) }};
    r2_t u = (r2_t){{ 1.0 + dabrandom(-xmax,+xmax), dabrandom(-ymax,+ymax) }};
    r2_t q; r2_add(&p, &u, &q);
    hr2_pmap_t M = hr2_pmap_similarity_from_two_points(&p, &q, FALSE);
    hpot_normalize_matrix(&(M.inv));
    hpot_print_map("original", &M);

    /* Test encoding: */
    int32_t nv = 4;
    double rad[nv], y[nv];
    rad[0] = xmax;
    rad[1] = ymax;
    rad[2] = xmax; 
    rad[3] = ymax; 
    hr2_pmap_opt_similarity_encode(&M, rad, y);
    hpot_print_encoding(nv, rad, y);
    demand(fabs(y[0] - p.c[0]/rad[0]) < 1.0e-14, "encode failed (0)");
    demand(fabs(y[1] - p.c[1]/rad[1]) < 1.0e-14, "encode failed (1)");
    demand(fabs(y[2] - (u.c[0]-1.0)/rad[2]) < 1.0e-14, "encode failed (2)");
    demand(fabs(y[3] - u.c[1]/rad[3]) < 1.0e-14, "encode failed (3)");

    /* Test decoding: */
    hr2_pmap_t N = hr2_pmap_identity();
    hr2_pmap_opt_similarity_decode(y, rad, &N);
    hpot_normalize_matrix(&(N.inv));  /* Assumes {N.dir} is normalized. */
    hpot_print_map("decoded", &N);
    hpot_check_decoded_map(&M, &N);

    fprintf(stderr, "< --- %s -----------------------------------------------------\n", __FUNCTION__);
  }

void hpot_test_enc_dec_affine(void)
  { 
    fprintf(stderr, "> --- %s -----------------------------------------------------\n", __FUNCTION__);
   
    /* Choose ranges and the test matrix: */
    double xmax = 10.0, ymax = 20.0;
    r2_t o = (r2_t){{ dabrandom(-xmax,+xmax), dabrandom(-ymax,+ymax) }};
    r2_t u = (r2_t){{ 1.0 + dabrandom(-xmax,+xmax), dabrandom(-ymax,+ymax) }};
    r2_t v = (r2_t){{ dabrandom(-xmax,+xmax), 1.0 + dabrandom(-ymax,+ymax) }};
    r2_t p; r2_add(&o, &u, &p);
    r2_t q; r2_add(&o, &v, &q);
    hr2_pmap_t M = hr2_pmap_aff_from_three_points(&o, &p, &q);
    hpot_normalize_matrix(&(M.inv));
    hpot_print_map("original", &M);

    /* Test encoding: */
    int32_t nv = 6;
    double rad[nv], y[nv];
    rad[0] = xmax;
    rad[1] = ymax;
    rad[2] = xmax; 
    rad[3] = ymax; 
    rad[4] = xmax; 
    rad[5] = ymax; 
    hr2_pmap_opt_affine_encode(&M, rad, y);
    hpot_print_encoding(nv, rad, y);
    demand(fabs(y[0] - o.c[0]/rad[0]) < 1.0e-14, "encode failed (0)");
    demand(fabs(y[1] - o.c[1]/rad[1]) < 1.0e-14, "encode failed (1)");
    demand(fabs(y[2] - (u.c[0] - 1)/rad[2]) < 1.0e-14, "encode failed (2)");
    demand(fabs(y[3] - u.c[1]/rad[3]) < 1.0e-14, "encode failed (3)");
    demand(fabs(y[4] - v.c[0]/rad[4]) < 1.0e-14, "encode failed (4)");
    demand(fabs(y[5] - (v.c[1] - 1)/rad[5]) < 1.0e-14, "encode failed (5)");

    /* Test decoding: */
    hr2_pmap_t N = hr2_pmap_identity();
    hr2_pmap_opt_affine_decode(y, rad, &N);
    hpot_normalize_matrix(&(N.inv));  /* Assumes {N.dir} is normalized. */
    hpot_print_map("decoded", &N);
    hpot_check_decoded_map(&M, &N);

    fprintf(stderr, "< --- %s -----------------------------------------------------\n", __FUNCTION__);
  }

void hpot_test_enc_dec_generic(void)
  { 
    fprintf(stderr, "> --- %s -----------------------------------------------------\n", __FUNCTION__);
   
    /* Choose element ranges: */
    r3x3_t R = hpot_throw_elem_var_matrix();
    fprintf(stderr, "template matrix = \n");
    r3x3_gen_print(stderr, &R, "%12.7f", "    [ ","\n      "," ]\n", "[ "," "," ]");
    
    /* Choose the test matrix: */
    hr2_pmap_t M = hr2_pmap_identity();
    double *Mp[9]; /* Pointers to variable elems of {M.dir}. */
    double rad[9]; /* The non-zero elems of {R}. */
    int32_t nv; /* Count of nonzero elems in {R}. */
    hpot_get_var_elems(&(M.dir), &R, Mp, rad, &nv);
    assert(nv <= 8); /* Shoul be ensured by {hpot_throw_elem_var_matrix}. */
    for (int32_t kv = 0; kv < nv; kv++)
      { double *Mij = Mp[kv];
        (*Mij) = dabrandom(-rad[kv], +rad[kv]);
      }
    r3x3_inv(&(M.dir), &(M.inv));
    hpot_normalize_matrix(&(M.inv));
    hpot_print_map("original", &M);

    /* Test encoding: */
    double y[nv];
    hr2_pmap_opt_generic_encode(&M, &R, y);
    hpot_print_encoding(nv, rad, y);
    int32_t kv = 0;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++) 
          { double Rij = R.c[i][j];
            if (Rij != 0)
              { double Dij = M.dir.c[i][j] - (i == j ? 1 : 0);
                demand(fabs(y[kv] - Dij/rad[kv]) < 1.0e-14, "encode failed");
                kv++;
              }
          }
      }

    /* Test decoding: */
    hr2_pmap_t N = hr2_pmap_identity();
    hr2_pmap_opt_generic_decode(y, &R, &N);
    hpot_normalize_matrix(&(N.inv));  /* Assumes {N.dir} is normalized. */
    hpot_print_map("decoded", &N);
    hpot_check_decoded_map(&M, &N);

    fprintf(stderr, "< --- %s -----------------------------------------------------\n", __FUNCTION__);
  }

void hpot_normalize_matrix(r3x3_t *A)
  { double norm = r3x3_norm(A);
    double head = fabs(A->c[0][0]);
    double scale = (head > 0.01*norm ? head : norm); 
    if ((scale == 1.0) || (scale < 1.0e-200)) { return; }
    for(int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++) 
          { A->c[i][j] /= scale; }
      }
  }

void hpot_check_decoded_map(hr2_pmap_t *M, hr2_pmap_t *N)
  { for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++) 
          { demand(fabs(N->dir.c[i][j] - M->dir.c[i][j]) < 1.0e-14, "decode failed (dir)");
            demand(fabs(N->inv.c[i][j] - M->inv.c[i][j]) < 1.0e-14, "decode failed (inv)");
          }
      }
  }

void hpot_print_encoding(int32_t nv, double rad[], double y[])
  { 
    fprintf(stderr, "  encoding and scaling factors:\n");
    fprintf(stderr, "    ");
    for (int32_t kv = 0; kv < nv; kv++)
      { fprintf(stderr, " %+12.8f", y[kv]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "    ");
    for (int32_t kv = 0; kv < nv; kv++)
      { fprintf(stderr, " %+12.8f", rad[kv]); }
    fprintf(stderr, "\n");
  }

void hpot_one(hpot_options_t *o, char *method, bool_t verbose)
  { 
    fprintf(stderr, "> --- %s -----------------------------------------------------\n", __FUNCTION__);
    fprintf(stderr, "testing with method %s\n", method);
    
    bool_t debug_maps = FALSE; /* TRUE to print every probe point. */
    
    auto double F2_mismatch(hr2_pmap_t *A);
      /* A goal function that returns the squared mismatch between a map
        {*A} and the Aopt {*Aopt}. Also prints the map if
        {debug_maps} is true. */
    
    hr2_pmap_t *Aopt = &(o->optimum);
    r3x3_t *R = &(o->deviation);
    r3x3_gen_print(stderr, R, "%12.7f", "R = [ ","\n      "," ]\n", "[ "," "," ]");
    
    /* Choose the initial guess {Aini}: */
    hr2_pmap_t Aini;
    hpot_choose_initial_guess(Aopt, R, &Aini);
       
    /* Choose two orthogonal perturbation vectors in {\RR^6} for plotting: */
    r3x3_t U, V;
    hpot_choose_plot_directions(R, &U, &V); 
    
    /* Print raw function value and bias term for Aopt map: */
    double F2opt = F2_mismatch(Aopt);
    hpot_debug_map("actual optimum", TRUE, Aopt, NULL, NULL, R, F2opt);
    
    /* Print raw function value and bias term for initial map: */
    double F2ini = F2_mismatch(&Aini);
    hpot_debug_map("initial guess ", TRUE, &Aini, NULL, Aopt, R, F2ini);
    
    /* Plot the goal function in the neighborhood of the initial guess: */ 
    fprintf(stderr, "Plotting goal function around initial point...\n");
    hpot_plot_goal(o->outPrefix, method, F2_mismatch, &Aini, &U, &V, o->ns);
    
    /* Call the optimizer: */
    fprintf(stderr, "optimizing...\n");
    debug_maps = TRUE;
    hr2_pmap_t Asol = Aini;  /* Computed affine map. */
    double F2sol;   /* Goal function at {Asol}. */
    if (strcmp(method, "quad") == 0)
      { double max_mismatch = 0.02;
        int32_t maxIter = 5;
        double max_mod = 0.25;
        hr2_pmap_opt_aff_quadratic(F2_mismatch, o->type, maxIter, max_mismatch, max_mod, &Asol, &F2sol);
      }
    else
      { demand(FALSE, "invalid method"); }
    
    /* Print raw function value and bias term for computed optimum: */
    F2sol = F2_mismatch(&Asol);
    hpot_debug_map("computed optimum", TRUE, &Asol, &Aini, Aopt, R, F2sol);

    fprintf(stderr, "< --- %s -----------------------------------------------------\n", __FUNCTION__);
    return;
    
    double F2_mismatch(hr2_pmap_t *A)
      { double d2 = hpot_mismatch_sqr_1(Aopt, A);
        if (debug_maps) 
          { hpot_debug_map("probe map", verbose, A, &Aini, Aopt, R, d2); }
        return d2;
      }
  }

double hpot_mismatch_sqr_1(hr2_pmap_t *A, hr2_pmap_t *B)
  { 
    double diff2 = hr2_pmap_diff_sqr(A, B);
    return diff2;
  }

double hpot_mismatch_sqr_2(hr2_pmap_t *A, hr2_pmap_t *B)
  { 
    hr2_pmap_t ABdif = hr2_pmap_inv_compose(A, B); /* The map {A^{-1} B}. */
    /* Compare the matrix of {ABdif} with the identity map, using a non-trivial metric: */
    double d2 = 0;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { double d = ABdif.dir.c[i][j] - (i == j ? 1.0 : 0.0);
            d2 += (1 + 0.1*i + 0.3*j)*d*d;
          }
      }
    return d2;
  }
  
void hpot_choose_initial_guess(hr2_pmap_t *Aopt, r3x3_t *R, hr2_pmap_t *Aini)
  {
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { double frac = dabrandom(-0.50, +0.50);
            Aini->dir.c[i][j] = Aopt->dir.c[i][j] + frac*R->c[i][j];
          }
      }
      r3x3_inv(&(Aini->dir), &(Aini->inv));
  }
    
void hpot_plot_goal
  ( char *outPrefix,
    char *method,
    hr2_pmap_opt_func_t *F2, 
    hr2_pmap_t *A,
    r3x3_t *U,
    r3x3_t *V,
    int32_t ns
  )
  {
    /* Show the {U,V} vectors: */
    r3x3_gen_print(stderr, &(A->dir), "%12.7f", "A = [ ","\n      "," ]\n", "[ "," "," ]");
    r3x3_gen_print(stderr, U,         "%12.7f", "U = [ ","\n      "," ]\n", "[ "," "," ]");
    r3x3_gen_print(stderr, V,         "%12.7f", "V = [ ","\n      "," ]\n", "[ "," "," ]");
    
    /* Sweep the {A,U,V} plane and Plot: */
    char *fname = NULL;
    asprintf(&fname, "%s_%s.dat", outPrefix, method);
    FILE *wr = open_write(fname, TRUE);
    
    hr2_pmap_t B; /* Probe map. */
    fprintf(stderr, "\n");
    for (int32_t iu = -ns; iu <= +ns; iu++)
      { double du = ((double)iu)/((double)ns);
        fprintf(stderr, ".");
        for (int32_t iv = -ns; iv <= +ns; iv++)
          { double dv = ((double)iv)/((double)ns);
            /* Compute the probe map {B}. */
            r3x3_mix(du, U, dv, V, &(B.dir)); 
            r3x3_add(&(A->dir), &(B.dir), &(B.dir));
            r3x3_inv(&(B.dir), &(B.inv));
            /* Evaluate the function and Plot: */
            double F2B = F2(&B);
            fprintf(wr, "%+9.6f %+9.6f  %12.6f\n", du, dv, F2B); 
          }
        /* Blank line between scanlines, for {gnuplot}: */
        fprintf(wr, "\n"); 
      }

    fprintf(stderr, "\n");
    fclose(wr);
    free(fname);
  }

void hpot_choose_plot_directions(r3x3_t *R, r3x3_t *U, r3x3_t *V)
  { 
    /* Identify the elements of {R} that are not zero, get the addresses {Up,Vp}: */
    int32_t ne = 0;
    double Re[9]; /* Nonzero elements of {R} are  {Re[0..ne-1]}. */
    double *Up[9], *Vp[9]; /* Pointers of elements of {U} and {V} where {R} is nonzero. */
    
    hpot_get_var_elems(U, R, Up, Re, &ne);
    hpot_get_var_elems(V, R, Vp, Re, &ne);
      
    /* Throw a random unit vector {ue[0..ne-1]}:*/
    double ue[ne];
    rn_throw_dir(ne, ue);
    
    /* Throw a unit vector {ve[0..ne-1]} orthogonal to {ue}: */
    double ve[ne];
    do
      { rn_throw_dir(ne, ve);
        rn_decomp(ne, ve, ue, NULL, ve);
      }
    while (rn_norm_sqr(ne, ve) < 1.0e-4);
    double vnorm = rn_dir(ne, ve, ve);
    assert(vnorm > 0);
    
    /* Adjust the lengths to fit in {R}: */
    /* Find the max scale {s} that keeps {s*(±ue±ve)} withing {±R}: */
    double s = +INF;
    double te[ne];
    for (int32_t du = -1; du <= +1; du += 2)
      { for (int32_t dv = -1; dv <= +1; dv += 2)
          { rn_mix(ne, (double)du, ue, (double)dv, ve, te);
            for (int32_t ie = 0; ie < ne; ie++)
              { double tei = fabs(te[ie]);
                double Rei = fabs(Re[ie]);
                assert (Rei > 0);
                if (tei > 1.0e-15*Rei)
                  { double si = Rei/tei;
                    if (si < s) { s = si; }
                  }
               }
          }
      }
    assert(isfinite(s));
    assert(s > 0);
   
    /* Insert {ue,ve} scaled by {s} into {U,V}: */
    r3x3_zero(U);
    r3x3_zero(V);
    for (int32_t ie = 0; ie < ne; ie++)
      { (*(Up[ie])) = s*ue[ie];
        (*(Vp[ie])) = s*ve[ie];
      }
  }

r3x3_t hpot_throw_elem_var_matrix(void)
  { 
    r3x3_t R;
    double R_max = 9.999;
    int32_t mv = int32_abrandom(1, 8); /* Target number of non-zero elements. */
    int32_t nv = 0; /* Elements of {R} set to non-zero. */
    int32_t ne = 0; /* Elements of {R} enumerated so far. */
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++) 
          { double prob = ((double)(mv - nv))/(9.0 - ne);
            if (drandom() < prob)
              { assert(nv < 8);
                R.c[i][j] = dabrandom(0.01, R_max);
                nv++;
              }
            else
              { R.c[i][j] = 0.0; }
            ne++;
          }
      }
    assert(ne == 9);
    assert(nv <= 8);
    return R;
  }

void hpot_get_var_elems
  ( r3x3_t *M, 
    r3x3_t *R, 
    double *Mp[],
    double Re[],
    int32_t *nvP
  )
  { 
    int32_t nv = 0;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { double Rij = R->c[i][j];
            demand(Rij >= 0, "elems of {R} must be non-negative");
            if (Rij != 0)
              { demand(nv < 8, "too many nonzero elems in {R}, max 8");
                Mp[nv] = &(M->c[i][j]);
                Re[nv] = Rij;
                nv++;
              }
          }
      }
    (*nvP) = nv;
  }
  
double hpot_diff_rel_sqr(r3x3_t *D, r3x3_t *R)
  {
    double sum_r2 = 0;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { double Dij = D->c[i][j];
            double Rij = R->c[i][j];
            if (Rij == 0)
              { assert(Dij == 0); }
            else
              { double rij = Dij/Rij; 
                sum_r2 += rij*rij;
              }
          }
      }
    return sum_r2;
  }

void hpot_show_disp
  ( hr2_pmap_t *A, 
    char *Bname,
    hr2_pmap_t *B, 
    r3x3_t *R
  )
  {
    bool_t debug = TRUE;
    
    fprintf(stderr, "  comparing with %s:\n", Bname);

    /* Difference between {dir} and {inv} arrays, normalized: */
    double diff2 = hr2_pmap_diff_sqr(A, B);

    /* Difference between {dir} arrays relative to {R}: */
    r3x3_t D;
    r3x3_sub(&(A->dir), &(B->dir), &D);
    if (debug) { r3x3_gen_print(stderr, &D, "%12.7f", "  D = [ ","\n        "," ]\n", "[ "," "," ]"); }
    double drel2 = hpot_diff_rel_sqr(&D, R);
    fprintf(stderr, "  diff = %12.7f drel = %12.7f", sqrt(diff2), sqrt(drel2));
    
    if (hr2_pmap_is_affine(A) && hr2_pmap_is_affine(B))
      { /* Difference between {A} and {B} over the unit circle: */
        double mism2 = hr2_pmap_aff_discr_sqr(A, B);
        fprintf(stderr, " aff mism = %12.7f", sqrt(mism2));
      }
    fprintf(stderr, "\n");
  }
  
void hpot_print_map(char *name, hr2_pmap_t *M)
  {
    fprintf(stderr, "  %s = \n", name);
    hr2_pmap_gen_print(stderr, M, "%12.7f", "    ", "[ ","  "," ]\n    ", "[ "," "," ]", "\n");
  }

void hpot_debug_map
  ( char *label,
    bool_t print, 
    hr2_pmap_t *A,
    hr2_pmap_t *Aini,
    hr2_pmap_t *Aopt, 
    r3x3_t *R,
    double F2A
  )
  { 
    fprintf(stderr, "%s\n", label);
    if (print)
      { hr2_pmap_gen_print
          ( stderr, A, 
            "%12.7f", "  A = \n    ", 
            "[ ","  "," ]\n    ",
            "[ "," "," ]",
            "\n"
          );
      }

    if (Aini != NULL) { hpot_show_disp(A, "Aini", Aini, R); }
    if (Aopt != NULL) { hpot_show_disp(A, "optimum", Aopt, R); }

    fprintf(stderr, "  F2(A) = %22.10f\n", F2A);
    fprintf(stderr, "\n");
    return;
   }

hpot_options_t *hpot_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, "HELP!");
    argparser_set_info(pp, "KNOW THYSELF");
    argparser_process_help_info_options(pp);
    
    hpot_options_t *o = notnull(malloc(sizeof(hpot_options_t)), "no mem");

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_get_keyword(pp, "-optimum");
    hpot_parse_next_r3x3(pp, &(o->optimum.dir));  /* Target map (expected optimum). */
    r3x3_inv(&(o->optimum.dir), &(o->optimum.inv));

    argparser_get_keyword(pp, "-deviation");
    hpot_parse_next_r3x3(pp, &(o->deviation));     /* Max elementwise deviation. */
    
    argparser_get_keyword(pp, "-nPlot");
    o->ns = (int32_t)argparser_get_next_uint(pp, 1, 500); /* Number of evaluations for plotting. */

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }

void hpot_parse_next_r3x3(argparser_t *pp, r3x3_t *R)
  {
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { R->c[i][j] = argparser_get_next_double(pp, -100.0, +100.0); }
      }
  }

