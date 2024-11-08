#define PROG_NAME "test_hr2_pmap_aff_from_point_pairs"
#define PROG_DESC "test of {hr2_pmap_aff_from_point_pairs.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-08 11:19:45 by stolfi */ 
/* Created on 2020-07-11 by J. Stolfi, UNICAMP */
/* Based on {test_align.c} by J. Stolfi, UNICAMP */

#define test_hr2_pmap_aff_from_point_pairs_COPYRIGHT \
  "Copyright © 2020  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

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

#include <hr2_pmap_aff_from_point_pairs.h>
    test_hr2_pmap_aff_from_point_pairs(TRUE); // verbose

typedef struct hpmat_options_t 
  { char *prefix;         /* Prefix for output files. */
    hr2_pmap_t Aopt;      /* Expected optimum. */
    r3x3_t R;             /* Max elementwise adjustment in direct array of map. */
    int32_t ns;           /* Number of evaluations for plotting. */
  } hpmat_options_t;
  /* Command line parameters.
    
    The matrix {R} is the search region radius: it determines the
    entries of the projective map that are to be adjusted, and the
    amount of adjustment. Namely, entry {A.dir[i][j]} can be ajusted by
    {±R[i][j]}. In particular, element {A.dir[i][j]} will be fixed at
    {Aopt.dir[i][j]} if {R[i][j]} is zero. The elements of {R} must be
    non-negative. Because of the homoegneity property, at least one of
    the elements of {R} must be zero. */

void hpmat_one(hpmat_options_t *o, char *method, bool_t verbose);
  /* Tests affine map adjustment algorithm {method} ("quad" etc). */

void test_hr2_pmap_aff_from_point_pairs(bool_t verbose);

void test_hr2_pmap_aff_from_point_pairs(bool_t verbose)
  {
    /* Decide how many data point pairs to use: */
    
    int32_t np;
    if (drandom() < 0.250)
      { np = int32_abrandom(0,3); }
    else
      { np = int32_abrandom(4,20); }
     
    if (verbose) { fprintf(stderr, "--- hr2_pmap_aff_from_point_pairs  np = %d ---\n", np); }
     
    double eps = (np <= 3 ? 0.0 : 0.01); /* Typical perturbation size. */
    verbose = verbose | (eps == 0.0);

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
        M = hr2_pmap_translation(&v, +1, +1);
      }
    else if (np == 2)
      { /* M rotation, scaling, and translation: */
        double ang = dabrandom(0.0, 2.0)*M_PI;
        double scale = dabrandom(0.5, 2.0);
        hr2_pmap_t R = hr2_pmap_rotation_and_scaling(ang, scale);
        r2_t v; r2_throw_cube(&v);
        hr2_pmap_t T = hr2_pmap_translation(&v, +1, +1);
        M = hr2_pmap_compose(&R, &T);
      }
    else
      { throw_aff_map(&M); }
    
    /* Map the points through {M}, and add perturbation: */
    r2_t q[np];
    for (int32_t ip = 0; ip < np; ip++) 
      { q[ip] = hr2_pmap_r2_point(&(p[ip]), &M);
        if (eps > 0.0)
          { r2_t d; r2_throw_cube(&d);
            r2_mix(1.0, &(q[ip]), eps, &d, &(q[ip]));
          }
      }
    
    /* Compute the map: */
    hr2_pmap_t N = hr2_pmap_aff_from_point_pairs(np, p, q, w);
    
    if (debug)
      { print_pmap("M", &M);
        print_pmap("N", &N);
      }
    /* Compare the maps with {hr2_pmap_aff_discr_sqr}: */
    double mis2 = hr2_pmap_aff_discr_sqr(&M, &N);
    double mis = sqrt(mis2); 
    if (eps == 0.0)
      { check_num_eps("mis", mis, 0.0, 0.0000001, "hr2_pmap_aff_from_point_pairs failed"); }
    else 
      { /* Compare the maps by testing with given points: */
        assert (np > 0);
        double sum_d2A = 0.0;
        double sum_d2B = 0.0;
        for (int32_t ip = 0; ip < np; ip++)
          { r2_t pi = p[ip]; /* Given source point. */
            r2_t qi = q[ip]; /* Given destination point. */
            r2_t piA = hr2_pmap_r2_point(&pi, &M); /* Where the ideal map {M} sent it. */
            double d2A = r2_dist_sqr(&piA, &qi);
            sum_d2A += d2A;
            r2_t piB = hr2_pmap_r2_point(&pi, &N); /* Where the fitted map {N} sent it. */
            double d2B = r2_dist_sqr(&piB, &qi);
            sum_d2B += d2B;
          }
        double dA = sqrt(sum_d2A/np); /* RMS error of original data points rel to {M}. */
        double dB = sqrt(sum_d2B/np); /* RMS error of original data points rel to {N}. */
        if (verbose) { fprintf(stderr, "rms error M = %12.7f N = %12.7f\n", dA, dB); }
        double bad = fmax(0.0, dB - dA); /* Positive if {N} is worse than {M}, 0 if better. */
        check_num_eps("bad", bad, 0.0, np*eps*eps, "hr2_pmap_aff_from_point_pairs failed");
      }
  }

void hpmat_choose_initial_guess(hr2_pmap_t *Aopt, r3x3_t *R, hr2_pmap_t *Aini);
  /* Stores into {*Aini} the initail guess for the affine map
    to be adjusted, chosen at random in the box {Aopt±(R/2)}.  
    At least one of them must be zero. */

void hpmat_show_disp
  ( hr2_pmap_t *A, 
    char *Bname, 
    hr2_pmap_t *B,
    r3x3_t *R
  );
  /* Shows the discrepancy beweeen maps {A} and {B}.
    
    The arrays are compared with {hr2_pmap_diff_sqr(A,B)} and with
    {hpmat_diff_rel_sqr(D,R)} where {D=A.dir-B.dir}. If {print} is true
    also prints {D}. If the maps {A} and {B} are affine, they are also
    compared with {hr2_pmap_aff_discr_sqr(A,B)}. */

void hpmat_debug_map
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

void hpmat_plot_goal
  ( char *prefix,
    char *method,
    hr2_pmap_aff_from_point_pairs_func_t *F2, 
    hr2_pmap_t *A,
    r3x3_t *U,
    r3x3_t *V,
    int32_t ns
  );
  /* Writes a file "{prefix}_{method}.dat" with a random 2D slice of the goal
    function {F2(B)}, where {B} ranges in the neighborhood of {A}, in
    the plane definde by the displacements {±U} and {±V}.
    
    Specifically, generates the affine map as {B = A + u*U + v*V[k]},
    where {u} is {tu/ns} for {tu} in {-ns..+ns}, and {v} is {tv/ns} for
    {tv} in {-ns..+ns}. Plots {F2(B)} as a function of {u} and {v}. */

void hpmat_choose_plot_directions
  ( r3x3_t *R, 
    r3x3_t *U, 
    r3x3_t *V
  );
  /* Chooses two random matrices {U,V} that are orthogonal and of equal
    length when viewed as vectors of {\RR^9}. They have
    non-zero elements only where {R} itself is non-zero. Their norms are
    adjusted so that the rectangle {{a*U + b*V : a,b \in [-1_+1]}} lies
    within the box {±R}. */
  
double hpmat_mismatch_sqr_1(hr2_pmap_t *A, hr2_pmap_t *B);
  /* A squared mismatch function between {*A} and {*B},
    defined as {hr2_pmap_diff_sqr(A,B)}. */
  
double hpmat_mismatch_sqr_2(hr2_pmap_t *A, hr2_pmap_t *B);
  /* A squared mismatch function between {*A} and {*B},
    defined a positive definite quadratic funtion of 
    the six elements of {A^{-1} B - I} where {I} is the 
    identity.  This metric is used instead of the simple squared
    Euclidean distance (sum of squares of element differences)
    because the latter would make elementwise optimzation too easy. */
  
double hpmat_diff_rel_sqr(r3x3_t *D, r3x3_t *R);
  /* Computes the sum of squared elements of {D}, each divided by the
    corresponding element of {R}.  The elements that are zero in {R} 
    must be zero in {D} too, and are ignored. */

hpmat_options_t *hpmat_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */
  
void hpmat_parse_next_r3x3(argparser_t *pp, r3x3_t *R);
  /* Parses a 3x3 matrix from the command line, as 9 numbers 
    {R[0][0] R[0][1] R[0][2]  R[1][0] ... R[2][2]}. */

int32_t main(int32_t argn, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srand(4615*417);
    bool_t verbose = TRUE;
    
    hpmat_options_t *o = hpmat_parse_options(argc, argv);

    hpmat_one(o, "quad", verbose);

    return 0;
  }
  
void hpmat_one(hpmat_options_t *o, char *method, bool_t verbose)
  { 
    fprintf(stderr, "testing with method %s\n", method);
    
    bool_t debug_maps = FALSE; /* TRUE to print every probe point. */
    
    auto double F2_mismatch(hr2_pmap_t *A);
      /* A goal function that returns the squared mismatch between a map
        {*A} and the optimum {*Aopt}. Also prints the map if
        {debug_maps} is true. */
    
    hr2_pmap_t *Aopt = &(o->Aopt);
    r3x3_t *R = &(o->R);
    r3x3_gen_print(stderr, R, "%12.7f", "R = [ ","\n      "," ]\n", "[ "," "," ]");
    
    /* Choose the initial guess {Aini}: */
    hr2_pmap_t Aini;
    hpmat_choose_initial_guess(Aopt, R, &Aini);
       
    /* Choose two orthogonal perturbation vectors in {\RR^6} for plotting: */
    r3x3_t U, V;
    hpmat_choose_plot_directions(R, &U, &V); 
    
    /* Print raw function value and bias term for optimum map: */
    double F2opt = F2_mismatch(Aopt);
    hpmat_debug_map("actual optimum", TRUE, Aopt, NULL, NULL, R, F2opt);
    
    /* Print raw function value and bias term for initial map: */
    double F2ini = F2_mismatch(&Aini);
    hpmat_debug_map("initial guess ", TRUE, &Aini, NULL, Aopt, R, F2ini);
    
    /* Plot the goal function in the neighborhood of the initial guess: */ 
    fprintf(stderr, "Plotting goal function around initial point...\n");
    hpmat_plot_goal(o->prefix, method, F2_mismatch, &Aini, &U, &V, o->ns);
    
    /* Call the optimizer: */
    fprintf(stderr, "optimizing...\n");
    debug_maps = TRUE;
    hr2_pmap_t Asol = Aini;  /* Computed affine map. */
    double F2sol;   /* Goal function at {Asol}. */
    if (strcmp(method, "quad") == 0)
      { double tol = 0.02;
        hr2_pmap_aff_from_point_pairs_quad(F2_mismatch, R, tol, &Asol, &F2sol);
      }
    else
      { demand(FALSE, "invalid method"); }
    
    /* Print raw function value and bias term for computed optimum: */
    F2sol = F2_mismatch(&Asol);
    hpmat_debug_map("computed optimum", TRUE, &Asol, &Aini, Aopt, R, F2sol);

    return;
    
    double F2_mismatch(hr2_pmap_t *A)
      { double d2 = hpmat_mismatch_sqr_1(Aopt, A);
        if (debug_maps) 
          { hpmat_debug_map("probe map", verbose, A, &Aini, Aopt, R, d2); }
        return d2;
      }
  }

double hpmat_mismatch_sqr_1(hr2_pmap_t *A, hr2_pmap_t *B)
  { 
    double diff2 = hr2_pmap_diff_sqr(A, B);
    return diff2;
  }

double hpmat_mismatch_sqr_2(hr2_pmap_t *A, hr2_pmap_t *B)
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
  
void hpmat_choose_initial_guess(hr2_pmap_t *Aopt, r3x3_t *R, hr2_pmap_t *Aini)
  {
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { double frac = dabrandom(-0.50, +0.50);
            Aini->dir.c[i][j] = Aopt->dir.c[i][j] + frac*R->c[i][j];
          }
      }
      r3x3_inv(&(Aini->dir), &(Aini->inv));
  }
    
void hpmat_plot_goal
  ( char *prefix,
    char *method,
    hr2_pmap_aff_from_point_pairs_func_t *F2, 
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
    asprintf(&fname, "%s_%s.dat", prefix, method);
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

void hpmat_choose_plot_directions(r3x3_t *R, r3x3_t *U, r3x3_t *V)
  { 
    /* Identify the elements of {R} that are not zero, get the addresses {Up,Vp}: */
    int32_t ne = 0;
    double Re[9]; /* Nonzero elements of {R} are  {Re[0..ne-1]}. */
    double *Up[9], *Vp[9]; /* Pointers of elements of {U} and {V} where {R} is nonzero. */
    
    hr2_pmap_aff_from_point_pairs_get_var_elems(U, R, Up, Re, &ne);
    hr2_pmap_aff_from_point_pairs_get_var_elems(V, R, Vp, Re, &ne);
      
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

double hpmat_diff_rel_sqr(r3x3_t *D, r3x3_t *R)
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

void hpmat_show_disp
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
    double drel2 = hpmat_diff_rel_sqr(&D, R);
    fprintf(stderr, "  diff = %12.7f drel = %12.7f", sqrt(diff2), sqrt(drel2));
    
    if (hr2_pmap_is_affine(A) && hr2_pmap_is_affine(B))
      { /* Difference between {A} and {B} over the unit circle: */
        double mism2 = hr2_pmap_aff_discr_sqr(A, B);
        fprintf(stderr, " aff mism = %12.7f", sqrt(mism2));
      }
    fprintf(stderr, "\n");
  }

void hpmat_debug_map
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

    if (Aini != NULL) { hpmat_show_disp(A, "Aini", Aini, R); }
    if (Aopt != NULL) { hpmat_show_disp(A, "Aopt", Aopt, R); }

    fprintf(stderr, "  F2(A) = %22.10f\n", F2A);
    fprintf(stderr, "\n");
    return;
   }

hpmat_options_t *hpmat_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, "HELP!");
    argparser_set_info(pp, "KNOW THYSELF");
    argparser_process_help_info_options(pp);
    
    hpmat_options_t *o = notnull(malloc(sizeof(hpmat_options_t)), "no mem");

    argparser_get_keyword(pp, "-outPrefix");
    o->prefix = argparser_get_next(pp);

    argparser_get_keyword(pp, "-optimum");
    hpmat_parse_next_r3x3(pp, &(o->Aopt.dir));  /* Target map (expected optimum). */
    r3x3_inv(&(o->Aopt.dir), &(o->Aopt.inv));

    argparser_get_keyword(pp, "-deviation");
    hpmat_parse_next_r3x3(pp, &(o->R));     /* Max elementwise deviation. */
    
    argparser_get_keyword(pp, "-nPlot");
    o->ns = (int32_t)argparser_get_next_uint(pp, 1, 500); /* Number of evaluations for plotting. */

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }

void hpmat_parse_next_r3x3(argparser_t *pp, r3x3_t *R)
  {
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { R->c[i][j] = argparser_get_next_double(pp, -100.0, +100.0); }
      }
  }

