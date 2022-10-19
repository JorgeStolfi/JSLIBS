#define PROG_NAME "hr2_pmap_adjust_test"
#define PROG_DESC "test of {hr2_pmap_adjust.h}"
#define PROG_VERS "1.0"

/* Last edited on 2022-03-01 12:05:43 by stolfi */ 
/* Created on 2020-07-11 by J. Stolfi, UNICAMP */
/* Based on {test_align.c} by J. Stolfi, UNICAMP */

#define hr2_pmap_adjust_test_COPYRIGHT \
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

#include <hr2_pmap_adjust.h>

typedef struct hpmat_options_t 
  { char *prefix;         /* Prefix for output files. */
    hr2_pmap_t Aopt;      /* Expected optimum. */
    r3x3_t R;             /* Max elementwise adjustment in direct array of map. */
    int32_t ns;           /* Number of evaluations for plotting. */
  } hpmat_options_t;
  /* Command line parameters. */


void hpmat_test_hr2_pmap_from_many_pairs(bool_t verbose);

void hpmat_one(hpmat_options_t *o, char *method);
  /* Tests affine map adjustment algorithm {method} ("quad" etc). */

void hpmat_choose_initial_guess(hr2_pmap_t *Aopt, r3x3_t *R, hr2_pmap_t *Aini);
  /* Stores into {*Aini} the initail guess for the affine map
    to be adjusted, chosen at random in the box {Aopt±(R/2)}. 
    The elements of the region radius {R} must be non-negative. */

void hpmat_show_disp
  ( hr2_pmap_t *A, 
    double f2sol,
    hr2_pmap_t *Aini, 
    hr2_pmap_t *Aopt, 
    r3x3_t *R
  ); 
  /* Shows the discrepancy beweeen the map {A} against known optimum {Aopt} (if not NULL).
    The solution is correct if every affine map {A[i]} differs from {Aopt[i]}
    Also checks that the displacements from {Aini[i]}
    to {A[i]} are within the search {±R[i]}. */

void hpmat_debug_map
  ( char *label, 
    hr2_pmap_t *A, 
    hr2_pmap_t *Aini, 
    hr2_pmap_t *Aopt, 
    r3x3_t *R,
    double f2p
  );
  /* Prints the affine map {*A} and the corresponding raw goal function value {f2p}.
    If {Aini} is not {NULL}, prints also the difference between 
    {*A} and the corresponding {*Aini}.  Ditto for the difference to {Aopt},
    if {Aopt} is not {NULL}. */

void hpmat_plot_goal
  ( char *prefix,
    char *method,
    hr2_pmap_adjust_func_t *f2, 
    hr2_pmap_t *A,
    r3x3_t *U,
    r3x3_t *V,
    int32_t ns
  );
  /* Writes a file "{prefix}_{method}.dat" with a random 2D slice of the goal
    function {f2(B)}, where {B} ranges in the neighborhood of {A}, in
    the plane definde by the displacements {±U} and {±V}.
    
    Specifically, generates the affine map as {B = A + u*U + v*V[k]},
    where {u} is {tu/ns} for {tu} in {-ns..+ns}, and {v} is {tv/ns} for
    {tv} in {-ns..+ns}. Plots {f2(B)} as a function of {u} and {v}. */

void hpmat_choose_plot_directions
  ( r3x3_t *R, 
    r3x3_t *U, 
    r3x3_t *V
  );
  /* Chooses two random affine maps {U,V} that are orthogonal and of equal length when viewed as 
    vectors of {\RR^6}. Their norms are adjusted so that they are within the box {±R}. */
  
double hpmat_mismatch_sqr(hr2_pmap_t *A, hr2_pmap_t *B);
  /* A squared mismatch function between {*A} and {*B},
    defined a positive definite quadratic funtion of 
    the six elements of {A^{-1} B - I} where {I} is the 
    identity.  This metric is used instead of the simple squared
    Euclidean distance (sum of squares of element differences)
    because the latter would make elementwise optimzation too easy. */

hpmat_options_t *hpmat_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */
  
void hpmat_parse_next_affine_map(argparser_t *pp, hr2_pmap_t *A);
  /* Parses an affine map from the command line, as 6 numbers 
    {m[0][0] m[0][1] m[1][0] m[1][1] d[0] d[1]},
    where {m = A->dir} and {d = A->disp}. */

int32_t main(int32_t argn, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srand(4615*417);
    bool_t verbose = TRUE;
    
    hpmat_options_t *o = hpmat_parse_options(argc, argv);

    hpmat_one(o, "quad");
    /* TEST: hr2_pmap_t hr2_pmap_from_many_pairs */
    hpmat_test_hr2_pmap_from_many_pairs(verbose);

    return 0;
  }
  
void hpmat_one(hpmat_options_t *o, char *method)
  { 
    fprintf(stderr, "testing with method %s\n", method);
    
    bool_t debug_maps = FALSE; /* TRUE to print every probe point. */
    
    auto double f2_mismatch(hr2_pmap_t *A);
      /* A goal function that returns the squared mismatch between a map
        {*A} and the optimum {*Aopt}. Also prints the map if
        {debug_maps} is true. It does not use the square Frobenius norm,
        because that would make elementwise optimization too easy. */
    
    hr2_pmap_t *Aopt = &(o->Aopt);
    r3x3_t *R = &(o->R);
    
    /* Choose the initial guess {Aini}: */
    hr2_pmap_t Aini;
    hpmat_choose_initial_guess(Aopt, R, &Aini);
       
    /* Choose two orthogonal perturbation vectors in {\RR^6} for plotting: */
    r3x3_t U, V;
    hpmat_choose_plot_directions(R, &U, &V); 
    
    /* Print raw function value and bias term for optimum map: */
    double f2opt = f2_mismatch(Aopt);
    hpmat_debug_map("actual optimum", Aopt, NULL, NULL, R, f2opt);
    
    /* Print raw function value and bias term for initial map: */
    double f2ini = f2_mismatch(&Aini);
    hpmat_debug_map("initial guess ", &Aini, NULL, Aopt, R, f2ini);
    
    /* Plot the goal function in the neighborhood of the initial guess: */ 
    fprintf(stderr, "Plotting goal function around initial point...\n");
    hpmat_plot_goal(o->prefix, method, f2_mismatch, &Aini, &U, &V, o->ns);
    
    /* Call the optimizer: */
    fprintf(stderr, "optimizing");
    debug_maps = TRUE;
    hr2_pmap_t Asol = Aini;  /* Computed affine map. */
    double f2sol;   /* Goal function at {Asol}. */
    if (strcmp(method, "quad") == 0)
      { double tol = 0.02;
        hr2_pmap_adjust_quad(f2_mismatch, R, tol, &Asol, &f2sol);
      }
    else
      { demand(FALSE, "invalid method"); }
    
    /* Print raw function value and bias term for computed optimum: */
    f2sol = f2_mismatch(&Asol);
    hpmat_debug_map("computed optimum", &Asol, &Aini, Aopt, R, f2sol);
    
    hpmat_show_disp(&Asol, f2sol, &Aini, Aopt, R);

    return;
    
    double f2_mismatch(hr2_pmap_t *A)
      { double d2 = hpmat_mismatch_sqr(Aopt, A);
        if (debug_maps) 
          { hpmat_debug_map("probe map   ", A, &Aini, Aopt, R, d2); }
        return d2;
      }
  }
  
void hpmat_test_hr2_pmap_from_many_pairs(bool_t verbose)
  { fprintf(stderr, "!! {hr2_pmap_from_many_pairs} NOT TESTED\n"); }

double hpmat_mismatch_sqr(hr2_pmap_t *A, hr2_pmap_t *B)
  { 
    hr2_pmap_t ABdif = hr2_pmap_inv_comp(A, B); /* The map {A^{-1} B}. */
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
    hr2_pmap_adjust_func_t *f2, 
    hr2_pmap_t *A,
    r3x3_t *U,
    r3x3_t *V,
    int32_t ns
  )
  {
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
            double f2B = f2(&B);
            fprintf(wr, "%+9.6f %+9.6f  %12.6f\n", du, dv, f2B); 
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
    double Re[8]; /* Nonzero elements of {R} are  {Re[0..ne-1]}. */
    double *Up[8], *Vp[8]; /* Pointers of elements of {U} and {V} where {R} is nonzero. */
    
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { double *Rp = &(R->c[i][j]);
            assert((*Rp) >= 0.0);
            if ((*Rp) != 0)
              { demand(ne < 8, "at least one element of the array must be fixed");
                Re[ne] = (*Rp);
                Up[ne] = &(U->c[i][j]); /* Pointer to element of {*U}. */
                Vp[ne] = &(V->c[i][j]); /* Pointer to element of {*V}. */
                ne++;
              }
          }
      }
      
    /* Throw a random vector {ue[0..ne-1]}:*/
    double ue[ne];
    rn_throw_dir(ne, ue);
    
    /* Throw a vector {ve[0..ne-1]} orthogonal to {ue}: */
    double ve[ne];
    do
      { rn_throw_dir(ne, ve);
        rn_decomp(ne, ve, ue, NULL, ve);
      }
    while (rn_norm_sqr(ne, ve) < 1.0e-4);
    
    /* Adjust the lengths to be equal to the length of {R}: */
    double Rnorm = rn_norm(ne, Re);
    double unorm = rn_norm(ne, ue);
    double vnorm = rn_norm(ne, ve);
    
    /* Insert {ue,ve} into {U,V}: */
    for (int32_t ie = 0; ie < ne; ie++)
      { (*(Up[ie])) = Rnorm*ue[ie]/unorm;
        (*(Vp[ie])) = Rnorm*ve[ie]/vnorm;
      }
  }

void hpmat_show_disp
  ( 
    hr2_pmap_t *A, 
    double f2sol,
    hr2_pmap_t *Aini,
    hr2_pmap_t *Aopt, 
    r3x3_t *R
  )
  {
    auto void compare_to_map(char *tag, hr2_pmap_t *Aref);
      /* Prints the RMS absolute difference {A - Aref}, and 
        the RMS difference {A - Aref} relative to {R}.  
        Ignores any element of {A,Aref} if the corresponding
        element of {R} is zero. */
    
    fprintf(stderr, "checking the solution\n");
    
    if (Aopt != NULL) { compare_to_map("optimum", Aopt); }
    if (Aini != NULL) { compare_to_map("initial", Aini); }
    
    fprintf(stderr, "\n");
    
    void compare_to_map(char *tag, hr2_pmap_t *B)
      { 
        double dabs2, drel2;
        hr2_pmap_dist_sqr(A, B, R, &dabs2, &drel2);
        double dabs = sqrt(dabs2);
        double drel = sqrt(drel2);
        fprintf(stderr, "  disp from %s:", tag);
        fprintf(stderr, " abs = %10.6f rel = %10.6f\n", dabs, drel);
      }
  }

void hpmat_debug_map
  ( char *label, 
    hr2_pmap_t *A, 
    hr2_pmap_t *Aini,
    hr2_pmap_t *Aopt, 
    r3x3_t *R,
    double f2A
  )
  { 
    auto void show_disp(char *tag, hr2_pmap_t *C, r3x3_t *R);
      /* Shows the difference {A} to {C}, abslolute and relative. */

    fprintf(stderr, "%s\n", label);
    fprintf(stderr, "  A =        ");
    hr2_pmap_gen_print(stderr, A, "%12.7f", "%12.7f", "[ "," "," ]","[ "," "," ]"," + "); 
    fprintf(stderr, "\n");
    if (Aini != NULL) { show_disp("Aini", Aini, R); }
    if (Aopt != NULL) { show_disp("Aopt", Aopt, R); }
    fprintf(stderr, "  f2(A) = %22.10f\n", f2A);
    fprintf(stderr, "\n");
    return;

    void show_disp(char *tag, hr2_pmap_t *Aref, r3x3_t *r) 
      {
        fprintf(stderr, "    = %-4s + ", tag);
        hr2_pmap_t Adif;
        r3x3_sub(&(A->dir), &(Aref->dir), &(Adif.dir));
        r2_sub(&(A->disp), &(Aref->disp), &(Adif.disp));
        hr2_pmap_gen_print(stderr, &Adif, "%12.7f", "%12.7f", "[ "," "," ]","[ "," "," ]"," + "); 
        double dabs2, drel2;
        hr2_pmap_disp_sqr(A, Aref, R, &dabs2, &drel2);
        fprintf(stderr, " dabs = %12.7f drel = %12.7f", sqrt(dabs2), sqrt(drel2));
        double mis2 = hr2_pmap_aff_mismatch_sqr(A, Aref);
        fprintf(stderr, " mis = %12.7f\n", sqrt(mis2));
      }
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
    hpmat_parse_next_affine_map(pp, &(o->Aopt));  /* Target map (expected optimum). */

    argparser_get_keyword(pp, "-deviation");
    hpmat_parse_next_affine_map(pp, &(o->R));     /* Max emementwise deviation. */
    
    argparser_get_keyword(pp, "-nPlot");
    o->ns = (int32_t)argparser_get_next_uint(pp, 1, 500); /* Number of evaluations for plotting. */

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }

void hpmat_parse_next_affine_map(argparser_t *pp, hr2_pmap_t *A)
  {
    for (int32_t i = 0; i < 2; i++)
      { for (int32_t j = 0; j < 2; j++)
          { A->dir.c[i][j] = argparser_get_next_double(pp, -100.0, +100.0); }
      }
    for (int32_t j = 0; j < 2; j++)
      { A->disp.c[j] = argparser_get_next_double(pp, -100.0, +100.0); }
  }

