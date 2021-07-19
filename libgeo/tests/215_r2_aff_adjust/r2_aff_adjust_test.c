#define PROG_NAME "r2_aff_adjust_test"
#define PROG_DESC "test of {r2_aff_adjust.h}"
#define PROG_VERS "1.0"

/* Last edited on 2020-10-16 03:39:48 by jstolfi */ 
/* Created on 2020-07-11 by J. Stolfi, UNICAMP */
/* Based on {test_align.c} by J. Stolfi, UNICAMP */

#define r2_aff_adjust_test_COPYRIGHT \
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
#include <rn.h>
#include <r2x2.h>
#include <bool.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>
#include <assert.h>

#include <r2_aff_adjust.h>

typedef struct r2aat_options_t 
  { char *prefix;       /* Prefix for output files. */
    r2_aff_map_t Aopt;  /* Target map (expected optimum). */
    r2_aff_map_t R;     /* Max emementwise deviation. */
    int32_t ns;         /* Number of evaluations for plotting. */
  } r2aat_options_t;
  /* Command line parameters. */

void r2aat_one(r2aat_options_t *o, char *method);
  /* Tests affine map adjustment algorithm {method} ("quad" etc). */

void r2aat_choose_initial_guess(r2_aff_map_t *Aopt, r2_aff_map_t *R, r2_aff_map_t *Aini);
  /* Stores into {*Aini} the initail guess for the affine map
    to be adjusted, chosen at random in the box {Aopt±(R/2)}. 
    The elements of the region radius {R} must be non-negative. */

void r2aat_show_disp
  ( r2_aff_map_t *A, 
    double f2sol,
    r2_aff_map_t *Aini, 
    r2_aff_map_t *Aopt, 
    r2_aff_map_t *R
  ); 
  /* Shows the discrepancy beweeen the map {A} against known optimum {Aopt} (if not NULL).
    The solution is correct if every affine map {A[i]} differs from {Aopt[i]}
    Also checks that the displacements from {Aini[i]}
    to {A[i]} are within the search {±R[i]}. */

void r2aat_debug_map
  ( char *label, 
    r2_aff_map_t *A, 
    r2_aff_map_t *Aini, 
    r2_aff_map_t *Aopt, 
    r2_aff_map_t *R,
    double f2p
  );
  /* Prints the affine map {*A} and the corresponding raw goal function value {f2p}.
    If {Aini} is not {NULL}, prints also the difference between 
    {*A} and the corresponding {*Aini}.  Ditto for the difference to {Aopt},
    if {Aopt} is not {NULL}. */

void r2aat_plot_goal
  ( char *prefix,
    char *method,
    r2_aff_adjust_func_t *f2, 
    r2_aff_map_t *A,
    r2_aff_map_t *U,
    r2_aff_map_t *V,
    int32_t ns
  );
  /* Writes a file "{prefix}_{method}.dat" with a random 2D slice of the goal
    function {f2(B)}, where {B} ranges in the neighborhood of {A}, in
    the plane definde by the displacements {±U} and {±V}.
    
    Specifically, generates the affine map as {B = A + u*U + v*V[k]},
    where {u} is {tu/ns} for {tu} in {-ns..+ns}, and {v} is {tv/ns} for
    {tv} in {-ns..+ns}. Plots {f2(B)} as a function of {u} and {v}. */

void r2aat_choose_plot_directions
  ( r2_aff_map_t *R, 
    r2_aff_map_t *U, 
    r2_aff_map_t *V
  );
  /* Chooses two random affine maps {U,V} that are orthogonal and of equal length when viewed as 
    vectors of {\RR^6}. Their norms are adjusted so that they are within the box {±R}. */
  
double r2aat_mismatch_sqr(r2_aff_map_t *A, r2_aff_map_t *B);
  /* A squared mismatch function between {*A} and {*B},
    defined a positive definite quadratic funtion of 
    the six elements of {A^{-1} B - I} where {I} is the 
    identity.  This metric is used instead of the simple squared
    Euclidean distance (sum of squares of element differences)
    because the latter would make elementwise optimzation too easy. */

r2aat_options_t *r2aat_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */
  
void r2aat_parse_next_affine_map(argparser_t *pp, r2_aff_map_t *A);
  /* Parses an affine map from the command line, as 6 numbers 
    {m[0][0] m[0][1] m[1][0] m[1][1] d[0] d[1]},
    where {m = A->mat} and {d = A->disp}. */

int32_t main(int32_t argn, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srand(4615*417);
    
    r2aat_options_t *o = r2aat_parse_options(argc, argv);

    r2aat_one(o, "quad");
    return 0;
  }
  
void r2aat_one(r2aat_options_t *o, char *method)
  { 
    fprintf(stderr, "testing with method %s\n", method);
    
    bool_t debug_maps = FALSE; /* TRUE to Arint every Arobe Aoint. */
    
    auto double f2_mismatch(r2_aff_map_t *A);
      /* A goal function that returns the squared mismatch {*A} to the
        optimum {*Aopt}. Also prints the map if {debug_maps} is true.
        It does not use the square Frobenisu norm,
        becaus that would make elementwise optimzation too easy. */
    
    r2_aff_map_t *Aopt = &(o->Aopt);
    r2_aff_map_t *R = &(o->R);
    
    /* Choose the initial guess {Aini}: */
    r2_aff_map_t Aini;
    r2aat_choose_initial_guess(Aopt, R, &Aini);
       
    /* Choose two orthogonal perturbation vectors in {\RR^6} for plotting: */
    r2_aff_map_t U, V;
    r2aat_choose_plot_directions(R, &U, &V); 
    
    /* Print raw function value and bias term for optimum map: */
    double f2opt = f2_mismatch(Aopt);
    r2aat_debug_map("actual optimum", Aopt, NULL, NULL, R, f2opt);
    
    /* Print raw function value and bias term for initial map: */
    double f2ini = f2_mismatch(&Aini);
    r2aat_debug_map("initial guess ", &Aini, NULL, Aopt, R, f2ini);
    
    /* Plot the goal function in the neighborhood of the initial guess: */ 
    fprintf(stderr, "Plotting goal function around initial point...\n");
    r2aat_plot_goal(o->prefix, method, f2_mismatch, &Aini, &U, &V, o->ns);
    
    /* Call the optimizer: */
    fprintf(stderr, "optimizing");
    debug_maps = TRUE;
    r2_aff_map_t Asol = Aini;  /* Computed affine map. */
    double f2sol;   /* Goal function at {Asol}. */
    if (strcmp(method, "quad") == 0)
      { double tol = 0.02;
        r2_aff_adjust_quad(f2_mismatch, R, tol, &Asol, &f2sol);
      }
    else
      { demand(FALSE, "invalid method"); }
    
    /* Print raw function value and bias term for computed optimum: */
    f2sol = f2_mismatch(&Asol);
    r2aat_debug_map("computed optimum", &Asol, &Aini, Aopt, R, f2sol);
    
    r2aat_show_disp(&Asol, f2sol, &Aini, Aopt, R);

    return;
    
    double f2_mismatch(r2_aff_map_t *A)
      { double d2 = r2aat_mismatch_sqr(Aopt, A);
        if (debug_maps) 
          { r2aat_debug_map("probe map   ", A, &Aini, Aopt, R, d2); }
        return d2;
      }
  }
  
double r2aat_mismatch_sqr(r2_aff_map_t *A, r2_aff_map_t *B)
  { 
    /* Compute the map {A^{-1} B}: */
    r2_aff_map_t Ainv, Adif;
    r2_aff_map_invert(A, &Ainv);
    r2_aff_map_compose(&Ainv, B, &Adif);
    /* Compare the matrix of {Adif} with the identity map: */
    double d2 = 0;
    for (int32_t i = 0; i < 2; i++)
      { for (int32_t j = 0; j < 2; j++)
          { double d = Adif.mat.c[i][j] - (i == j ? 1.0 : 0.0);
            d2 += d*d;
          }
      }
    /* Compare the displacement with a skew distance squared metric: */
    double u = Adif.disp.c[0] + Adif.disp.c[1];
    double v = Adif.disp.c[0] - Adif.disp.c[1]; 
    d2 += (u*u + 5*v*v);
    return d2;
  }
  
void r2aat_choose_initial_guess(r2_aff_map_t *Aopt, r2_aff_map_t *R, r2_aff_map_t *Aini)
  {
    for (int32_t j = 0; j < 2; j++)
      { for (int32_t i = 0; i < 2; i++)
          { double frac = dabrandom(-0.50, +0.50);
            Aini->mat.c[i][j] = Aopt->mat.c[i][j] + frac*R->mat.c[i][j];
          }
        double frac = dabrandom(-0.50, +0.50);
        Aini->disp.c[j] = Aopt->disp.c[j] + frac*R->disp.c[j];
      }
  }
    
void r2aat_plot_goal
  ( char *prefix,
    char *method,
    r2_aff_adjust_func_t *f2, 
    r2_aff_map_t *A,
    r2_aff_map_t *U,
    r2_aff_map_t *V,
    int32_t ns
  )
  {
    /* Sweep the {A,U,V} plane and Plot: */
    char *fname = NULL;
    asprintf(&fname, "%s_%s.dat", prefix, method);
    FILE *wr = open_write(fname, TRUE);
    
    r2_aff_map_t B; /* Probe map. */
    fprintf(stderr, "\n");
    for (int32_t iu = -ns; iu <= +ns; iu++)
      { double du = ((double)iu)/((double)ns);
        fprintf(stderr, ".");
        for (int32_t iv = -ns; iv <= +ns; iv++)
          { double dv = ((double)iv)/((double)ns);
            /* Compute the probe map {B}. */
            r2x2_mix(du, &(U->mat), dv, &(V->mat), &(B.mat)); 
            r2x2_add(&(A->mat), &(B.mat), &(B.mat));
            r2_mix(du, &(U->disp), dv, &(V->disp), &(B.disp)); 
            r2_add(&(A->disp), &(B.disp), &(B.disp));
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

void r2aat_choose_plot_directions(r2_aff_map_t *R, r2_aff_map_t *U, r2_aff_map_t *V)
  { 
    /* Identify the elements of {R}  that are not zero, get the addresses of those elems in {U,V}: */
    int32_t ne = 0;
    double Re[6]; /* Nonzero elements of {R} are  {Re[0..ne-1]}. */
    double *Up[6], *Vp[6]; /* Pointers of elements of {U} and {V} where {R} is nonzero. */
    for (int32_t s = 0; s < 6; s++) 
      { double *Rp = r2_aff_map_elem_addr(R, s); /* Pointer to element of {*R}. */
        assert((*Rp) >= 0.0);
        if ((*Rp) != 0)
          { Re[ne] = (*Rp);
            Up[ne] = r2_aff_map_elem_addr(U, s); /* Pointer to element of {*A}. */
            Vp[ne] = r2_aff_map_elem_addr(V, s); /* Pointer to element of {*A}. */
            ne++;
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
    r2x2_zero(&(U->mat)); r2_zero(&(U->disp));
    r2x2_zero(&(V->mat)); r2_zero(&(V->disp));
    for (int32_t ie = 0; ie < ne; ie++)
      { (*(Up[ie])) = Rnorm*ue[ie]/unorm;
        (*(Vp[ie])) = Rnorm*ve[ie]/vnorm;
      }
  }

void r2aat_show_disp
  ( 
    r2_aff_map_t *A, 
    double f2sol,
    r2_aff_map_t *Aini,
    r2_aff_map_t *Aopt, 
    r2_aff_map_t *R
  )
  {
    auto void compare_to_map(char *tag, r2_aff_map_t *Aref);
      /* Prints the RMS absolute difference {A - Aref}, and 
        the RMS difference {A - Aref} relative to {R}.  
        Ignores any element of {A,Aref} if the corresponding
        element of {R} is zero. */
    
    fprintf(stderr, "checking the solution\n");
    
    if (Aopt != NULL) { compare_to_map("optimum", Aopt); }
    if (Aini != NULL) { compare_to_map("initial", Aini); }
    
    fprintf(stderr, "\n");
    
    void compare_to_map(char *tag, r2_aff_map_t *B)
      { 
        double dabs2, drel2;
        r2_aff_map_disp_sqr(A, B, R, &dabs2, &drel2);
        double dabs = sqrt(dabs2);
        double drel = sqrt(drel2);
        fprintf(stderr, "  disp from %s:", tag);
        fprintf(stderr, " abs = %10.6f rel = %10.6f\n", dabs, drel);
      }
  }

void r2aat_debug_map
  ( char *label, 
    r2_aff_map_t *A, 
    r2_aff_map_t *Aini,
    r2_aff_map_t *Aopt, 
    r2_aff_map_t *R,
    double f2A
  )
  { 
    auto void show_disp(char *tag, r2_aff_map_t *C, r2_aff_map_t *R);
      /* Shows the difference {A} to {C}, abslolute and relative. */

    fprintf(stderr, "%s\n", label);
    fprintf(stderr, "  A =        ");
    r2_aff_map_gen_print(stderr, A, "%12.7f", "%12.7f", "[ "," "," ]","[ "," "," ]"," + "); 
    fprintf(stderr, "\n");
    if (Aini != NULL) { show_disp("Aini", Aini, R); }
    if (Aopt != NULL) { show_disp("Aopt", Aopt, R); }
    fprintf(stderr, "  f2(A) = %22.10f\n", f2A);
    fprintf(stderr, "\n");
    return;

    void show_disp(char *tag, r2_aff_map_t *Aref, r2_aff_map_t *r) 
      {
        fprintf(stderr, "    = %-4s + ", tag);
        r2_aff_map_t Adif;
        r2x2_sub(&(A->mat), &(Aref->mat), &(Adif.mat));
        r2_sub(&(A->disp), &(Aref->disp), &(Adif.disp));
        r2_aff_map_gen_print(stderr, &Adif, "%12.7f", "%12.7f", "[ "," "," ]","[ "," "," ]"," + "); 
        double dabs2, drel2;
        r2_aff_map_disp_sqr(A, Aref, R, &dabs2, &drel2);
        fprintf(stderr, " dabs = %12.7f drel = %12.7f", sqrt(dabs2), sqrt(drel2));
        double mis2 = r2_aff_map_mismatch_sqr(A, Aref);
        fprintf(stderr, " mis = %12.7f\n", sqrt(mis2));
      }
   }

r2aat_options_t *r2aat_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, "HELP!");
    argparser_set_info(pp, "KNOW THYSELF");
    argparser_process_help_info_options(pp);
    
    r2aat_options_t *o = notnull(malloc(sizeof(r2aat_options_t)), "no mem");

    argparser_get_keyword(pp, "-outPrefix");
    o->prefix = argparser_get_next(pp);

    argparser_get_keyword(pp, "-optimum");
    r2aat_parse_next_affine_map(pp, &(o->Aopt));  /* Target map (expected optimum). */

    argparser_get_keyword(pp, "-deviation");
    r2aat_parse_next_affine_map(pp, &(o->R));     /* Max emementwise deviation. */
    
    argparser_get_keyword(pp, "-nPlot");
    o->ns = (int32_t)argparser_get_next_uint(pp, 1, 500); /* Number of evaluations for plotting. */

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }

void r2aat_parse_next_affine_map(argparser_t *pp, r2_aff_map_t *A)
  {
    for (int32_t i = 0; i < 2; i++)
      { for (int32_t j = 0; j < 2; j++)
          { A->mat.c[i][j] = argparser_get_next_double(pp, -100.0, +100.0); }
      }
    for (int32_t j = 0; j < 2; j++)
      { A->disp.c[j] = argparser_get_next_double(pp, -100.0, +100.0); }
  }

