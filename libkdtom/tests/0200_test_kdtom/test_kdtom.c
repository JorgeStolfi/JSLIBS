#define PROG_NAME "test_kdtom"
#define PROG_DESC "Test the k-d-tree rep of multi-dim sample arrays"
#define PROG_VERS "1.0"

#define tkdt_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"

/* Last edited on 2021-07-03 13:22:57 by jstolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    " argparser_help_info_HELP

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program tests the basic functions of the {libkdtom} library.\n" \
  "\n" \
  "OPTIONS\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  The Taj Mahal.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2021-06-25 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  " 2025-06-25 J. Stolfi: Created.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " tkdt_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>
#include <ppv_array.h>

#include <kdtom.h>
#include <kdtom_array.h>
#include <kdtom_split.h>
#include <kdtom_const.h>

/* COMMAND-LINE OPTIONS */

typedef struct tkdt_options_t
  { int32_t dummy;
  } tkdt_options_t;

/* INTERNAL PROTOTYPES */

    /* The unit for all linear dimensions is the voxel side. */ 
    
tkdt_options_t *tkdt_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {tkdt_options_t}. */

void tkdt_paint_bullseye(ppv_array_t *A, double ctr[], double R);
  /* Paints a suitable test pattern into {A}. Namely, two thin spherical
    shells where the original values of {A} are preserved,
    constant values outside, between, and inside the shells,
    and a different value at the center. The pattern will 
    be centered at the point {ctr}, and the outer radius of the
    outer shell will be {R}. */
    
void tkdt_choose_array_size(ppv_dim_t d, ppv_size_t sz[]);
  /* Fills {sz[0..d-1]} with suitable values, not all equal. */

ppv_array_t *tkdt_make_array(ppv_dim_t d, ppv_size_t sz[], ppv_nbits_t bps);
  /* Creates an array {A} with the given attributes {d,sz[0..d-1],bps},
    and fills it with a suitable test pattern. */

void tkdt_check_tree(kdtom_t *T);
  /* Prints data about the tree {T}. */

void tkdt_plot(kdtom_t *T);
  /* Plots the array {T}.  Only if {T.d} is 2. */

void tkdt_do_tests(tkdt_options_t *o,ppv_dim_t d, ppv_nbits_t bps);
void tkdt_test_array_node(ppv_array_t *A);
void tkdt_test_grinding(ppv_array_t *A);

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    tkdt_options_t *o = tkdt_parse_options(argc, argv);
    
    for (ppv_dim_t d = 0; d <= 16; d++)
      for (ppv_nbits_t bps = 0; bps <= 1;  bps = (ppv_nbits_t)((5*bps+4)/3))
        { tkdt_do_tests(o, d, bps); }
      
    return 0;
  }
  
void tkdt_do_tests(tkdt_options_t *o, ppv_dim_t d, ppv_nbits_t bps)
  { 
    fprintf(stderr, "Testing d = %u bps = %u ...\n", d, bps);
    ppv_size_t sz[d];
    tkdt_choose_array_size(d, sz);
    ppv_array_t *A = tkdt_make_array(d, sz, bps);
    tkdt_test_array_node(A);
    tkdt_test_grinding(A);
    return;
 }

ppv_array_t *tkdt_make_array(ppv_dim_t d, ppv_size_t sz[], ppv_nbits_t bps)
  {
    ppv_nbits_t bpw = ppv_best_bpw(bps);
    /* ppv_sample_t maxval = (ppv_sample_t)((1 << bps) - 1); */
    ppv_array_t *A = ppv_array_new(d, sz, bps, bpw);
    
    if (bps == 0) { /* Array must be all zeros anyway: */ return A; }
    
    /* Fill the array with random samples: */
    srandom(4615);
    ppv_throw_noise(A); 

    if (d == 0) { /* Only one sample: */  return A; }

    /* Compute min size {szmin} and cords {ctr} of domain center: */
    ppv_size_t szmin = ppv_MAX_SIZE;
    ppv_size_t szmax = 0;
    double ctr[d]; /* Center of bullseye. */
    for (ppv_axis_t k = 0; k < d; k++)
      { ppv_size_t szk = A->size[k];
        if (szk < szmin) { szmin = szk; }
        if (szk > szmax) { szmax = szk; }
        ctr[k] = 0.5*(double)szk;
      }
    if (szmax == 0) { /* Empty aray: */ return A; }
    if (szmin >= 3) 
      { /* Paint the test pattern: */
        double R = 0.50*(double)szmin; /* Radius of max inscribed ball. */
        tkdt_paint_bullseye(A, ctr, R);
      }
    return A;
  }
    
void tkdt_test_grinding(ppv_array_t *A)
  {
    kdtom_t *Tg = kdtom_grind_array(A);
    tkdt_check_tree(Tg);
    return;
  }

void tkdt_test_array_node(ppv_array_t *A)
  {
    kdtom_array_t *Ta = kdtom_array_make(A);
    tkdt_check_tree((kdtom_t *)Ta);
    return;
  }
  
void tkdt_check_tree(kdtom_t *T)
  {
    size_t rec_bytes = kdtom_bytesize(T, FALSE);
    size_t tot_bytes = kdtom_bytesize(T, TRUE);
    fprintf(stderr, "byte sizes: root record = %lu  whole tree = %lu\n", rec_bytes, tot_bytes);
    
    if (T->d == 2) { tkdt_plot(T); }
    return;
  }
    
void tkdt_paint_bullseye(ppv_array_t *A, double ctr[], double R)
  { 
    ppv_dim_t d = A->d;
    ppv_nbits_t bps = A->bps;
    
    /* Bullseye raddi (may be negative): */
    double R3 = R;           /* Outer radius of outer noise ring. */
    double R2 = R3-sqrt(d);  /* Inner radius of outer noise ring. */
    double R1 = 0.5*R2;      /* Outer radius of inner noise ring. */
    double R0 = R1-sqrt(d);  /* Inner radius of inner noise ring. */

    ppv_sample_t maxval = (ppv_sample_t)((1 << bps) - 1);
    assert(maxval > 0);
    ppv_sample_t valA = 0;
    ppv_sample_t valB = (ppv_sample_t)(maxval+1)/2;
    assert(valA != valB);

    auto bool_t bullpaint(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC);
      /* Sets the area within {R0} and outside {R3} to {valA}.
        Sets the area between {R1} and {R2} to {valB}.
        Leaves other voxels unchanged. */
    ppv_enum(bullpaint, FALSE, A, NULL, NULL);

    /* Invert the center sample: */
    ppv_index_t ctrix[d];
    for (ppv_axis_t k = 0; k < d; k++) 
      { ctrix[k] = A->size[k]/2;
        assert(ctrix[k] < A->size[k]);
      }
    ppv_sample_t ctrval = ppv_get_sample(A, ctrix);
    assert(ctrval <= maxval);
    ppv_sample_t newval = (maxval - ctrval < ctrval ? 0 : maxval);
    ppv_set_sample(A, ctrix, newval);
    return;

    /* Internal proc implementations: */

    bool_t bullpaint(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC)
      {
        /* Compute squared distance from pixel center to domain center: */
        double dist2 = 0;
        for (ppv_axis_t k = 0; k < d; k++)
          { double dk = ((double)ix[k]) + 0.5 - ctr[k];
            dist2 += dk*dk;
          }
        double dist = sqrt(dist2);
        if ((dist <= R0) || (dist >= R3))
          { ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pA, valA); }
        else if ((dist >= R1) && (dist <= R2))
          { ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pA, valB); }
        else
          { /* Leave sample unchaged. */ }
        return FALSE;
      }
  }

void tkdt_choose_array_size(ppv_dim_t d, ppv_size_t sz[])
  {
    ppv_sample_count_t npos = 1;
    if (d > 0)
      { /* Choose relative sizes: */
        double rsz_prod = 1;
        double rsz[d];
        for (ppv_axis_t k = 0; k <d; k++)
          { rsz[k] = pow(3.0, ((double)k)/d);
            rsz_prod *= rsz[k];
          }
        double rsz_avg = pow(rsz_prod, 1.0/d);
        /* Compute actual sizes: */
        double xnpos = 1000000.0;         /* Expected number of samples. */
        double xsize = pow(xnpos, 1.0/d); /* Expected size per axis. */
        fprintf(stderr, "array size:");
        for (ppv_axis_t k = 0; k <d; k++)
          { sz[k] = (ppv_size_t)floor(xsize*rsz[k]/rsz_avg);
            fprintf(stderr, " %lu", sz[k]);
            npos *= sz[k];
          }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "array has %lu samples\n", npos);
    return;
  }

tkdt_options_t *tkdt_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    tkdt_options_t *o = (tkdt_options_t *)malloc(sizeof(tkdt_options_t)); 
    
    /* Parse keyword parameters: */
    
    /* argparser_get_keyword(pp, "-size"); */
    /* o->size = (int32_t)argparser_get_next_int(pp, tkdt_size_MIN, tkdt_size_MAX); */

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }

void tkdt_plot(kdtom_t *T)
  {
    demand(T->d == 2, "Can't plot this");
    fprintf(stderr, "Pretend that I am plotting...\n");    
    fprintf(stderr, "*plot*, *plot*, *plit*, *plut*...\n");
    return;
  }
