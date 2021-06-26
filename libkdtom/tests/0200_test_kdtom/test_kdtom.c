#define PROG_NAME "test_kdtom"
#define PROG_DESC "Test the k-d-tree rep of multi-dim sample arrays"
#define PROG_VERS "1.0"

#define test_kdtom_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"

/* Last edited on 2021-06-25 19:50:53 by jstolfi */

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
  "  " test_kdtom_C_COPYRIGHT ".\n" \
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

/* COMMAND-LINE OPTIONS */

typedef struct test_kdtom_options_t
  { int32_t dummy;
  } test_kdtom_options_t;

/* INTERNAL PROTOTYPES */

    /* The unit for all linear dimensions is the voxel side. */ 
    
test_kdtom_options_t *test_kdtom_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {test_kdtom_options_t}. */

void test_kdtom_paint_bullseye(ppv_array_t *A);
  /* Fills {A} with a suitable test pattern. Namely, two spherical
    shells of random values, with constant values between the shells,
    and a different value at the center.  Unless {A} is too small 
    along some axis, in which case the whole array is filled with 
    random values. */

ppv_pos_count_t  test_kdtom_choose_array_size(ppv_dim_t d, ppv_size_t sz[]);
  /* Fills {sz[0..d-1]} with suitable values, not all equal. */

void test_kdtom_do_tests(ppv_dim_t d, ppv_nbits_t bps);

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    test_kdtom_options_t *o = test_kdtom_parse_options(argc, argv);
    
    for (ppv_dim_t d = 0; d <= 16; d++)
      for (ppv_nbits_t bps = 0; bps <= 1;  bps = (5*bps+4)/3)
        { test_kdtom_do_tests(d, bps); }
      
    return 0;
  }
  
void test_kdtom_do_tests(ppv_dim_t d, ppv_nbits_t bps)
  { 
    fprintf(stderr, "Testing d = %u bps = %u ...\n", d, bps);
    ppv_size_t sz[d];
    ppv_pos_count_t npos = test_kdtom_choose_array_size(d, sz);
    ppv_nbits_t bpw = ppv_best_bpw(bps);
    ppv_sample_t maxval = (ppv_sample_t)((1 << bps) - 1);
    ppv_array_t *A = ppv_array_new(d, sz, bps, bpw);
    
    test_kdtom_paint_bullseye(A);
    
    kdtom_array_t *Ta = kdtom_array_make(A);
    /* ??? Check it ??? */
    
    kdtom_t *Tg = kdtom_grind(A);
    /* ??? Check it ??? */
    /* ??? Plot it ??? */

    return;
  }
  
void test_kdtom_paint_bullseye(ppv_array_t *A);
  {
    ppv_dim_t d = A->d;
    ppv_nbits_t bps = A->bps;
    
    if (bps == 0) { /* Array must be all zeros anyway: */ return; }
    
    /* Fill {A} with noise samples: */
    srandom(4615);
    ppv_throw_t(A);
    
    if (d == 0) { /* Only one sample: */ return; }
    
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
    
    if (szmax == 0) { /* Empty aray: */ return; }
        
    if (szmin >= 3) 
      { 
        /* Bullseye raddi (may be negative): */
        double Rmax = 0.50*(double)szmin; /* Radius of max inscribed ball. */
        double R3 = 0.75*Rmax;   /* Outer radius of outer noise ring. */
        double R2 = R3-sqrt(d);  /* Inner radius of outer noise ring. */
        double R1 = 0.5*R2;      /* Outer radius of inner noise ring. */
        double R0 = R1-sqrt(d);  /* Inner radius of inner noise ring. */

        ppv_sample_t maxval = (ppv_sample_t)((1 << bps) - 1);
        assert(maxval > 0);
        ppv_sample_t valA = 0;
        ppv_sample_t valB = (ppv_sample_t)(maxval+1)/2;
        assert(valA != valB);
        auto bool_t bullpaint(const ix[], sample_pos_t pA, sample_pos_t pB, sample_pos_t pC);
          /* Sets the area within {R0} and outside {R3} to {valA}.
            Sets the area between {R1} and {R2} to {valB}.
            Leaves other voxels unchanged. */
        ppv_enum(bullpaint, FALSE, A, NULL, NULL);

        /* Invert the center pixel: */
        ppv_index_t ctrix[d];
        for (ppv_axis_t k = 0; k < d; k++) 
          { ctrix[k] = A->size[k]/2;
            assert(ctrix[k] < A->size[k]);
          }
       ppv_sample_t ctrval = ppv_get_sample(A, ctrix);
       assert(ctrval <= maxval)
       ppv_sample_t newval = (maxval - ctrval < ctrval ? 0 : maxval);
       ppv_set_sample(A, ix, newval);
    }
        
    return;
    
    /* Internal proc implementations: */
    
    bool_t bullpaint(const ix[], sample_pos_t pA, sample_pos_t pB, sample_pos_t pC)
      {
        /* Compute squared distance from pixel center to domain center: */
        double dist2 = 0;
        for (ppv_axis_t k = 0; k < d; k++)
          { double dk = ((double)ix[k]) + 0.5 - ctr[k];
            dist2 += dk*dk;
          }
        double dist = sqrt(dist2);
        if ((dist <= R0) || (dist >= R3))
          { ppv_set_sample_at_pos(A->el, pA, A->bps, A->bpw, valA); }
        else if ((dist >= R1) && (dist <= R2))
          { ppv_set_sample_at_pos(A->el, pA, A->bps, A->bpw, valB); }
        else
          { /* Leave sample unchaged. */ }
        return FALSE;
      }
  }

test_kdtom_options_t *test_kdtom_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    test_kdtom_options_t *o = (test_kdtom_options_t *)malloc(sizeof(test_kdtom_options_t)); 
    
    /* Parse keyword parameters: */
    
    /* argparser_get_keyword(pp, "-size"); */
    /* o->size = (int32_t)argparser_get_next_int(pp, test_kdtom_size_MIN, test_kdtom_size_MAX); */

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }
