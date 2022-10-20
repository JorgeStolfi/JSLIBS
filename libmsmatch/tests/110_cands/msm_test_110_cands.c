#define PROG_NAME "msm_test_110_cands"
#define PROG_DESC "test of pairing tools and file formats"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 11:13:03 by stolfi */

#define msm_test_110_cands_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  {LENGTH_X} \\\n" \
  "  {LENGTH_Y}  \\\n" \
  "  {N_CANDS} \\\n" \
  "  {OUT_NAME}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program generates {N_CANDS} random candidates" \
  " between two abstract sequences with lengths" \
  " {LENGTH_X} and {LENGTH_Y}, writes them to a disk file, and plots" \
  " them.  Then reads the disk file back in, and compares" \
  " its contents with the original list.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  All output files will have names starting with {OUT_NAME}.\n" \
  "\n" \
  "OPTIONS\n" \
  "  None.\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  msm_match(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created on 29/dec/2006 by J. Stolfi.\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " msm_test_110_cands_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include <affirm.h>
#include <float_image.h>
#include <uint16_image.h>
#include <argparser.h>

#include <msm_basic.h>
#include <msm_seq_desc.h>
#include <msm_pairing.h>
#include <msm_cand.h>
#include <msm_cand_vec.h>
#include <msm_image.h>
#include <msm_image_tools.h>

#define msm_MAX_LEVELS 20
  /* Maximum levels to consider above level 0. */
  
#define msm_MIN_SEQ_SIZE 20
  /* Minimum number of datums in a sequence at coarsest level. */

typedef struct msm_options_t 
  { int32_t lengthX;    /* length of X sequence. */
    int32_t lengthY;    /* length of Y sequence. */
    int32_t nCands;     /* Number of candidates to generate. */
    /* Output parameters: */
    char *outName;  /* Output file name prefix (minus extensions). */
  } msm_options_t;
  
int32_t main(int32_t argc, char**argv);

msm_options_t *msm_get_options(int32_t argc, char**argv);
  /* Parses the command line options, packs 
    them into a {msm_options_t} record. */

msm_cand_vec_t msm_generate_cands(msm_seq_desc_t *xp, msm_seq_desc_t *yp, msm_options_t *o);
  /* Generates a list of {o->nCands} pairings for two abstract
    sequences {xp,yp}. */
  
void msm_compare_cand_vecs(msm_cand_vec_t *cdva, msm_cand_vec_t *cdvb);
  /* Checks whether the candidate lists {cdva} and {cdvb} contains the same
    candidates, with the same scores, in the same order.  */

int32_t main(int32_t argc, char**argv)
  { 
    msm_options_t *o = msm_get_options(argc, argv);
    
    /* Make descriptors for the two sequences: */
    int32_t wrad0 = 5; /* Filter kernel radius for level 0-->1 (arbitrary). */
    int32_t wrad1 = 3; /* Filter kernel radius for higher levels (arbitrary). */
    int8_t estep = 3; /* Arbitrary. */
    int32_t skip = wrad0 + wrad1*((1 << estep) - 1); /* Skip for {level}. */
    msm_seq_desc_t x = msm_seq_desc_make( 0, "x", FALSE, o->lengthX - 2*skip, estep, skip);
    msm_seq_desc_t y = msm_seq_desc_make( 1, "y", TRUE,  o->lengthY - 2*skip, estep, skip);
    
    fprintf(stderr, "generating cands ...\n");
    msm_cand_vec_t cdv = msm_generate_cands(&x, &y, o);

    /* Write and plot candidates: */
    fprintf(stderr, "writing cands ...\n");
    msm_cand_vec_write_named(&cdv, o->outName, "-cd", ".cdv");
    msm_image_cand_vec_write_named(&cdv, &x, &y,TRUE, o->outName, "-cd");

    /* Read the candidates back in and compare: */
    fprintf(stderr, "re-reading cands ...\n");
    msm_cand_vec_t cdr = msm_cand_vec_read_named(o->outName, "-cd", ".cdv");
    msm_compare_cand_vecs(&cdv, &cdr);
    
    return 0;
  }
    
void msm_compare_cand_vecs(msm_cand_vec_t *cdva, msm_cand_vec_t *cdvb)
  { demand(cdva->ne == cdvb->ne, "candidate counts differ"); 
    int32_t ncd = cdva->ne;
    int32_t k;
    for (k = 0; k < ncd; k++)
      { msm_cand_t *cda = &(cdva->e[k]);
        msm_cand_t *cdb = &(cdvb->e[k]);
        msm_cand_equivalent(cda, cdb, TRUE);
      }
  }

msm_cand_vec_t msm_generate_cands(msm_seq_desc_t *xp, msm_seq_desc_t *yp, msm_options_t *o)
  { /* Define the max candidate length {maxlen}: */
    int32_t maxlen = (xp->size < yp->size ? xp->size : yp->size);
    int32_t minlen = (maxlen + 9)/10;
    /* Generate the candidates: */
    msm_cand_vec_t cdv = msm_cand_vec_new(o->nCands);
    int32_t ncd = 0; /* Start storing candidates here. */
    msm_cand_vec_throw
      ( o->nCands,
        xp, 
        yp,
        minlen, maxlen,
        /*atomProb*/ 0.5, 
        /*diagProb*/ 0.2,
        /*skipProb*/ 0.1,
        &cdv,
        &ncd
      );
    msm_cand_vec_trim(&cdv, ncd);
    return cdv;
  }

msm_options_t *msm_get_options(int32_t argc, char**argv)
  { 
    msm_options_t *o = (msm_options_t *)notnull(malloc(sizeof(msm_options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    argparser_skip_parsed(pp);
    
    o->lengthX = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX);
    
    o->lengthY = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX);
    
    o->nCands = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX);
    o->outName = argparser_get_next(pp);

    argparser_finish(pp);
    
    return o;
  }
