#define PROG_NAME "msm_test_110_cands"
#define PROG_DESC "test of pairing tools and file formats"
#define PROG_VERS "1.0"

/* Last edited on 2008-01-28 17:26:33 by hcgl */

#define msm_test_110_cands_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  {LENGTH_X} {CIRC_X} \\\n" \
  "  {LENGTH_Y} {CIRC_Y} \\\n" \
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
  " them.  Then reads the disck file back in, and compares" \
  " its contents with the original list.\n" \
  "\n" \
  "  The boolean arguments {CIRC_X} and {CIRC_Y} (\"F\" for FALSE, \"T\" for TRUE)" \
  " indicate whether the corresponding sequence is circular or open." \
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
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include <affirm.h>
#include <float_image.h>
#include <float_pnm_image.h>
#include <jspnm_image.h>
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
  { int lengthX;    /* length of X sequence. */
    bool_t circX;   /* TRUE iff X sequence is circular. */
    int lengthY;    /* length of Y sequence. */
    bool_t circY;   /* TRUE iff Y sequence is circular. */
    int nCands;     /* Number of candidates to generate. */
    /* Output parameters: */
    char *outName;  /* Output file name prefix (minus extensions). */
  } msm_options_t;
  
int main(int argc, char**argv);

msm_options_t *msm_get_options(int argc, char**argv);
  /* Parses the command line options, packs 
    them into a {msm_options_t} record. */

msm_cand_vec_t msm_generate_cands(msm_seq_desc_t *xp, msm_seq_desc_t *yp, msm_options_t *o);
  /* Generates a list of {o->nCands} pairings for two abstract
    sequences {xp,yp}. */
  
void msm_compare_cand_vecs(msm_cand_vec_t *cdva, msm_cand_vec_t *cdvb);
  /* Checks whether the candidate lists {cdva} and {cdvb} contains the same
    candidates, with the same scores, in the same order.  */

int main(int argc, char**argv)
  { 
    msm_options_t *o = msm_get_options(argc, argv);
    
    /* Make descriptors for the two sequences: */
    int level = 3; /* Arbitrary. */
    msm_seq_desc_t x = msm_seq_desc_make( 0, "x", level, o->lengthX, o->circX );
    msm_seq_desc_t y = msm_seq_desc_make( 1, "y", level, o->lengthY, o->circY );
    
    fprintf(stderr, "generating cands ...\n");
    msm_cand_vec_t cdv = msm_generate_cands(&x, &y, o);

    /* Write and plot candidates: */
    fprintf(stderr, "writing cands ...\n");
    msm_cand_vec_write_named(&cdv, o->outName, "-cd");
    msm_image_cand_vec_write_named(&cdv, &x, &y, o->outName, "-cd");

    /* Read the candidates back in and compare: */
    fprintf(stderr, "re-reading cands ...\n");
    msm_cand_vec_t cdr = msm_cand_vec_read_named(o->outName, "-cd");
    msm_compare_cand_vecs(&cdv, &cdr);
    
    return 0;
  }
    
void msm_compare_cand_vecs(msm_cand_vec_t *cdva, msm_cand_vec_t *cdvb)
  { demand(cdva->ne == cdvb->ne, "candidate counts differ"); 
    int ncd = cdva->ne;
    int k;
    for (k = 0; k < ncd; k++)
      { msm_cand_t *cda = &(cdva->e[k]);
        msm_cand_t *cdb = &(cdvb->e[k]);
        msm_cand_equivalent(cda, cdb, TRUE);
      }
  }

msm_cand_vec_t msm_generate_cands(msm_seq_desc_t *xp, msm_seq_desc_t *yp, msm_options_t *o)
  { /* Define the max candidate length {maxlen}: */
    int maxlen = (xp->npos > yp->npos ? xp->npos : yp->npos);
    int minlen = (maxlen + 9)/10;
    /* Generate the candidates: */
    msm_cand_vec_t cdv = msm_cand_vec_throw
      ( o->nCands,
        xp, 
        yp,
        minlen, maxlen,
        /*circProb*/ 0.1, 
        /*atomProb*/ 0.5, 
        /*diagProb*/ 0.2,
        /*skipProb*/ 0.1
      );
    return cdv;
  }

msm_options_t *msm_get_options(int argc, char**argv)
  { 
    msm_options_t *o = (msm_options_t *)notnull(malloc(sizeof(msm_options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    argparser_skip_parsed(pp);
    
    o->lengthX = argparser_get_next_int(pp, 0, INT_MAX);
    o->circX = argparser_get_next_bool(pp);
    
    o->lengthY = argparser_get_next_int(pp, 0, INT_MAX);
    o->circY = argparser_get_next_bool(pp);
    
    o->nCands = argparser_get_next_int(pp, 0, INT_MAX);
    o->outName = argparser_get_next(pp);

    argparser_finish(pp);
    
    return o;
  }
