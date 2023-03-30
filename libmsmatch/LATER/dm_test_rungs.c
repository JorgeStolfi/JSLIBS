#define PROG_NAME "dm_test_rungs"
#define PROG_DESC "test of procedures that deal with rungs"
#define PROG_VERS "1.0"

/* Last edited on 2023-03-18 11:36:02 by stolfi */

#define dm_test_rungs_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  {LENGTH_X} {CIRC_X} \\\n" \
  "  {LENGTH_Y} {CIRC_Y} \\\n" \
  "  {N_TESTS} \\\n" \
  "  {OUT_NAME}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program bla bla bla.\n" \
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
  "  dm_match(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created on 11/dec/2007 by J. Stolfi.\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " dm_test_rungs_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  

#include <dm_basic.h>
#include <stdint.h>

#include <msm_image.h>
#include <msm_rung.h>

#include <affirm.h>
#include <float_image.h>
#include <uint16_image.h>
#include <argparser.h>

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#define dm_MAX_LEVELS 20
  /* Maximum levels to consider above level 0. */
  
#define dm_MIN_SEQ_SIZE 20
  /* Minimum number of datums in a sequence at coarsest level. */

typedef struct options_t 
  { int32_t lengthX;       /* length of X sequence. */
    bool_t circX;      /* TRUE iff X sequence is circular. */
    int32_t lengthY;       /* length of Y sequence. */
    bool_t circY;      /* TRUE iff Y sequence is circular. */
    /* Output parameters: */
    char *outName;     /* Output file name prefix (minus extensions). */
  } options_t;
  
int32_t main(int32_t argc, char**argv);

options_t *msm_get_options(int32_t argc, char**argv);
  /* Parses the command line options, packs 
    them into a {options_t} record. */

dm_cand_vec_t msm_generate_cands(options_t *o);
  /* Generates a list of {o->nPairings} candidates for the two
    {dm_seq_t}s. Assumes that the first sequence has {o->lengthX}
    elements, and is circular iff {o->circX} is TRUE. Ditto for
    {o->lengthY,o->circY}. */
  
int32_t main(int32_t argc, char**argv)
  { 
    options_t *o = msm_get_options(argc, argv);
    
    fprintf(stderr, "generating candidates ...\n");
    dm_cand_vec_t cdv = msm_generate_cands(o);

    /* Write and plot candidates: */
    fprintf(stderr, "writing candidates ...\n");
    msm_cand_vec_write_named(&cdv, NULL, "-cd");
    msm_cand_vec_image_write_named(&cdv, o->lengthX, o->circX, o->lengthY, o->circY, NULL, "-cd");

    return 0;
  }
    
dm_cand_vec_t msm_generate_cands(options_t *o)
  { /* Define the max candidate length {maxlen}: */
    int32_t nx = o->lengthX;
    int32_t ny = o->lengthY;
    int32_t maxlen = (nx > ny ? nx : ny);
    int32_t minlen = (maxlen + 9)/10;
    /* Generate the candidates: */
    dm_cand_vec_t cdv = msm_cand_vec_throw
      ( o->nPairings,
        0, "a", nx, o->circX, 
        1, "b", ny, o->circY,
        minlen, maxlen,
        /*circProb*/ 0.1, 
        /*atomProb*/ 0.5, 
        /*diagProb*/ 0.2,
        /*skipProb*/ 0.1
      );
    return cdv;
  }

options_t *msm_get_options(int32_t argc, char**argv)
  { 
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    argparser_skip_parsed(pp);
    
    o->lengthX = argparser_get_next_int(pp, 0, INT32_MAX);
    o->circX = argparser_get_next_bool(pp);
    
    o->lengthY = argparser_get_next_int(pp, 0, INT32_MAX);
    o->circY = argparser_get_next_bool(pp);
    
    o->nPairings = argparser_get_next_int(pp, 0, INT32_MAX);
    o->outName = argparser_get_next(pp);

    argparser_finish(pp);
    
    return o;
  }
