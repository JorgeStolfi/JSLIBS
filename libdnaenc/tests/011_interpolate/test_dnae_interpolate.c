#define PROG_NAME "test_dnae_interpolate"
#define PROG_DESC "test of DNA signal interpolating routines"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-31 14:38:17 by stolfi */

#define test_dnae_interpolate_C_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -step {STEP} \\\n" \
  "  {SEQ_FILE_NAME} \\\n" \
  "  {OUT_NAME}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program reads a" \
  " DNA/RNA sequence from the file named \"{SEQ_FILE_NAME}.eqs\", converts it" \
  " to numeric format, and outputs interpolated versions of that" \
  " sequence.  The input file should be in the \".eqs\" format.\n" \
  "\n" \
  "INPUT FILE FORMAT\n" \
  "  " dnae_seq_file_format_INFO "\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  All output files will have names starting with {OUT_NAME}.\n" \
  "\n" \
  "OPTIONS\n" \
  "\n" \
  "  -step {STEP} \n" \
  "    This mandatory argument specifies the interpolation step, which should be a negative power of 2.\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  dm_match(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created on 2014-07-16 by J. Stolfi.\n" \
  "MODIFICATION HISTORY\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_dnae_interpolate_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <argparser.h>
#include <affirm.h>
#include <vec.h>
#include <wt_table.h>

#include <msm_multi.h>
#include <msm_test_tools.h>

#include <dnae_nucleic.h>
#include <dnae_seq.h>
#include <dnae_test_tools.h>

typedef struct options_t 
  { char *seqFileName; /* Name of file with DNA/RNA sequence (minus ".bas" ext). */
    double step;       /* Interpolation step. */
    char *outName;     /* Output file name prefix (minus extensions). */
  } options_t;
  
int32_t main(int32_t argc, char**argv);

options_t *dm_get_options(int32_t argc, char**argv);
  /* Parses the command line options, packs 
    them into a {options_t} record. */

int32_t main(int32_t argc, char**argv)
  { 
    options_t *o = dm_get_options(argc, argv);

    fprintf(stderr, "reading input sequence ...\n");
    dnae_seq_t s = dnae_seq_read_named(o->seqFileName, "", ".eqs");
    fprintf(stderr, "sequence has %d sample datums\n", dnae_seq_num_datums(&s));
    
    fprintf(stderr, "trying to interpolate with step %8.6f ...\n", o->step);
    double step = o->step;
    int8_t ek = 0;
    while (step <= 0.5) { ek--; step *= 2; }
    demand(fabs(step - 1.0) < 1.0e-14, "step must be a non-positive power of 2");
    dnae_seq_t r = dnae_seq_interpolate(&s, ek);
    
    fprintf(stderr, "writing and plotting the interpolated sequence ...\n");
    bool_t plot = TRUE;
    double hSize = 180.0; /* mm */
    double vSize = 60.0; /* mm */
    double fontSize = 8.0; /* pt */
    dnae_test_tools_seq_write_and_plot_named
      ( &r, "Interpolated Sequence", o->outName, "-ot", 
        plot, hSize, vSize, fontSize, -1
      );

    /* ...and we are done: */
    fprintf(stderr, "done.\n");
    return 0;
  }

options_t* dm_get_options(int32_t argc, char**argv)
  { options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    argparser_get_keyword(pp, "-step");
    o->step = argparser_get_next_double(pp, 1.0/256.0, 1.000);

    argparser_skip_parsed(pp);
    
    o->seqFileName = argparser_get_next(pp);
    
    o->outName = argparser_get_next(pp);

    argparser_finish(pp);
    
    return o;
  }
