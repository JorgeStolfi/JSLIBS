#define PROG_NAME "test_dnae_seqio"
#define PROG_DESC "test of DNA signal reading/writing routines"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-31 09:33:29 by stolfi */

#define test_dnae_seqio_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  {SEQ_FILE} \\\n" \
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
  " DNA/RNA sequence from {SEQ_FILE}, converts it to" \
  " numeric format with subsampling factor {NSUB}, and" \
  " outputs the numeric version of that" \
  " sequence.  The program then reads" \
  " that file and checks whether the data is preserved.  The" \
  " input file should be in \".seq\" format" \
  " (see the INPUT FILE FORMAT section below).\n" \
  "\n" \
  "INPUT FILE FORMAT\n" \
  dnae_nucleic_file_format_INFO "\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  All output files will have names starting with {OUT_NAME}.\n" \
  "\n" \
  "OPTIONS\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  dm_match(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created on 16/dec/2006 by J. Stolfi.\n" \
  "MODIFICATION HISTORY\n" \
  "  2014-06-09 J. Stolfi: renamed {test_dnae_coding}.\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_dnae_seqio_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>

#include <dnae_seq.h>
#include <dnae_nucleic.h>

typedef struct options_t 
  { char *seqName;     /* Name of file with DNA/RNA sequence (minus ".bas" ext). */
    char *outName;     /* Output file name prefix (minus extensions). */
  } options_t;
  
int32_t main(int32_t argc, char**argv);

options_t *dm_get_options(int32_t argc, char**argv);
  /* Parses the command line options, packs 
    them into a {options_t} record. */

int32_t main(int32_t argc, char**argv)
  { 
    options_t *o = dm_get_options(argc, argv);

    /* Arbitrary header data: */
    msm_seq_id_t id = 37;
    char *name = "hxxx";
    bool_t rev = FALSE;

    fprintf(stderr, "reading input DNA/RNA sequence ...\n");
    dnae_seq_t seq = dnae_seq_read_from_nucleic_file_named(o->seqName, "", ".bas", id, name, rev);

    /* Create output filename: */
    fprintf(stderr, "writing numeric DNA/RNA sequence...\n");
    FILE *wr = open_write_tag_ext(o->outName, "", ".eqs", TRUE);
    dnae_seq_write(wr, &seq);
    fclose(wr);
    
    fprintf(stderr, "reading numeric DNA/RNA sequence back in...\n");
    FILE *rd = open_read_tag_ext(o->outName, "", ".eqs", TRUE);
    dnae_seq_t srd = dnae_seq_read(rd);
    fclose(rd);
    
    fprintf(stderr, "checking...\n");
    assert(srd.sd.id == seq.sd.id);
    assert(strcmp(srd.sd.name, seq.sd.name) == 0);
    assert(srd.sd.rev == seq.sd.rev);
    assert(srd.sd.estep == seq.sd.estep);
    assert(srd.sd.skip == seq.sd.skip);
    assert(strcmp(srd.cmt, seq.cmt) == 0);
    assert(srd.dv.ne == seq.dv.ne);
    int32_t k;
    for (k = 0; k < dnae_CHANNELS; k++)
      { assert(seq.sfac.f[k] == srd.sfac.f[k]); 
        int32_t i;
        for (i = 0; i < seq.dv.ne; i++)
          { dnae_sample_enc_t ssmp = dnae_seq_get_sample_enc(&seq, i, k); 
            dnae_sample_enc_t tsmp = dnae_seq_get_sample_enc(&srd, i, k);
            assert(ssmp == tsmp);
          }
      }
    
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
    
    argparser_skip_parsed(pp);
    
    o->seqName = argparser_get_next(pp);
    
    o->outName = argparser_get_next(pp);

    argparser_finish(pp);
    
    return o;
  }
