#define PROG_NAME "test_dnae_coding"
#define PROG_DESC "test of DNA signal encoding/decoding routines"
#define PROG_VERS "1.0"

/* Last edited on 2014-06-14 02:52:28 by stolfilocal */

#define test_dnae_coding_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " {NSAMPLES} {OUT_NAME}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program checks the signal encoding and decoding" \
  " routines with {NSAMPLES} samples.\n" \
  "\n" \
  "  All output files will have names starting with {OUT_NAME}.\n" \
  "\n" \
  "OPTIONS\n" \
  "  None.\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  dm_test_filter(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created on 16/dec/2006 by J. Stolfi as {dm_test_005_coding}.\n" \
  "MODIFICATION HISTORY\n" \
  "  2014-06-09 J. Stolfi: renamed {test_dnae_coding}.\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_dnae_coding_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <jsrandom.h>
#include <jsfile.h>
#include <affirm.h>
#include <argparser.h>

#include <dnae_sample.h>
#include <dnae_seq.h>

typedef struct options_t 
  { int nSamples;      /* Average number of test samples per bin. */
  } options_t;
  
int main(int argc, char**argv);

options_t *dm_get_options(int argc, char**argv);
  /* Parses the command line options, packs 
    them into a {options_t} record. */

int main(int argc, char**argv)
  { 
    options_t *o = dm_get_options(argc, argv);
    int NS = o->nSamples;
    
    double scale = 1.00; 
    
    fprintf(stderr, "  code value        recodings\n");
    fprintf(stderr, "------ ------------ --------------\n");
    
    /* Encoding consistency check: */
    int i;
    for (i = dnae_sample_enc_VALID_MIN; i < dnae_sample_enc_VALID_MAX; i++)
      { dnae_sample_enc_t ev = (dnae_sample_enc_t)i;
        double fv = dnae_sample_decode(ev, scale);
        dnae_sample_enc_t tva = dnae_sample_encode(fv - 1.0e-8*(1 + fabs(fv)), scale);
        dnae_sample_enc_t tvb = dnae_sample_encode(fv, scale);
        dnae_sample_enc_t tvc = dnae_sample_encode(fv + 1.0e-8*(1 + fabs(fv)), scale);
        fprintf(stderr, "%+06d %12.8f %+06d %+06d %+06d\n", (int)ev, fv, (int)tva, (int)tvb, (int)tvc);
        assert((ev == tva) && (ev == tvb) && (ev == tvc));
      }
    
    /* Encoding error test: */
    double sumd = 0;
    double sumd2 = 0;
    int NV = dnae_sample_enc_VALID_MAX - dnae_sample_enc_VALID_MIN + 1;
    int hist[NV];
    double hsumd[NV];
    double hsumd2[NV];
    for (i = 0; i < NV; i++) { hist[i] = 0; hsumd[i] = 0; hsumd2[i] = 0; }
    int k;
    for (k = 0; k < NS*NV; k++)
      { /* Generate a random normal-distributed value {dv}: */
        double rv = dgaussrand();
        /* Encode it: */
        dnae_sample_enc_t ev = dnae_sample_encode(rv, scale);
        int i = ((int)ev) - ((int)dnae_sample_enc_VALID_MIN);
        assert((i >= 0) && (i < NV));
        /* Decode the encoded value: */
        double dv = dnae_sample_decode(ev, scale);
        /* Compute the encoding/decoding error {dif}: */
        double dif = dv - rv;
        /* Accumulate the value and squared value of error: */
        sumd += dif;
        sumd2 += dif*dif;
        /* Histogram count and accumulate error in bin: */
        hist[i] += 1;
        hsumd[i] += dif;
        hsumd2[i] += dif*dif;
      }
    /* Compute and print the average and root-mean-square error: */
    double avg_err = sumd/NS;
    double rms_err = sqrt(sumd2/NS);
    fprintf(stderr, "avg error = %+11.7f\n", avg_err);
    fprintf(stderr, "rms error = %11.7f\n", rms_err);
    
    /* Print to stdout the encoded value histogram and errors per bin: */
    fprintf(stdout, "#   code        value count    prob   avg error   rms error\n");
    fprintf(stdout, "# ------ ------------ ----- ------- ----------- -----------\n");
    for (i = 0; i < NV; i++)
      { dnae_sample_enc_t ev = (dnae_sample_enc_t)(i + dnae_sample_enc_VALID_MIN);
        double dv = dnae_sample_decode(ev, scale);
        double prob_h = ((double)(hist[i]))/NS;
        double avg_err_h = hsumd[i]/hist[i];
        double rms_err_h = sqrt(hsumd2[i]/hist[i]);
        fprintf(stdout, "  %+06d %12.8f %5d %7.5f %+11.7f %11.7f\n", ev, dv, hist[i], prob_h, avg_err_h, rms_err_h);
      }
    fflush(stdout);

    /* ...and we are done: */
    fprintf(stderr, "done.\n");
    return 0;
  }
  
options_t* dm_get_options(int argc, char**argv)
  { options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    argparser_skip_parsed(pp);
    
    o->nSamples = (int)argparser_get_next_int(pp, 1, INT_MAX);

    argparser_finish(pp);
    
    return o;
  }
