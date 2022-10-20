#define PROG_NAME "msm_test_120_mapping"
#define PROG_DESC "test of mapping candidates between scales"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 11:13:13 by stolfi */

#define msm_test_120_mapping_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  {LENGTH_X}  \\\n" \
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
  "  This program generates {N_CANDS} candidate pairings" \
  " between two abstract biosequences with lengths" \
  " {LENGTH_X} and {LENGTH_Y}, reduced by" \
  " a factor {2^M} for some {M}.  It then maps" \
  " those candidates at succesively finer scales, interpolating" \
  " the pairings as needed to keep them atomic.\n" \
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
  "  This program was created on 17/dec/2006 by J. Stolfi.\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " msm_test_120_mapping_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  

#define _GNU_SOURCE
#include <math.h>
#include <stdint.h>
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

#define msm_MAX_MAX_LEVEL 15
  /* Maximum level to consider. */
  
#define msm_MIN_SEQ_SIZE 10
  /* Minimum number of datums in a sequence at coarsest level. */

typedef struct msm_options_t 
  { int32_t lengthX;    /* length of sequence {A}. */
    int32_t lengthY;    /* length of sequence {B}. */
    int32_t nCands;     /* Number of candidates to generate. */
    /* Output parameters: */
    char *outName;  /* Output file name prefix (minus extensions). */
  } msm_options_t;
  
int32_t main(int32_t argc, char**argv);

msm_options_t *msm_get_options(int32_t argc, char**argv);
  /* Parses the command line options, packs 
    them into a {msm_options_t} record. */

msm_cand_vec_t msm_fake_initial_cands
  ( int32_t nCands,
    msm_seq_desc_t *seq0,
    msm_seq_desc_t *seq1
  );
  /* Generates a list of {nCands} candidates between the sequences {seq0,seq1}. */
  
int32_t main(int32_t argc, char**argv)
  { 
    int32_t M = msm_MAX_MAX_LEVEL;
    
    msm_options_t *o = msm_get_options(argc, argv);

    fprintf(stderr, "computing sequence descriptors for scales %d .. %d ...\n", 0, M-1);
    msm_seq_desc_t sv0[M+1];
    msm_seq_desc_t sv1[M+1];
    int32_t level;
    for (level = 0; level <= M; level++) 
      { if (level == 0)
          { /* Original sequences: */
            sv0[level] = msm_seq_desc_make(100, "X", FALSE, o->lengthX, 0, 0);
            sv1[level] = msm_seq_desc_make(101, "Y", TRUE,  o->lengthY, 0, 0);
          }
        else
          { /* Filtered sequences: */
            /* Choose an arbitrary kernel radius {rwtb}: */
            int32_t rwtb = 3 + (level % 3);
            /* Choose an arbitrary resampling order {ek}: */
            int8_t ek = (int8_t)(level == 1 ? -2 : ((level % 4) == 2 ? -1 : +1));
            /* Simulate the filtering: */
            fprintf(stderr, "simulating filtering from level %2d to level %2d:", level-1, level);
            fprintf(stderr, "  kernel radius = %3d  resampling order = %+3d\n\n", rwtb, ek);
            /* Trim ends where filter don't fit: */
            msm_seq_desc_t tmp0 = msm_seq_desc_trim(&(sv0[level-1]), rwtb, rwtb);
            msm_seq_desc_t tmp1 = msm_seq_desc_trim(&(sv1[level-1]), rwtb, rwtb);
            fprintf(stderr, "    tmp 0 = "); msm_seq_desc_write(stderr, "(", &tmp0, 4, 1, 6, ")\n");
            fprintf(stderr, "    tmp 1 = "); msm_seq_desc_write(stderr, "(", &tmp1, 4, 1, 6, ")\n\n");
            /* Virtual resampling: */
            sv0[level] = msm_seq_desc_resample(&(tmp0), ek);
            sv1[level] = msm_seq_desc_resample(&(tmp1), ek);
          }
        msm_seq_desc_t *seq0 = &(sv0[level]);
        msm_seq_desc_t *seq1 = &(sv1[level]);
        fprintf(stderr, "  seq 0 = "); msm_seq_desc_write(stderr, "(", seq0, 4, 1, 6, ")\n");
        fprintf(stderr, "  seq 1 = "); msm_seq_desc_write(stderr, "(", seq1, 4, 1, 6, ")\n\n");
        fprintf(stderr, "level %2d  seq lengths %3d %3d", level, seq0->size, seq1->size);
        fprintf(stderr, " esteps %+3d %+3d", seq0->estep, seq1->estep);
        fprintf(stderr, " skips %3d %3d", seq0->skip, seq1->skip);
        fprintf(stderr, "\n\n");
      }

    fprintf(stderr, "reducing max level to ensure nonzero strings...\n");
    while ((M > 0) && ((sv0[M].size < msm_MIN_SEQ_SIZE) || (sv1[M].size < msm_MIN_SEQ_SIZE))) { M--; }
    fprintf(stderr, "  max useful level = %d\n", M);
    
    /* The current candidate list: */
    msm_cand_vec_t cdv = msm_cand_vec_new(0);

    fprintf(stderr, "multiscale mapping ...\n");
    for (level = M; level >= 0; level--)
      { 
        /* Scale reduction factor from original strings: */
        int32_t scale = 1 << level;

        fprintf(stderr, "\n");
        fprintf(stderr, "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n");
        fprintf(stderr, "LEVEL %d - SCALE %d\n", level, scale);

        /* Descriptors for the two sequences at the current level: */
        msm_seq_desc_t *seq0 = &(sv0[level]);
        msm_seq_desc_t *seq1 = &(sv1[level]);

        /* Make filename prefix for current level: */
        char *levName = NULL;
        asprintf(&levName, "%s-%02d", o->outName, level);
        
        if (level == M)
          { /* Set {cdv} to an arbitrary set of initial candidates: */
            fprintf(stderr, "finding initial candidates at level %2d ...\n", level);
            fprintf(stderr, "  seq 0 = "); msm_seq_desc_write(stderr, "(", seq0, 4, 1, 6, ")\n");
            fprintf(stderr, "  seq 1 = "); msm_seq_desc_write(stderr, "(", seq1, 4, 1, 6, ")\n");
            cdv = msm_fake_initial_cands(o->nCands, seq0, seq1);
          }
        else
          { /* Map the candidates in {cdv} to the next finer scale,yield {cdm}: */
            msm_seq_desc_t *ant0 = &(sv0[level+1]);
            msm_seq_desc_t *ant1 = &(sv1[level+1]);
            fprintf(stderr, "mapping from level %2d to level %2d ...\n", level+1, level);
            fprintf(stderr, "  seq 0 = "); 
            msm_seq_desc_write(stderr, "(", ant0, 4, 1, 6, ")");
            fprintf(stderr, " --> "); 
            msm_seq_desc_write(stderr, "(", seq0, 4, 1, 6, ")\n");
            fprintf(stderr, "  seq 1 = "); 
            msm_seq_desc_write(stderr, "(", ant1, 4, 1, 6, ")");
            fprintf(stderr, " --> "); 
            msm_seq_desc_write(stderr, "(", seq1, 4, 1, 6, ")\n");
            msm_cand_vec_t cdm = msm_cand_vec_map(&cdv, seq0, seq1);
            msm_cand_vec_free(&cdv);

            /* Write and plot the mapped candidates {cdm} for this level: */
            fprintf(stderr, "writing mapped candidates ...\n");
            msm_cand_vec_write_named(&cdm, levName, "-m", ".cdv");
            msm_image_cand_vec_write_named(&cdm, seq0, seq1, TRUE, levName, "-m");

            /* Make the candidates strictly increasing on both sides, yielding {cdi}: */
            msm_cand_vec_t cdi = msm_cand_vec_make_increasing(&cdm, 1, 1);
            msm_cand_vec_free(&cdm);
            
            /* Eliminate steps that increase 2 or more on both sides, save in {cdv}: */
            cdv = msm_cand_vec_interpolate(&cdi);
            msm_cand_vec_free(&cdi);
          }
          
        /* Write and plot the official candidates for this level: */
        fprintf(stderr, "writing official candidates ...\n");
        msm_cand_vec_write_named(&cdv, levName, "-v", ".cdv");
        fprintf(stderr, "plotting the official candidates ...\n");
        msm_image_cand_vec_write_named(&cdv, seq0, seq1, TRUE, levName, "-v");

        fprintf(stderr, "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n");
        fprintf(stderr, "\n");
      }

    return 0;
  }

msm_cand_vec_t msm_fake_initial_cands
  ( int32_t nCands,
    msm_seq_desc_t *seq0,
    msm_seq_desc_t *seq1
  )
  { 
    /* Define the max pairing length {maxlen}: */
    int32_t n0 = seq0->size;
    int32_t n1 = seq1->size;
    int32_t maxlen = (n0 < n1 ? n0 : n1);
    int32_t minlen = (maxlen + 9)/10;
    /* Generate the candidates: */
    msm_cand_vec_t cdv = msm_cand_vec_new(nCands);
    int32_t ncd = 0; /* Start storing candidates here. */
    msm_cand_vec_throw
      ( nCands,
        seq0, seq1,
        minlen, maxlen,
        /*unitProb*/ 0.2,
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
