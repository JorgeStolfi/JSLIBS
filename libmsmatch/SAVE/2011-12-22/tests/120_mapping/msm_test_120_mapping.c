#define PROG_NAME "msm_test_120_mapping"
#define PROG_DESC "test of mapping candidates between scales"
#define PROG_VERS "1.0"

/* Last edited on 2008-02-01 12:48:34 by hcgl */

#define msm_test_120_mapping_C_COPYRIGHT \
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
  "  This program generates {N_CANDS} candidate pairings" \
  " between two abstract biosequences with lengths" \
  " {LENGTH_X} and {LENGTH_Y}, reduced by" \
  " a factor {2^M} for some {M}.  It then maps" \
  " those candidates at succesively finer scales, interpolating" \
  " the pairings as needed to keep them atomic.\n" \
  "\n" \
  "  The boolean arguments {CIRC_X} and {CIRC_Y}" \
  " (\"F\" for FALSE, \"T\" for TRUE)" \
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
  { int lengthX;    /* length of sequence {A}. */
    bool_t circX;   /* TRUE iff sequence {A} is circular. */
    int lengthY;    /* length of sequence {B}. */
    bool_t circY;   /* TRUE iff sequence {B} is circular. */
    int nCands;     /* Number of candidates to generate. */
    /* Output parameters: */
    char *outName;  /* Output file name prefix (minus extensions). */
  } msm_options_t;
  
int main(int argc, char**argv);

msm_options_t *msm_get_options(int argc, char**argv);
  /* Parses the command line options, packs 
    them into a {msm_options_t} record. */

void msm_compute_multiscale_lengths(int ns, bool_t circ, int maxLevel, int nwt[], int nfs[]);
 /* Computes the lengths of sequences produces by {msm_MAX_LEVELS}
   successive applications of filtering and downsampling on a sequence
   with {ns} elements.
   
   Assumes that the fitering from level {r} to level {r+1} uses a
   weight table with {nwt[r]} entries. If {circ} is TRUE, assumes that
   the original sequence and all its refined versions are circular.
   
   The computed lengths are returned in {nfs[0..msm_MAX_LEVELS]}, with
   {nfs[0] = ns}. */

msm_cand_vec_t msm_fake_initial_cands
  ( int nCands,
    int level,
    msm_seq_desc_t *xp,
    msm_seq_desc_t *yp
  );
  /* Generates a list of {nCands} candidates between the sequences
    {xp,yp}, which must have the same {level}. */

int msm_seq_desc_filtered_size(int nFine, bool_t circ, int nwtb);
  /* Computes the number of distinct matching positions in the
    filtered and downsampled version of a sequence with {nFine}
    positions and circularity {circ}. */
  
int main(int argc, char**argv)
  { 
    int M = msm_MAX_LEVELS;
    
    msm_options_t *o = msm_get_options(argc, argv);

    fprintf(stderr, "choosing the filter window size for scales %d .. %d ...\n", 0, M-1);
    int nwt[M]; /* Number of weights in filter kernel at each scale. */
    int level;
    for (level = 0; level < M; level++) { nwt[level] = 7 + 2*(level % 3); }
    
    fprintf(stderr, "computing seq lengths at scales %d .. %d ...\n", 0, M);
    int nfx[M + 1]; msm_compute_multiscale_lengths(o->lengthX, o->circX, M, nwt, nfx);
    int nfy[M + 1]; msm_compute_multiscale_lengths(o->lengthY, o->circY, M, nwt, nfy);
    for (level = 0; level <= M; level++)
      { int nx = nfx[level], ny = nfy[level]; 
        fprintf(stderr, "  level %2d  seq lengths %3d %3d", level, nx, ny);
        if (level < M) 
          { int nw = nwt[level];
            fprintf(stderr, " filter window width %3d", nw);
          }
        fprintf(stderr, "\n");
      }

    fprintf(stderr, "reducing max level to ensure nonzero strings...\n");
    while ((M > 0) && ((nfx[M] < msm_MIN_SEQ_SIZE) || (nfy[M] < msm_MIN_SEQ_SIZE))) { M--; }
    fprintf(stderr, "  max level = %d\n", M);
    
    /* The current candidate list: */
    msm_cand_vec_t cdv = msm_cand_vec_new(0);

    fprintf(stderr, "multiscale mapping ...\n");
    for (level = M; level >= 0; level--)
      { 
        /* Scale reduction factor from original strings: */
        int scale = 1 << level;

        fprintf(stderr, "\n");
        fprintf(stderr, "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n");
        fprintf(stderr, "LEVEL %d - SCALE %d\n", level, scale);

        int nx = nfx[level];
        int ny = nfy[level];
        int nw = nwt[level];
        
        /* Make descriptors for the two sequences at the current level: */
        msm_seq_desc_t x = msm_seq_desc_make(0, "x", level, nx, o->circX);
        msm_seq_desc_t y = msm_seq_desc_make(1, "y", level, ny, o->circY);

        /* Make filename prefix for current level: */
        char *levName = NULL;
        asprintf(&levName, "%s-%02d", o->outName, level);
        
        if (level == M)
          { /* Set {cdv} to an arbitrary set of initial candidates: */
            fprintf(stderr, "finding initial candidates at level %d ...\n", level);
            fprintf(stderr, " seq 0 length = %d\n", nx);
            fprintf(stderr, " seq 1 length = %d\n", ny);
            cdv = msm_fake_initial_cands(o->nCands, level, &x, &y);
          }
        else
          { /* Map the candidates in {cdv} to the next finer scale: */
            fprintf(stderr, "mapping from level %d to level %d ...\n", level+1, level);
            fprintf(stderr, "  sequence 0 length %d -> %d\n", nfx[level+1], nx);
            fprintf(stderr, "  sequence 1 length %d -> %d\n", nfy[level+1], ny);
            fprintf(stderr, "  filter window size %d\n", nw);
            msm_cand_vec_t cdm = msm_cand_vec_map_to_finer
              ( &cdv, nx, ny, x.circ, y.circ, nw );

            /* Reclaim previous cand list {cdv}: */
            /* msm_cand_vec_free(&cdv); */

            /* Write and plot the mapped candidates for this level: */
            fprintf(stderr, "writing mapped candidates ...\n");
            msm_cand_vec_write_named(&cdm, levName, "-m");
            msm_image_cand_vec_write_named(&cdm, &x, &y, levName, "-m");

            /* Interpolate the mapped candidates: */
            cdv = msm_cand_vec_interpolate(&cdm);

            /* Reclaim the mapped non-inerpolated candidates {cdm}: */
            /* msm_cand_vec_free(&cdm); */
          }
          
        /* Write and plot the official candidates for this level: */
        fprintf(stderr, "writing official candidates ...\n");
        msm_cand_vec_write_named(&cdv, levName, "-v");
        fprintf(stderr, "plotting the official candidates ...\n");
        msm_image_cand_vec_write_named(&cdv, &x, &y, levName, "-v");

        fprintf(stderr, "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n");
        fprintf(stderr, "\n");
      }

    return 0;
  }
    
void msm_compute_multiscale_lengths(int ns, bool_t circ, int maxLevel, int nwt[], int nfs[])
  { int level; 
    nfs[0] = ns;
    for (level = 1; level <= maxLevel; level++) 
      { int nFine = nfs[level-1];
        int nwtb = nwt[level-1];
        int nCoarse = msm_seq_desc_filtered_size(nFine, circ, nwtb);
        if (nCoarse < 0) { nCoarse = 0; }
        nfs[level] = nCoarse;
      }
  }

msm_cand_vec_t msm_fake_initial_cands
  ( int nCands,
    int level,
    msm_seq_desc_t *xp,
    msm_seq_desc_t *yp
  )
  { 
    /* Paranoia: */
    demand(xp->level == level, "inconsistent {level} of {xp}");
    demand(yp->level == level, "inconsistent {level} of {yp}");

    /* Define the max pairing length {maxlen}: */
    int nx = xp->npos;
    int ny = yp->npos;
    int maxlen = (nx > ny ? nx : ny);
    int minlen = (maxlen + 9)/10;
    /* Generate the candidates: */
    msm_cand_vec_t cdv = msm_cand_vec_throw
      ( nCands,
        xp, yp,
        minlen, maxlen,
        /*circProb*/ 0.1, 
        /*unitProb*/ 0.2,
        /*diagProb*/ 0.2,
        /*skipProb*/ 0.1
      );
    return cdv;
  }

int msm_seq_desc_filtered_size(int nFine, bool_t circ, int nwtb)
  { return (nFine - (circ ? 0 : nwtb - 1) + 1)/2; }

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
