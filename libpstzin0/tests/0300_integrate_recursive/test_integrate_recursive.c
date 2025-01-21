#define PROG_NAME "test_integrate_recursive"
#define PROG_DESC "test of {pst_integrate_iterative}"
#define PROG_VERS "1.0"

#define slope_to_height_C_COPYRIGHT "Copyright © 2005 by the State University of Campinas (UNICAMP)"

/* Last edited on 2025-01-18 12:30:18 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -maxIter {MAX_ITER} ] \\\n" \
  "    [ -convTol [CONV_TOL} ] \\\n" \
  "    [ -topoSort ] \\\n" \
  "    [ -compareZ {REF_Z_FNI_NAME} ] \\\n" \
  "    [ -verbose ] \\\n" \
  "    [ -debugG ] [ -debugZ ] [ -debugSys ] \\\n" \
  "    [ -debugIter {DBG_STEP} ] \\\n" \
  "    -slopes {IG_FNI_NAME} \\\n" \
  "    -weights {IW_FNI_NAME} \\\n" \
  "    -outPrefix {PREFIX} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "\n" \
  "OUTPUT FILES\n" \
  PROG_INFO_OUT_FILES "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  fni_to_pnm(1), pnm_to_fni(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2005-08-15 by Jorge Stolfi, UNICAMP.\n" \
  "MODIFICATION HISTORY\n" \
  "  All by J. Stolfi, UNICAMP, unless noted otherwise.\n" \
  "\n" \
  "  2008-11-08 added the weight map option.\n" \
  "  2008-11-09 added the reference height map.\n" \
  "  2010-04-26 options \"-newStyle\", \"-avgWidth\".\n" \
  "  2010-04-28\n" \
  "    * Computes a {Z} confidence map {OW} from the {G} confidence map {IW}.\n" \
  "    * \"-debugZ\"   implies output of Z and error maps at various levels and iterations.\n" \
  "    * \"-compareZ\" also writes one-line files with error summaries.\n" \
  "    * Explicit \"-slopes\" and \"-weights\" options.\n" \
  "    * Removed \"-debugPrefix\" added \"-outPrefix\".\n" \
  "    * The debug prefix is {PREFIX} plus \"-dbg\".\n" \
  "    * Output {OZ} goes to file, not to standard output.\n" \
  "  2010-05-04 Excludes zero-weight pixels from system.\n" \
  "  2025-01-05 Removed \"-newStyle\", \"-avgWidth\", \"-harmonicAvgW\".\n" \
  "  2025-01-05 Added \"-topoSort\".\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " slope_to_height_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define DEFAULT_MAX_ITER 100000
  /* Default max iterations per level. */
  
#define DEFAULT_CONV_TOL 0.0000005
  /* Default convergence threshold per level. */
  
#define DEFAULT_AVG_WIDTH 2
  /* Default kernel window size for map reduction. */
  
#define stringify(x) strngf(x)
#define strngf(x) #x

#define PROG_INFO_DESC \
  "  The program reads a /slope map/ {IG}, consisting of" \
  " the {X} and {Y} derivatives of some terrain surface {Z(X,Y)}; and" \
  " computes the corresponding /height map/ {OZ}.  Optionally, the" \
  " program may also compares the computed height map {OZ} with" \
  " a reference height map {RZ} provided by the user.\n" \
  "\n" \
  "  The input slope map is a two-channel float image {IG}.  Each" \
  " pixel {IG[X,Y]} of this image is" \
  " the gradient of the {Z} function, averaged" \
  " over the unit square cell (/pixel/) whose lower left corner" \
  " is the point {(X,Y)}.\n" \
  "\n" \
  "  The program may also read an optional /weight map/ {IW}, a" \
  " single-channel image that specifies the reliability" \
  " of the slope data.  Specifically, {IW[X,Y]} is the reliability" \
  " for the slope datum {IG[X,Y]}.  A zero weight means that the slope" \
  " datum should be ignored.\n" \
  "\n" \
  "  The output height map {OZ} is a single-channel" \
  " image containing the {Z} values of the terrain most" \
  " compatible with the given slopes.  It includes an arbitrary" \
  " global integration constant for each maximal part of the" \
  " domain which is connected by paths with non-zero weight.  This" \
  " image has one more column and one more row than the input slope maps.\n" \
  "\n" \
  "  If the \"-compareZ {REF_Z_NAME}\" option is given, the program" \
  " also outputs an error map {EZ}, defined as the difference between" \
  " the height map computed by the program and a reference map" \
  " read from the FNI-format file {REF_Z_NAME}.  The error map" \
  " is adjusted to have zero mean.  Any {Z} values that are" \
  " next to missing or unreliable slope data are excluded from the comparison.\n" \
  "\n" \
  "  If the \"-debugZ\" option is given, the program" \
  " also writes out the intermediate (reduced-scale)" \
  " solutions found, for debugging.  If \"-debugIter {DBG_STEP}\" is" \
  " given and {DBG_STEP} is nonzero, these diagnostics are written" \
  " every that many iterations.  The options" \
  " \"-debugG\" and \"-debugSys\" cause additional" \
  " debugging data to be written---namely, the reduced slope" \
  " maps and the linear systems solved at each scale," \
  " respectively."
  
#define PROG_INFO_OUT_FILES \
  "  {PREFIX}-Z.fni\n" \
  "    The computed height map {OZ}.\n" \
  "\n" \
  "  {PREFIX}-eZ.fni\n" \
  "    The error map {EZ} for the cmputed map {OZ}.\n" \
  "\n" \
  "  {PREFIX}-eZ.txt\n" \
  "    A one-line summary of the error, with the format\n" \
  "\n" \
  "     {level} {NX} {NY} {iter} {change} {sRZ} {sOZ} {sEZ} {sre}.\n" \
  "\n" \
  "    where {NX,NY} is the image size, {iter} is the number of" \
  " iterations performed at level 0, {change} is the change in" \
  " the last iteration, {sRZ} and {sOZ} are the deviations of" \
  " the {RZ} and {OZ} maps, {sEZ} is the RMS value of {EZ}, and" \
  " {sre} is the relative error.\n" \
  "\n" \
  "  The following files are produced for each level" \
  " of the recursion (except the topmost one in some" \
  " cases).  The {LEVEL} part is the recursion" \
  " level (0 = original scale), with 2 digits.  They are" \
  " requested by the \"-debugG\", \"-debugZ\", and \"-debugSys\" options," \
  " some of them in conjunction with \"-compareZ\".\n" \
  "\n" \
  "  {PREFIX}-dbg-{LEVEL}-dZ.fni\n" \
  "    The input slopes {IG}, reduced to level {LEVEL}.\n" \
  "\n" \
  "  {PREFIX}-dbg-{LEVEL}-W.fni\n" \
  "    The input slope weight mask {IW}, reduced to level {LEVEL}.\n" \
  "\n" \
  "  {PREFIX}-dbg-{LEVEL}-rZ.fni\n" \
  "    The reference height map {RZ}, reduced to level {LEVEL}.\n" \
  "\n" \
  "  {PREFIX}-dbg-{LEVEL}-Z.fni\n" \
  "    The computed height map {OZ} at level {LEVEL}.\n" \
  "\n" \
  "  {PREFIX}-dbg-{LEVEL}-U.fni\n" \
  "    The confidence mask {OW} for the heights at level {LEVEL}.\n" \
  "\n" \
  "  {PREFIX}-dbg-{LEVEL}-eZ.fni\n" \
  "    The height error map {EZ} for {OZ} and {RZ} at level {LEVEL}.\n" \
  "\n" \
  "  {PREFIX}-dbg-{LEVEL}-eZ.txt\n" \
  "    The error summary for level {LEVEL}.\n\n" \
  "\n" \
  "  {PREFIX}-dbg-{LEVEL}-S.sys\n" \
  "    The linear system solved at level {LEVEL}.\n" \
  "\n" \
  "  The following files may be produced at the beginning, during, and" \
  " after the end of the iteration.  The {ITER} part is" \
  " the number of iterations done, with 9 digits.  They are" \
  " requested by the \"-debugG\", \"-debugZ\", and \"-debugSys\" options," \
  " some of them with \"-compareZ\", together with a positive \"-debugIter\".\n" \
  "\n" \
  "  {PREFIX}-dbg-{LEVEL}-{ITER}-Z.fni\n" \
  "    The computed height map {OZ} at level {LEVEL} after {ITER} iterations.\n" \
  "\n" \
  "  {PREFIX}-dbg-{LEVEL}-{ITER}-eZ.fni\n" \
  "    The height error map {EZ} for {OZ} and {RZ} at level {LEVEL} after {ITER} iterations.\n" \
  "\n" \
  "  {PREFIX}-dbg-{LEVEL}-{ITER}-eZ.txt\n" \
  "    The height error summary for level {LEVEL} after {ITER} iterations."

#define PROG_INFO_OPTS \
  "  -slopes {IG_FNI_NAME}\n" \
  "    This mandatory argument specifies the file name" \
  " containing the input slope map {IG}, which must be a two-channel image in FNI format" \
  " (see float_image.h).  If {IG_FNI_NAME}" \
  " is \"-\", the slope map is read from standard input.\n" \
  "\n" \
  "  -weights {IW_FNI_NAME}\n" \
  "    This optional argument specifies the file containing the input" \
  " slope weight map {IW}.  It must be a single-channel image with non-negative" \
  " values, preferably between 0 and 1, and with the same size as {IG}.  Slope" \
  " samples with weight 0 are ignored.  If {IW_FNI_NAME}" \
  " is \"-\", the slope weight map is read from standard input.  If this argument is omitted, the" \
  " program assumes that all slope samples have weight 1.\n" \
  "\n" \
  "  -outPrefix {PREFIX}\n" \
  "    Specifies the common prefix for all output" \
  " file names.  Mandatory." \
  "\n" \
  "  -compareZ {REF_Z_FNI_NAME}\n" \
  "    If this argument is present, the program" \
  " also outputs an error map, defined as the difference between" \
  " the height map computed by the program and a reference map" \
  " read from the FNI-format file \"{REF_Z_NAME}\".  The reference" \
  " map must have either the same size as the input slope map, or one extra column and" \
  " one extra row.  In the first case, the computed height map" \
  " will be interpolated before the comparison.  The error map, with" \
  " the same size as the reference map, is written" \
  " to the file \"{ERR_Z_FNI_NAME}\".  If this option is" \
  " omitted, the reference map is not read and the error file" \
  " is not written.\n" \
  "\n" \
  "  -maxIter {MAX_ITER}\n" \
  "  -convTol {CONV_TOL}\n" \
  "    These optional parameters specify the stopping criterion for the" \
  " iteration.  The iteration will stop when the maximuum change in any" \
  " height field is less than {CONV_TOL}, or after {MAX_ITER} iterations, whichever" \
  " happens first. The defaults are {MAX_ITER = " stringify(DEFAULT_MAX_ITER) "}," \
  " {CONV_TOL = " stringify(DEFAULT_CONV_TOL) "}.\n" \
  "\n" \
  "  -topoSort\n" \
  "    This optional flag causes the program to solve the equations in order of" \
  " increasing equation weight.  If omitted, the equations are solved in row-by-row order.\n" \
  "\n" \
  "  -verbose\n" \
  "    Causes the program to write various diagnostics to stderr, in" \
  " particular one line per iteration of the Gauss-Seidel solver.\n" \
  "\n" \
  "  -debugG\n" \
  "    Causes the program to write the reduced versions" \
  " of the slope maps used at intermediate scales," \
  " in human-readable format.  The reduced weight maps are written too.\n" \
  "\n" \
  "  -debugZ\n" \
  "    Causes the program to write the solutions" \
  " (height maps) obtained at intermediate scales," \
  " in human-readable formats.  The height confidence maps {OW} are" \
  " also written.  If \"compareZ\" was" \
  " also given, the program also writes the error" \
  " map and error summaries for these reduced solutions.\n" \
  "\n" \
  "  -debugSys\n" \
  "    Causes the program to write the linear systems" \
  " used at intermediate scales.\n" \
  "\n" \
  "  -debugIter {DBG_STEP}\n" \
  "    If this optional argument is given and {DBG_STEP} is" \
  " not zero, the heights {OZ}, the error maps {EZ}, and the" \
  " error summaries are written out after every {DBG_STEP} iterations, as" \
  " well as at the beginning (0 iterations) and after the end of the" \
  " iterative solving.  If {DBG_STEP} is larger than {MAX_ITER} only" \
  " the initial and final states are written.  If {DBG_STEP} is zero (the default)," \
  " no \"...-{ITER}-...\" files are written."

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <values.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <jsstring.h>
#include <jsfile.h>
#include <argparser.h>
#include <float_image.h>
#include <float_image_mscale.h>
#include <float_image_expand_by_one.h>

#include <pst_slope_map.h>
#include <pst_imgsys.h>
#include <pst_height_map.h>

#include <pst_integrate_recursive.h>

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  {
    char *slopeFileName;  /* File name of input slope file, or "-" for stdin. */
    char *weightFileName; /* File name of weight slope file, or "-" for stdin (NULL if none). */
    char *refZFileName;   /* File name of input reference height map (NULL if none). */
    char *outPrefix;      /* Output file name prefix. */
    uint32_t maxIter;     /* Max iterations per level. */
    double convTol;       /* Convergence threshold. */
    bool_t topoSort;      /* TRUE solves the equations in order of increasing weight. */
    /* Debugging and analysis */
    bool_t verbose;       /* TRUE to print various diagnostics. */
    bool_t debugZ;        /* TRUE writes out the intermediate Z images. */
    bool_t debugG;        /* TRUE writes out the intermediate slope and weight arrays. */
    bool_t debugSys;      /* TRUE writes the linear system(s). */
    uint32_t debugIter;   /* Frequency for debugging output during iteration, or 0 if none. */
  } options_t;

float_image_t *read_fni_file(char *fileName);
  /* Reads a FNI image from file "{fileName}" (which should include the extension ".fni").
    If {fileName} is "-", reads from standard input. */

void write_fni_file(float_image_t *I, char *fileName, int32_t level);
  /* Writes the float image {I} in a human-readable format, to a file
    called "{fileName}" (which should include the extension ".fni").
    If {fileName} is "-", writes the image to standard output.

    If {level} is non-negative, diagnostic messages are indented by
    {2*level+2} spaces. */

options_t *parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

void compute_and_write_height_map(options_t *o, float_image_t *G, float_image_t *W, float_image_t *R);
  /* Computes the height map {OZ} from the gradient map {G} and its
    reliability map {W}. If {W} is null, assumes all weights are 1.
    If {R} is not null, compares {OZ} with {R} and writes the error
    map {EZ=OZ-R} and the error sumary file. Depending on the options
    {o} also writes these outputs at each level of the multiscale
    recursion, and possibly also after certain iterations at each
    level. */

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    fprintf(stderr, "Reading the slope map {G}:\n");
    float_image_t *G;  /* Input gradient map. */
    G = read_fni_file(o->slopeFileName);
    
    float_image_t *W; /* Input (gradient) reliability weight map. */
    if (o->weightFileName == NULL)
      { W = NULL; }
    else
      { fprintf(stderr, "Reading the slope reliability weight map {W}:\n");
        W = read_fni_file(o->weightFileName);
      }
    
    float_image_t *R; /* Reference height map. */
    if (o->refZFileName == NULL)
      { R = NULL; }
    else
      { fprintf(stderr, "Reading the reference height map {R}:\n");
        R = read_fni_file(o->refZFileName);
      }
    
    fprintf(stderr, "Computing the height map {OZ}:\n");
    compute_and_write_height_map(o, G, W, R);
    
    float_image_free(G); G = NULL;
    float_image_free(W); W = NULL;
    float_image_free(R); R = NULL;

    fprintf(stderr, "Done!\n");
    return 0;
  }
  
#define MAX_LEVEL 20
  /* Max level expected in recursion; that is {log_2} of max image width or height. */
  
void compute_and_write_height_map(options_t *o, float_image_t *G, float_image_t *R)
  {
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(G, &NC_G, &NX_G, &NY_G);
    if (o->verbose) { fprintf(stderr, "slope map size = %d × %d\n", NX_G, NY_G); }
    demand(NC_G == 3, "gradient map must have 3 channels");
    
    /* Height map size: */
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;
    if (o->verbose) { fprintf(stderr, "height map size = %d × %d\n", NX_Z, NY_Z); }
    
    /* Check the reference height map {R}: */
    float_image_t *E = NULL; /* Reference height map {R} expanded if needed. */
    if (R != NULL)
      { int32_t NC_R, NX_R, NY_R;
        float_image_get_size(R, &NC_R, &NX_R, &NY_R);
        if (o->verbose) { fprintf(stderr, "reference height map size = %d × %d\n", NX_R, NY_R); }
        if ((NX_R == NX_Z) && (NY_R == NY_Z))
          { E = R; }
        else if ((NX_R == NX_G) && (NY_R == NY_G))
          { E = float_image_expand_by_one(R, NULL); }
        else
          { demand(FALSE, "wrong ref Z size"); }
      }

    float_image_t *Z = float_image_new(1, NX_Z, NY_Z);
    float_image_t *U = float_image_new(1, NX_Z, NY_Z);

    char *debugPrefix = txtcat(o->outPrefix, "-deb"); /* Prefix for debug file names. */

    /* Data for each mscale level, saved by the reporting procs: */
    float_image_t *ms_E[MAX_LEVEL+1]; /* Version of {E} reduced for each level. */
    for (int32_t k = 0; k <= MAX_LEVEL; k++) { ms_E[k] = NULL; }

    auto void reportData(int32_t level, float_image_t *cur_G, float_image_t *cur_W, float_image_t *cur_Z);
    auto void reportSys(int32_t level, pst_imgsys_t *cur_S, float_image_t *cur_U);
    auto void reportHeights(int32_t level, int32_t iter, double change, bool_t final, float_image_t *cur_Z, float_image_t *cur_U);
      /* These procedures are called at various times during the recursive
         multiscale integration. See {pst_integrate_recursive}.
         The {reportHeights} procedure will write all the requested intermediate
         {Z} and {E} images and the error summaries, as well as 
         the final ones. */
    
    void reportData(int32_t level, float_image_t *cur_G, float_image_t *cur_Z)
      {
        /* This procedure is called once per level, on the way down.*/
        
        assert((level >= 0) && (level <= MAX_LEVEL));
        int32_t indent = (level < -1 ? 0 : 2*level+2);
        
        int32_t NX_cur_G, NY_cur_G;
        float_image_get_size(cur_G, NULL, &NX_cur_G, &NY_cur_G);
        assert(NC_cur_G == 3);
        
        if (level == 0) { assert(cur_Z == Z); }
        int32_t NX_cur_Z, NY_cur_Z;
        float_image_get_size(cur_Z, NULL, &NX_cur_Z, &NY_cur_Z);
        assert((NX_cur_Z == NX_cur_G + 1) && (NY_cur_Z == NY_cur_G + 1));
        
        /* Compute the reference height map {ms_E[level]} for this level: */
        if (E != NULL)
          { if (level == 0)
              { ms_E[level] = E; }
            else
              { if (o->verbose) { fprintf(stderr, "%*sshrinking reference height map ...\n", indent, ""); }
                uint32_t avgWidth = 2;
                ms_E[level] = pst_height_map_shrink(ms_E[level-1]);
              }
            int32_t NX_cur_E, NY_cur_E;
            float_image_get_size(ms_E[level], NULL, &NX_cur_E, &NY_cur_E);
            if (o->verbose) { fprintf(stderr, "%*sreference map size =  %d × %d\n", indent, "", NX_cur_E, NY_cur_E); }
            assert((NX_cur_E = NX_cur_Z) && (NY_cur_E == NY_cur_Z));
          }
        float_image_t *cur_E = ms_E[level];
        
        if (o->debugG)
          { /* Write out the scaled slope and weight maps: */
            float_image_mscale_write_file(cur_G, debugPrefix, level, -1, "dZ");
            if (cur_W != NULL)  { float_image_mscale_write_file(cur_W, debugPrefix, level, -1, "W"); }
            if (cur_E != NULL)  { float_image_mscale_write_file(cur_E, debugPrefix, level, -1, "R"); }
          }
      }
      
    void reportSys(int32_t level, pst_imgsys_t *cur_S, float_image_t *cur_U)
      {
        /* This procedure is called once per level on the way up. */
        assert((level >= 0) && (level <= MAX_LEVEL));
        int32_t indent = (level < -1 ? 0 : 2*level+2);
        
        if (cur_S == NULL)
          { fprintf(stderr, "%*sno equation system was generated at this level\n", indent, ""); }
        else
          { if (o->debugSys) { pst_imgsys_write_report(cur_S, debugPrefix, level, "S"); }
            pst_imgsys_extract_system_eq_tot_weight_image(cur_S, cur_U, 0.0);
            float_image_mscale_write_file(cur_U, debugPrefix, level, -1, "U");
          }
      }
    
    void reportHeights(int32_t level, int32_t iter, double change, bool_t final, float_image_t *cur_Z, float_image_t *cur_U)
      { 
        /* This procedure is called at least once per level on the way up, including once with {final=TRUE}. */
        assert((level >= 0) && (level <= MAX_LEVEL));
        assert(cur_Z != NULL);
        assert(cur_U != NULL);

        /* Decide whether to write the data, and how: */
        bool_t writeAsIter = (o->debugIter > 0) && (final || (iter % (int32_t)o->debugIter == 0));
        bool_t writeAsLevel = o->debugZ && final;
        bool_t writeAsFinal = (level == 0) && final;
        
        bool_t analyzeError = (R != NULL);
        
        if (writeAsIter)
          { bool_t writeImages = o->debugZ;
            bool_t writeError = writeImages && analyzeError;
            if (writeError) { assert(ms_E[level] != NULL); }
            pst_height_map_level_analyze_and_write
              ( debugPrefix, level, iter, change, 
                cur_Z, ms_E[level], cur_U, 
                writeImages, writeError
              );
          }
        
        if (writeAsLevel)
          { bool_t writeImages = o->debugZ;
            bool_t writeError = writeImages && analyzeError;
            if (writeError) {  assert(ms_E[level] != NULL); }
            pst_height_map_level_analyze_and_write
              ( debugPrefix, level, -1, change, 
                cur_Z, ms_E[level], cur_U, 
                writeImages, writeError
              );
          }
        
        if (writeAsFinal)
          { bool_t writeImages = TRUE;
            bool_t writeError = writeImages && analyzeError;
            pst_height_map_level_analyze_and_write
              ( o->outPrefix, -1, -1, change, 
                cur_Z, ms_E[level], cur_U, 
                writeImages, writeError
              );
          }
      }

    /* Call recursive integrator: */
    
    bool_t keepNull = FALSE; /* ??? Should it be a parameter ??? */
    pst_integrate_recursive
      ( G, W, 
        keepNull, 
        Z, U,
        0,
        o->maxIter, o->convTol, o->topoSort,
        o->verbose,
        &reportData, 
        &reportSys, 
        o->debugIter,
        &reportHeights
      );
      
    /* Free working storage: */
    if ((E != NULL) && (E != R)) { float_image_free(E); }
    float_image_free(Z);
    float_image_free(U);
    free(debugPrefix);
    for (uint32_t k = 1; k <= MAX_LEVEL; k++) 
      { if (ms_E[k] != NULL) { float_image_free(ms_E[k]); }
      }
  }

float_image_t *read_fni_file(char *fileName)
  { demand(fileName != NULL, "file name not given");
    fprintf(stderr, "Reading %s ...", fileName);
    FILE *rd = open_read(fileName, FALSE);
    float_image_t *I = float_image_read(rd);
    if (rd != stdin) { fclose(rd); }
    fprintf(stderr, "\n");
    return I;
  }

void write_fni_file(float_image_t *I, char *fileName, int32_t indent)
  { demand(fileName != NULL, "file name not given");
    fprintf(stderr, "%*sWriting %s ...", indent, "", fileName);
    FILE* wr = open_write(fileName, FALSE);
    float_image_write(wr, I);
    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
    fprintf(stderr, "\n");
  }

options_t *parse_options(int32_t argc, char **argv)
  { argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 

    if (argparser_keyword_present(pp, "-maxIter"))
      { o->maxIter = (uint32_t)argparser_get_next_int(pp, 0, INT64_MAX); }
    else
      { o->maxIter = DEFAULT_MAX_ITER; }
    
    if (argparser_keyword_present(pp, "-convTol"))
      { o->convTol = argparser_get_next_double(pp, 0, DBL_MAX); }
    else
      { o->convTol = DEFAULT_CONV_TOL; }
    
    o->topoSort = argparser_keyword_present(pp, "-topoSort");
    
    if (argparser_keyword_present(pp, "-compareZ"))
      { o->refZFileName = argparser_get_next(pp); }
    else
      { o->refZFileName = NULL; }
    
    o->verbose = argparser_keyword_present(pp, "-verbose");
    
    o->debugZ = argparser_keyword_present(pp, "-debugZ");
    o->debugG = argparser_keyword_present(pp, "-debugG");
    o->debugSys = argparser_keyword_present(pp, "-debugSys");
    
    if (argparser_keyword_present(pp, "-debugIter"))
      { o->debugIter = (uint32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->debugIter = 0; }

    argparser_get_keyword(pp, "-slopes");
    o->slopeFileName = argparser_get_next_non_keyword(pp);
    
    if (argparser_keyword_present(pp, "-weights"))
      { o->weightFileName = argparser_get_next(pp); }
    else
      { o->weightFileName = NULL; }
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    
    argparser_finish(pp);
    
    return o;
  }
