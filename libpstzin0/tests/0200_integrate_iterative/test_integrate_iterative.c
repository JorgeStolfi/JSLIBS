#define PROG_NAME "test_integrate_iterative"
#define PROG_DESC "test of {pst_integrate_iterative}"
#define PROG_VERS "1.0"

#define slope_to_height_C_COPYRIGHT "Copyright © 2005 by the State University of Campinas (UNICAMP)"

/* Last edited on 2025-01-18 12:32:29 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -slopes {IG_FNI_NAME} \\\n" \
  "    -weights {IW_FNI_NAME} \\\n" \
  "    [ -keepNull {KEEPNULL} ] \\\n" \
  "    [ -initial {IZ_FNI_NAME} ] \\\n" \
  "    [ -maxIter {MAX_ITER} ] \\\n" \
  "    [ -convTol [CONV_TOL} ] \\\n" \
  "    [ -topoSort [TOPOSORT} ] \\\n" \
  "    [ -verbose ] \\\n" \
  "    [ -reportStep {REPORT_STEP} ] \\\n" \
  "    [ -compareZ {REF_Z_FNI_NAME} ] \\\n" \
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
  "  2008-11-08 by J. Stolfi, IC-UNICAMP: added the weight map option.\n" \
  "\n" \
  "  2008-11-09 by J. Stolfi, IC-UNICAMP: added the reference height map.\n" \
  "\n" \
  "  2010-04-28 by J. Stolfi, IC-UNICAMP:\n" \
  "    * Computes a {Z} confidence map {OW} from the {G} confidence map {IW}.\n" \
  "    * \"-compareZ\" also writes one-line files with error summaries.\n" \
  "    * Explicit \"-slopes\" and \"-weights\" options.\n" \
  "    * Removed \"-debugPrefix\" added \"-outPrefix\".\n" \
  "    * The debug prefix is {PREFIX} plus \"-dbg\".\n" \
  "    * Output {OZ} goes to file, not to standard output.\n" \
  "\n" \
  "  2010-05-04 by J. Stolfi, IC-UNICAMP:" \
  "    * Added \"-topoSort\".\n" \
  "\n" \
  "  2025-01-07 by J. Stolfi, IC-UNICAMP:" \
  "    * Option \"-keepNull\" to keep/excludes zero-weight equations.\n" \
  "    * Option \"-initial\" to give optional initial guess.\n" \
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
  "  If the \"-reportStep\" option is given, the program" \
  " also writes out the current height map at each" \
  " iteration.  If \"-reportStep {REPORT_STEP}\" is" \
  " given and {REPORT_STEP} is nonzero, these diagnostics are written" \
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
  "     {NX} {NY} {iter} {change} {sRZ} {sOZ} {sEZ} {sre}.\n" \
  "\n" \
  "    where {NX,NY} is the image size, {iter} is the number of" \
  " iterations performed, {change} is the change in" \
  " the last iteration, {sRZ} and {sOZ} are the deviations of" \
  " the {RZ} and {OZ} maps, {sEZ} is the RMS value of {EZ}, and" \
  " {sre} is the relative error.\n" \
  "\n" \
  "  The following files may be produced at the beginning, during, and" \
  " after the end of the iteration, if the \"-reportStep\" option" \
  " is given with non-zero alrgument.  The {ITER} part is" \
  " the number of iterations done, with 9 digits.  " \
  " some of them with \"-compareZ\", together with a positive \"-reportStep\".\n" \
  "\n" \
  "  {PREFIX}-dbg-{ITER}-Z.fni\n" \
  "    The computed height map {OZ}  after {ITER} iterations.\n" \
  "\n" \
  "  {PREFIX}-dbg-{ITER}-eZ.fni\n" \
  "    The height error map {EZ} for {OZ} and {RZ} after {ITER} iterations.\n" \
  "\n" \
  "  {PREFIX}-dbg-{ITER}-eZ.txt\n" \
  "    The height error summary after {ITER} iterations."

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
  "  -topoSort {TOPOSORT}\n" \
  "    This optional boolean argument specifies the order in which" \
  " equations are solved by the program.  The {TOPOSORT} may be 'T' or 1 to" \
  " solve in order of increasing equation weight, or 'F' of 0 to solve in" \
  " scan line orrder.  If omitted, the program assumes \"-topoSort F\".\n" \
  "\n" \
  "  -keepNull {KEEPNULL}\n" \
  "    This optional boolean argument specifies whether equations with total weight" \
  " zero (isolated height pixels) are to be" \
  " retained (if {KEEPNULL} is 'T' or 1) or excluded from the" \
  " linear system (if {KEEPNULL} is 'F' or 0).  If omitted, the" \
  " program assumes \"-keepNull F\".\n" \
  "\n" \
  "  -verbose\n" \
  "    Causes the program to write various diagnostics to stderr, in" \
  " particular one line per iteration of the Gauss-Seidel solver.\n" \
  "\n" \
  "  -reportStep {REPORT_STEP}\n" \
  "    If this optional argument is given and {REPORT_STEP} is" \
  " not zero, the heights {OZ}, the error maps {EZ}, and the" \
  " error summaries are written out after every {REPORT_STEP} iterations, as" \
  " well as at the beginning (0 iterations) and after the end of the" \
  " iterative solving.  If {REPORT_STEP} is larger than {MAX_ITER}, only" \
  " the initial and final states are written.  If" \
  " {REPORT_STEP} is zero (the default)," \
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
#include <pst_integrate.h>
#include <pst_integrate_iterative.h>

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  {
    char *slopes_fname;      /* File name of input slope file, or "-" for stdin. */
    char *weights_fname;     /* File name of weight slope file, or "-" for stdin (NULL if none). */
    char *initial_fname;     /* File name of weight slope file, or "-" for stdin (NULL if none). */
    char *compareZ_fname;    /* File name of input reference height map (NULL if none). */
    bool_t keepNull;         /* True to include equations for isolated height pixels. */
    uint32_t maxIter;        /* Max iterations. */
    double convTol;          /* Convergence threshold. */
    bool_t topoSort;         /* TRUE solves the equations in order of increasing weight. */
    char *outPrefix;         /* Output file name prefix. */
    /* Debugging and analysis */
    bool_t verbose;          /* TRUE to print various diagnostics. */
    uint32_t reportStep;      /* Frequency for debugging output during iteration, or 0 if none. */
  } options_t;

float_image_t *read_fni_file(char *fileName);
  /* Reads a FNI image from file "{fileName}" (which should include the extension ".fni").
    If {fileName} is "-", reads from standard input. */

void write_fni_file(float_image_t *I, char *fileName, int32_t indent);
  /* Writes the float image {I} in a human-readable format, to a file
    called "{fileName}" (which should include the extension ".fni").
    If {fileName} is "-", writes the image to standard output. 
    Diagnostic images are indented by {indent} spaces. */
    
void write_system_weights(pst_imgsys_t *S, char *fname);
  /* Writes the equation weights {S.eq[k].wtot} as a PNG image "{fname}". */

options_t *parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

void compute_and_write_height_map
  ( options_t *o,
    float_image_t *G,
    float_image_t *W,
    float_image_t *Z,
    float_image_t *U,
    float_image_t *RZ
  );
  /* Computes the height map {Z} from the gradient map {G} and its
    reliability map {W}.
    
    The map {G} must have two channels (intrepreted as {dZ/dX} and {dZ/dY}),
    and the map {W}, if not null, must have a single channel and the same
    size as {G}.  Sample {W[0,x,y]}, which must be non-negative,
    is interpreted as the reliability of {G[*,x,y]}; if zero, the value
    of {G[*x,y]} will be ignored.  If {W} is null, assumes all weights are 1.
    
    The image {Z} must have one channel, and one row and one col more than {G}
    and {W}.  On input, it gives the initial guess for the iterative method.
    On output, it will contain the computed height map. 
    
    If {o.reportStep} is not zero, writes the height map
    for the first {o.reportStep} iterations, and then 
    every {o.reportStep} iterations thereafter.
    
    If {RZ} is not null, it must be a single-channel image with the same
    size as {G} or {Z}. Whenever {Z} is reported, the procedure compares
    {Z} with {RZ} and writes the error map {EZ=Z-RZ} and the error
    sumary file. */

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    fprintf(stderr, "Reading the slope map {G}:\n");
    float_image_t *G;  /* Input gradient map. */
    G = read_fni_file(o->slopes_fname);
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(G, &NC_G, &NX_G, &NY_G);
    demand(NC_G == 3, "gradient map must have 3 channels");
    
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;
    float_image_t *Z; /* Input initial guess for height map. */
    if (o->initial_fname == NULL)
      { Z = float_image_new(1, NX_Z, NY_Z);
        float_image_fill_channel(Z, 0, 0.0);
      }
    else
      { fprintf(stderr, "Reading the intial height map {Z}:\n");
        Z = read_fni_file(o->initial_fname);
        float_image_check_size(Z, 2, NX_Z, NY_Z, "bad height map");
      }
    
    float_image_t *RZ; /* Reference height map. */
    if (o->compareZ_fname == NULL)
      { RZ = NULL; }
    else
      { fprintf(stderr, "Reading the reference height map {RZ}:\n");
        float_image_t *IRZ = read_fni_file(o->compareZ_fname);
        int32_t NC_R, NX_R, NY_R;
        float_image_get_size(IRZ, &NC_R, &NX_R, &NY_R);
        demand(NC_R == 1, "reference height map must have 1 channel");
        if ((NX_R == NX_Z) && (NY_R == NY_Z))
          { RZ = IRZ; }
        else if ((NX_R == NX_G) && (NY_R == NY_G))
          { RZ = float_image_expand_by_one(IRZ, NULL); float_image_free(IRZ); }
        else
          { demand(FALSE, "wrong ref height map {IRZ} size"); }
      }
    
    fprintf(stderr, "Computing the height map {Z}:\n");
    compute_and_write_height_map(o, G, W, Z, U, RZ);
    
    float_image_free(G); G = NULL;
    if (W != NULL) { float_image_free(W); W = NULL; }
    if (Z != NULL) { float_image_free(Z); Z = NULL; }
    if (RZ != NULL) { float_image_free(RZ); RZ = NULL; }

    fprintf(stderr, "Done!\n");
    return 0;
  }
  
void compute_and_write_height_map
  ( options_t *o,
    float_image_t *G,
    float_image_t *Z,
    float_image_t *RZ
  )
  {
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(G, &NC_G, &NX_G, &NY_G);
    demand(NC_G == 3, "gradient map must have 3 channels");
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;
    float_image_check_size(Z, 2, NX_Z, NY_Z, "bad height map");
    
    /* Check the reference height map {RZ}, expand it of needed: */
    float_image_t *EZ = NULL; /* Version of {RZ} to use when comparing with {Z}. */
    if (RZ != NULL)
      { int32_t NC_R, NX_R, NY_R;
        float_image_get_size(RZ, &NC_R, &NX_R, &NY_R);
        if ((NX_R == NX_Z) && (NY_R == NY_Z))
          { EZ = RZ; }
        else if ((NX_R == NX_G) && (NY_R == NY_G))
          { EZ = float_image_expand_by_one(RZ, NULL); }
        else
          { demand(FALSE, "wrong ref height map {RZ} size"); }
      }

    char *debugPrefix = txtcat(o->outPrefix, "-dbg"); /* Prefix for debug file names. */

    auto void reportSys(int32_t level, pst_imgsys_t *S, float_image_t *U);
      /* This procedure is called by {pst_integrate_iterative} once 
        to report the linear system created from the slope and weight maps. */

    auto void reportHeights(int32_t level, int32_t iter, double change, bool_t final, float_image_t *Z, float_image_t *U);
      /* This procedure is called by {pst_integrate_iterative} to write
         the current height map {Z} and the height error map {EZ}, and
         the error summary. */

    /* Call iterative integrator: */
    int32_t level = -1;
    pst_integrate_iterative
      ( G, W, o->keepNull, Z, U,
        o->maxIter, o->convTol, o->topoSort,
        o->verbose,
        level,
        &reportSys,
        o->reportStep,
        &reportHeights
      );
      
    /* Free working storage: */
    free(debugPrefix);
    if ((EZ != NULL) && (EZ != RZ)) { float_image_free(EZ); }
      
    return;
    
    void reportHeights(int32_t level, int32_t iter, double change, bool_t final, float_image_t *Z, float_image_t *U)
      { 
        assert(level == -1);
         
        bool_t writeImages = TRUE;
        bool_t writeError = (EZ != NULL);
        pst_height_map_level_analyze_and_write
          ( debugPrefix, level, (final ? -1 : iter), change, 
            Z, EZ, U, 
            writeImages, writeError
          );
      }
      
    void reportSys(int32_t level, pst_imgsys_t *S, float_image_t *U)
      { char *fileName_sys = float_image_mscale_file_name(debugPrefix, level, -1, "sys", "txt");
        FILE* wr_sys = fopen(fileName_sys, "wt");
        pst_imgsys_write(wr_sys, S);
        fclose(wr_sys);
        free(fileName_sys);
        
        char *fileName_U = float_image_mscale_file_name(debugPrefix, level, -1, "U", "fni");
        float_image_write_named(fileName_U, U);
        free(fileName_U);
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

    argparser_get_keyword(pp, "-slopes");
    o->slopes_fname = argparser_get_next_non_keyword(pp);
    
    if (argparser_keyword_present(pp, "-weights"))
      { o->weights_fname = argparser_get_next_non_keyword(pp); }
    else
      { o->weights_fname = NULL; }
    
    if (argparser_keyword_present(pp, "-initial"))
      { o->initial_fname = argparser_get_next_non_keyword(pp); }
    else
      { o->initial_fname = NULL; }

    if (argparser_keyword_present(pp, "-keepNull"))
      { o->keepNull = argparser_get_next_bool(pp); }
    else
      { o->keepNull = FALSE; }

    if (argparser_keyword_present(pp, "-maxIter"))
      { o->maxIter = (uint32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->maxIter = DEFAULT_MAX_ITER; }
    
    if (argparser_keyword_present(pp, "-convTol"))
      { o->convTol = argparser_get_next_double(pp, 0, DBL_MAX); }
    else
      { o->convTol = DEFAULT_CONV_TOL; }

    if (argparser_keyword_present(pp, "-topoSort"))
      { o->topoSort = argparser_get_next_bool(pp); }
    else
      { o->topoSort = FALSE; }
    
    if (argparser_keyword_present(pp, "-compareZ"))
      { o->compareZ_fname = argparser_get_next(pp); }
    else
      { o->compareZ_fname = NULL; }
    
    o->verbose = argparser_keyword_present(pp, "-verbose");
    
    if (argparser_keyword_present(pp, "-reportStep"))
      { o->reportStep = (uint32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->reportStep = 0; }
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    
    argparser_finish(pp);
    
    return o;
  }
