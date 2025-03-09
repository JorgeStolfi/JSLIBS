#define PROG_NAME "test_integrate_recursive"
#define PROG_DESC "test of {pst_integrate_recursive}"
#define PROG_VERS "1.0"

#define test_integrate_recursive_C_COPYRIGHT "Copyright © 2005 by the State University of Campinas (UNICAMP)"

/* Last edited on 2025-03-06 13:28:39 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    { -slopes {G_FNI_NAME} | -normals {N_FNI_NAME} } [ scale {G_SX} {G_SY} ] \\\n" \
  "    [ -reference {R_FNI_NAME} [ scale {R_SZ} ] \\\n" \
  "    [ -hints {H_FNI_NAME} [ scale {H_SZ} ] weight {H_WT} ] \\\n" \
  "    [ -initial {INIT_OPT} {INIT_NOISE} ] \\\n" \
  "    [ -clear {CLX_MIN} {CLX_MAX} {CLY_MIN} {CLY_MAX} ... ] \\\n" \
  "    [ -maxLevel {MAX_LEVEL} ] \\\n" \
  "    [ -convTol [CONV_TOL} ] \\\n" \
  "    [ -maxIter {MAX_ITER} ] \\\n" \
  "    [ -sortSys [TOPO_SORT} ] \\\n" \
  "    [ -verbose ] \\\n" \
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
  "  By J. Stolfi unless said otherwise.\n" \
  "\n" \
  "  2008-11-08 added the weight map option.\n" \
  "  2008-11-09 added the reference height map.\n" \
  "  2010-04-28 \"-reference\" also writes one-line files with error summaries.\n" \
  "  2010-04-28 Explicit \"-slopes\" option.\n" \
  "  2010-04-28 Removed \"-debugPrefix\" added \"-outPrefix\".\n" \
  "  2010-04-28 The debug prefix is {PREFIX} plus \"-dbg\".\n" \
  "  2010-04-28 Output {Z} goes to file, not to standard output.\n" \
  "\n" \
  "  2010-05-04 Added \"-sortSys\".\n" \
  "\n" \
  "  2010-05-04 Excludes zero-weight pixels from system.\n" \
  "  2025-01-05 Added \"-sortSys\".\n" \
  "  2025-01-07 Option \"-initial\" to give optional initial guess.\n" \
  "  2025-02-13 Added \"-hints\" option.\n" \
  "  2025-02-18 Made  \"-reference\" optional.\n" \
  "  2025-02-22 General cleanup and merge of {test_integrate_iterative}.\n" \
  "  2025-02-23 Added \"-maxLevel\".\n" \
  "  2025-02-18 Added \"-normals\" as alternative to \"-slopes\".\n" \
  "  2025-03-02 Created subfolder \"{outPrefix}-iters/\" for iteration reports.\n" \
  "  2025-03-03 Fixed file names in {PROG_INFO} and code.\n" \
  "  2025-03-04 Added the \"scale\" options for the {G,H,R} maps.\n" \
  "  2025-03-04 Added the \"-clear\" option.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_integrate_recursive_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define DEFAULT_MAX_LEVEL 30
  /* Anything greater than the {log_2} of the max image dimension. */

#define DEFAULT_MAX_ITER 100000
  /* Default max iterations per level. */

#define DEFAULT_CONV_TOL 0.0000005
  /* Default convergence threshold per level. */

#define stringify(x) strngf(x)
#define strngf(x) #x

#define PROG_INFO_DESC \
  "  The program reads a /slope map/ {G}, consisting of" \
  " the {X} and {Y} derivatives of some terrain surface {Z(X,Y)}; and" \
  " computes the corresponding /height map/ {Z} from that data.  Optionally, the" \
  " program reads a one- or two-channel /hints map/ {H} with an independent" \
  " estimate of the desired height map {Z}.  Optionally, the" \
  " program may also compare the computed height map {Z} with" \
  " a /reference height map/ {R} provided by the user.\n" \
  "\n" \
  "  Ideally the maps {H} and {R} must have either the same size" \
  " as the height map to be computed, that is, one col and one row more than the" \
  " slope (or normals) map.  However, for convenience, the program also" \
  " accepts, in each case, a map with the same size as the slope map. In this case, this" \
  " map will be expanded by one col and one row with {float_image_expand_by_one}.\n" \
  "\n" \
  "  Instead of a slope map, the program may read a /normal map/ {N} that" \
  " specifies the unit normal vector in each cell of the grid.  The" \
  " program converts {N} internally to a slope map {G}, and proceeds" \
  " as in that case.\n" \
  "\n" \
  "  " pst_integrate_recursive_INFO("H_WT","INITIAL_OPT","MAX_ITER","CONV_TOL","TOPO_SORT") "\n" \
  "\n" \
  "  If not specified by the user with the \"-maxLevel\" command" \
  " line argument, the recursion will stop the slope map is" \
  " reduced to a single pixel.\n" \
  "\n" \
  "  At some selected iterations, the program" \
  " will compute an error map {E}, defined as the difference between" \
  " the height map {Z} computed by the program and a reference height" \
  " map {R} specified by the \"-reference\" command line" \
  " argument.  At selected iterations, the error map {E}  is written" \
  " to the file \"{ERR_Z_FNI_NAME}\".  The" \
  " error map {E} is shifted to have zero mean.  Any {Z} values that are" \
  " not finite or have zero weight in {R}.\n" \
  "\n" \
  "  If the \"-reportStep\" option is given, the program" \
  " also writes out the current height map at each" \
  " iteration.  If \"-reportStep {REPORT_STEP}\" is" \
  " given and {REPORT_STEP} is nonzero, these diagnostics are written" \
  " every that many iterations.\n" \
  "\n" \
  "  The output height map {Z} is a single-channel" \
  " image containing the {Z} values of the terrain most" \
  " compatible with the given slopes.  It includes an arbitrary" \
  " global integration constant for each maximal part of the" \
  " domain which is connected by paths with non-zero weight.  This" \
  " image has one more column and one more row than the input" \
  " slope (or normal) map."
  
#define PROG_INFO_OUT_FILES \
  "  The following files will be produced for each level" \
  " of the recursion  The {LEVEL} is the recursion level formatted as '%02d'  The" \
  " input maps and the final solution are those with {LEVEL=\"00\"}:\n" \
  "\n" \
  "    {PREFIX}-{LEVEL}-beg-G.fni\n" \
  "      The input slopes {G}, reduced to level {LEVEL}.\n" \
  "\n" \
  "    {PREFIX}-{LEVEL}-beg-R.fni\n" \
  "      The reference height map {R}, reduced to level {LEVEL}.\n" \
  "\n" \
  "    {PREFIX}-{LEVEL}-beg-H.fni\n" \
  "      The hints height map {H}, reduced to level {LEVEL}.\n" \
  "\n" \
  "    {PREFIX}-{LEVEL}-beg-S.sys\n" \
  "      The linear system that was solved at level {LEVEL}.\n" \
  "\n" \
  "    {PREFIX}-{LEVEL}-beg-Z.fni\n" \
  "    {PREFIX}-{LEVEL}-end-Z.fni\n" \
  "      The initial and computed height map {Z} at level {LEVEL}.\n" \
  "\n" \
  "    {PREFIX}-{LEVEL}-beg-E.fni\n" \
  "    {PREFIX}-{LEVEL}-end-E.fni\n" \
  "      The initial and final height error map {E=Z-R} at" \
  " level {LEVEL}, shifted to zero mean.\n" \
  "\n" \
  "    {PREFIX}-{LEVEL}-beg-E.txt\n" \
  "    {PREFIX}-{LEVEL}-end-E.txt\n" \
  "      A summary of the initial and final error at" \
  " level {LEVEL} (see below).\n" \
  "\n" \
  " The following files too may be produced during the iterative solving at" \
  " each level, if the \"-reportStep\" option is given with non-zero alrgument.  The" \
  " {ITER} part is the number of iterations done, with 9 digits.   These files do not" \
  " include the initial state ({ITER} zero) and the final state, which are given above.\n" \
  "\n" \
  "  {PREFIX}-iters/it-{LEVEL}-{ITER}-Z.fni\n" \
  "    The computed height map {Z} at level {LEVEL} after {ITER} iterations.\n" \
  "\n" \
  "  {PREFIX}-iters/it-{LEVEL}-{ITER}-E.fni\n" \
  "    The height error map {E} for {Z} and {R} at level {LEVEL} after {ITER} iterations.\n" \
  "\n" \
  "  {PREFIX}-iters/it-{LEVEL}-{ITER}-E.txt\n" \
  "    The height error summary for level {LEVEL} after {ITER} iterations.\n" \
  "\n" \
  "ERROR SUMMARY FILES\n" \
   "  " pst_map_compare_error_summary_INFO("Z","R") ""

#define PROG_INFO_OPTS \
  "  -slopes {G_FNI_NAME} [ scale {G_SX} {G_SY} ]\n" \
  "  -normals {N_FNI_NAME} [ scale {G_SX} {G_SY} ]\n" \
  "    Exactly one of these two arguments must be given, to" \
  " specify the file name containing the input slope map {G}, or" \
  " the input normal map {N}.  In either case the file must" \
  " in the FNI format (see {float_image_read}).\n" \
  "\n" \
  "    The slope map must be a two- or three-channel image, containing" \
  " the derivatives {dZ/dX} and {dZ/dY} averaged in each cell of" \
  " the image domain grid, and an optional reliability" \
  " weight for that data (see pst_slope_map.h).\n" \
  "\n" \
  "    The normal map must be a three- or four-channel, containing" \
  " the coordinates {Nx,Ny,Nz} of the unint vector which is the mean" \
  " surface normal in each cell, and an optional reliability" \
  " weight for the same (see {pst_normal_map.h})." \
  "\n" \
  "    If the \"scale\" keyword follows the map name, the" \
  " coordinates {dZ/dX} and {dZ/dY} of each gradient (read from the" \
  " file or computed from the normal map) are multiplied by {G_SX}" \
  " and {G_SY}, respectively.  In particular, the {G_SY} parameter" \
  " should be negative if the direction of the {Y} axis assumed" \
  " in the map is opposite to that assumed by this program." \
  "\n" \
  "  -hints {H_FNI_NAME} [ scale {H_SZ} ] {H_WT}\n" \
  "    This optional argument specifies the file name" \
  " containing the independent estimated height map {H}, which must" \
  " be a one- or two-channel image in FNI format" \
  " (see pst_height_map.h).  The argument {H_WT} must be a number" \
  " between 0 and 1 that will be multiplied  by the reliability weight" \
  " of every height in {H}.  If the \"scale\" keyword is" \
  " present, the {Z} values in the map {H} will be multiplied" \
  " by {H_SZ}.\n" \
  "\n" \
  "  -reference {R_FNI_NAME} [ scale {R_SZ} ]\n" \
  "    This optional argument specifies the file that contains the" \
  " expecied height map.  If this the refernce map is not" \
  " specified, no error maps or error summaries will" \
  " be computed.  If the \"scale\" keyword is" \
  " present, the {Z} values in the map {R} will be multiplied" \
  " by {R_SZ}.\n" \
  "\n" \
  "  -initial {INIT_OPT} {INIT_NOISE}\n" \
  "    This optional argument specifies the initial guess of the height" \
  " map for the iterative solver, at the deepest level" \
  " (smallest map).  The {INIT_OPT} may be either \"zero\", to" \
  " mean all initial heights are zero, \"hints\" to mean they are to be" \
  " copied from the hints map {H}, or \"reference\", to mean that" \
  " they are to be copied from the reference height map {R}, both reduced" \
  " in size to the deepest level and scaled by the" \
  " respective \"scale\" factors.  In any case" \
  " case, if {INIT_NOISE} is nonzero, the intial map specified as" \
  " above will be perturbed before" \
  " the iteration starts.  The perturbation consists in adding to every" \
  " finite height value in the map {Z} a pseudorandom amount in the" \
  " interval {[-mag _ +mag]}, where {mag} is {INIT_NOISE} times" \
  " {hypot(rad, vmax-vmin)}, {rad} is half of the smallest" \
  " dimension (cols or rows) of the map {Z}, and {[vmin _ vmax]} the" \
  " original range of the heights in {Z} (ignoring non-finite values).  If this" \
  " argument is not specified, the default is \"-initial zero 0.0\".\n" \
  "\n" \
  "  -clear {CLX_MIN} {CLX_MAX} {CLY_MIN} {CLY_MAX}}\n" \
  "    Each occurrence of this optional argument specifies" \
  " a rectangle of the {G}, {H}, or {R} maps" \
  " where all samples should be replaced by {NAN}s with" \
  " weight zero before starting.  For the {G} map, that" \
  " applies to pixels {G[c,X,Y]} with {X} in" \
  " {CLX_MIN..CLX_MAX} and {Y} in {CLY_MIN..CLY_MAX}.  For" \
  " the {H} and {R} maps, the upper limits of the ranges" \
  " are {CLX_MAX+1} and {CLY_MAX+1}.\n" \
  "\n" \
  "  -maxLevel {MAX_LEVEL}\n" \
  "    This optional argument specifies the maximum recursion level.  If" \
  " not specified, or is too large, the recursion will continue until" \
  " the slope map {G} is reduced to a single pixel.  If {MAXLEVEL} is" \
  " zero, recursion will be suppressed, and the heights will" \
  " be computed directly from the given slope map.\n" \
  "\n" \
  "  -maxIter {MAX_ITER}\n" \
  "  -convTol {CONV_TOL}\n" \
  "    These optional parameters specify the stopping criterion for the" \
  " iteration.  The iteration will stop when the maximuum change in any" \
  " height field is less than {CONV_TOL}, or after {MAX_ITER} iterations, whichever" \
  " happens first. The defaults are {MAX_ITER = " stringify(DEFAULT_MAX_ITER) "}," \
  " {CONV_TOL = " stringify(DEFAULT_CONV_TOL) "}.\n" \
  "\n" \
  "  -sortSys {TOPO_SORT}\n" \
  "    This optional boolean argument specifies the order in which" \
  " equations are solved by the program.  The {TOPO_SORT} may be 'T' or 1 to" \
  " solve in order of increasing equation weight, or 'F' of 0 to solve in" \
  " scan line orrder.  If omitted, the program assumes \"-sortSys F\".\n" \
  "\n" \
  "  -verbose\n" \
  "    Causes the program to write various diagnostics to stderr, in" \
  " particular at selected iterations of the Gauss-Seidel solver.\n" \
  "\n" \
  "  -reportStep {REPORT_STEP}\n" \
  "    If this optional argument is given and {REPORT_STEP} is" \
  " not zero, the height map {Z}, the error maps {E}, and the" \
  " error summary are written out after every {REPORT_STEP} iterations, as" \
  " well as at the beginning (0 iterations) and after the end of the" \
  " iterative solving.  If {REPORT_STEP} is larger than {MAX_ITER} only" \
  " the initial and final states are written.  If {REPORT_STEP} is zero (the default)," \
  " no \"...-{ITER}-...\" files are written.\n" \
  "\n" \
  "  -outPrefix {PREFIX}\n" \
  "    Specifies the common prefix for all output" \
  " file names.  Mandatory."

#include <stdint.h>
#include <stdio.h>
#include <sys/stat.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <values.h>
#include <math.h>

#include <bool.h>
#include <i2.h>
#include <vec.h>
#include <jsstring.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <argparser.h>
#include <float_image.h>
#include <float_image_mscale.h>
#include <float_image_expand_by_one.h>

#include <pst_map.h>
#include <pst_slope_map.h>
#include <pst_normal_map.h>
#include <pst_height_map.h>
#include <pst_map_compare.h>
#include <pst_cell_map_clear.h>
#include <pst_vertex_map_clear.h>
#include <pst_imgsys.h>
#include <pst_integrate.h>
#include <pst_integrate_iterative.h>

#include <pst_integrate_recursive.h>

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { char *slopes_file;        /* File name of input slope map file, or NULL. */
    char *normals_file;       /* File name of input normal map file, or NULL. */
    float slopes_scaleX;      /* Scale factor for the {X} slope. */
    float slopes_scaleY;      /* Scale factor for the {Y} slope. */
    char *reference_file;     /* File name of input reference height map (NULL if none). */
    float reference_scale;    /* Scale factor for the reference heights. */
    char *hints_file;         /* File name of hints map (external height estimates) (NULL if none). */
    float hints_scale;        /* Scale factor for the external heights. */
    double hints_weight;      /* Weight multiplier for {H} weights. */
    i2_vec_t clear_min;       /* Lower corners of rectangles to be cleared. */
    i2_vec_t clear_max;       /* Upper corners of rectangles to be cleared. */
    char *outPrefix;          /* Output file name prefix. */
    char *initial_opt;        /* Initial option: "zero" or "hints" or "reference". */
    double initial_noise;     /* Relative magnitude of perturbations to the initial map. */
    uint32_t maxLevel;        /* Max recursion level. */
    uint32_t maxIter;         /* Max iterations per level. */
    double convTol;           /* Convergence threshold for top level. */
    bool_t sortSys;           /* TRUE solves the equations in order of increasing weight. */
    /* Debugging and analysis */
    bool_t verbose;           /* TRUE to print various diagnostics. */
    uint32_t reportStep;      /* Frequency for debugging output during iteration, or 0 if none. */
  } options_t;

float_image_t *tire_read_fni_file(char *fileName, int32_t NX, int32_t NY, int32_t wch, bool_t verbose);
  /* Reads a FNI image from file "{fileName}" (which should include the extension ".fni").
    If {fileName} is "-", reads from standard input. 
    
    If {NX} and {NY} are positive, checks that the map read from the file has 
    {NX} cols and {NY} rows.  As a convenience, if the map is one col
    and one row short of {NX} by {NY, enlarges it if {float_image_expand_by_one}.
    In this case, if {wch} is a valid channel index, that channel is treated as 
    a reliability weight channel for the other channels.
    
    IF {NX} and {NY} are negative, the map read from the file is returned
    without size check or resizing.  */

void tire_write_fni_file(float_image_t *I, char *fileName, int32_t indent);
  /* Writes the float image {I} in a human-readable format, to a file
    called "{fileName}" (which should include the extension ".fni").
    Diagnostic images are indented by {indent} spaces. */

void tire_compute_and_write_height_map
  ( options_t *o,
    float_image_t *G, /* Slope map. */
    float_image_t *H, /* External hints map, or {NULL}. */
    float_image_t *Z, /* Initial guess and final height map. */
    float_image_t *R  /* Refernce height map for comparisons, or {NULL}. */
  );
  /* Computes the height map {Z} from the gradient map {G} and optional
    hints map {H}.

    The map {G} must have two or three channels. Channels 0 and 1 are
    intrepreted as the slopes {dZ/dX} and {dZ/dY}). Channel 2, if it
    exists, is interpreted as the reliability weight of the
    corresponding slope data; otherwise the weight is assumed to be 1.

    The hints height map {H}, if not {NULL}, must have one or two
    channels, and one row and one col more than {G}. On input, it gives
    an independent estimate of the heights, and also the initial guess
    for the iterative system solver.

    The image {Z} must have two channels, and one row and one col more
    than {G}. On input, sample {Z[0,X,Y]} should contain the initial
    guess of the height at grid corner {(X,Y)} for the iterative solver;
    if not finite, it is taken as zero. The value of {Z[1,X,Y]} is
    ignored. On output, {Z[0,X,Y]} will be the computed height, and
    {Z[1,X,Y]} the corresponding estimated reliability weight.

    If {o.reportStep} is not zero, writes the height map {Z} for the
    first {o.reportStep} iterations, and then every {o.reportStep}
    iterations thereafter.

    If {R} is not null, it must be a one- or two-channel image with the same
    size as {G} or {Z}.  In this case, whenever {Z} is reported, the procedure computes
    {E=Z-R} with {}. */

options_t *tire_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

void tire_parse_scale(argparser_t *pp, float *scale1_P, float *scale2_P);
  /* If the next command line argument is the keyword "scale", parses
    it. Then for each {P} in {scale1_P} and {scale2_P}, if that {P}
    is not {NULL}, parses a float value and stores it in {*P}. 
    If the next command line argument is not "scale", for each
    {P} in {scale1_P} and {scale2_P}, sets {*P} to 1.0. */

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    options_t *o = tire_parse_options(argc, argv);

    float_image_t *G;  /* Input gradient map. */
    if (o->slopes_file != NULL)
      { fprintf(stderr, "reading the slope map {G} ...\n");
        G = tire_read_fni_file(o->slopes_file, -1, -1, 2, TRUE);
      }
    else if (o->normals_file != NULL)
      { fprintf(stderr, "reading the normal map {N} ...\n");
        float_image_t *N = tire_read_fni_file(o->normals_file, -1, -1, 3, TRUE);
        fprintf(stderr, "converting the normal map to a slope map {G} ...\n");
        double maxSlope = 1000.0; /* Should be enough... */
        G = pst_normal_map_to_slope_map(N, maxSlope);
        float_image_free(N);
      }
    else
      { assert(FALSE); }
    if (o->slopes_scaleX != 1.0) { float_image_rescale_samples(G, 0, 0,1, 0,o->slopes_scaleX); }
    if (o->slopes_scaleY != 1.0) { float_image_rescale_samples(G, 1, 0,1, 0,o->slopes_scaleY); }
    pst_map_ensure_pixel_consistency(G, 2);

    fprintf(stderr, "allocating the height map {Z} ...\n");
    int32_t NX_G, NY_G;
    float_image_get_size(G, NULL, &NX_G, &NY_G);
    int32_t NC_Z = 2;
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;
    float_image_t *Z = float_image_new(NC_Z, NX_Z, NY_Z);

    float_image_t *H; /* Input hints map. */
    if (o->hints_file != NULL)
      { fprintf(stderr, "reading the hints height map {H} ...\n");
        H = tire_read_fni_file(o->hints_file, NX_Z, NY_Z, 1, TRUE);
        if (o->hints_scale != 1.0) { float_image_rescale_samples(H, 0, 0,1, 0,o->hints_scale); }
        pst_map_ensure_pixel_consistency(H, 1);
      }

    float_image_t *R; /* Reference height map. */
    if (o->reference_file == NULL)
      { R = NULL; }
    else
      { fprintf(stderr, "reading the reference height map {R} ...\n");
        R = tire_read_fni_file(o->reference_file, NX_Z, NY_Z, 1, TRUE);
        if (o->reference_scale != 1.0) { float_image_rescale_samples(R, 0, 0,1, 0,o->reference_scale); }
        pst_map_ensure_pixel_consistency(R, 1);
      }

    tire_compute_and_write_height_map(o, G, H, Z, R);

    /* Free working storage: */
    if (G != NULL) { float_image_free(G); G = NULL; }
    if (H != NULL) { float_image_free(H); H = NULL; }
    if (Z != NULL) { float_image_free(Z); Z = NULL; }
    if (R != NULL) { float_image_free(R); R = NULL; }

    fprintf(stderr, "Done!\n");
    return 0;
  }

void tire_compute_and_write_height_map
  ( options_t *o,
    float_image_t *G,
    float_image_t *H,
    float_image_t *Z,
    float_image_t *R
  )
  {
    /* Get and check the slope map size: */
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(G, &NC_G, &NX_G, &NY_G);
    if (o->verbose) { fprintf(stderr, "slope map has %d channels, size = %d × %d\n", NC_G, NX_G, NY_G); }
    demand((NC_G == 2) || (NC_G == 3), "gradient map {G} must have 2 or 3 channels");

    /* Get and check the height map size: */
    int32_t NC_Z, NX_Z, NY_Z;
    float_image_get_size(Z, &NC_Z, &NX_Z, &NY_Z);
    if (o->verbose) { fprintf(stderr, "height map has %d channels, size = %d × %d\n", NC_Z, NX_Z, NY_Z); }
    demand(NC_Z == 2, "height map {Z} must have 2 channels");
    demand((NX_Z == NX_G+1) && (NY_Z == NY_G+1), "bad height map {Z}, size");
    pst_cell_map_clear(G, 2, &(o->clear_min), &(o->clear_max));

    /* Check the hints height map size: */
    if (H != NULL)
      { int32_t NC_H, NX_H, NY_H;
        float_image_get_size(H, &NC_H, &NX_H, &NY_H);
        if (o->verbose) { fprintf(stderr, "hints height map {H} has %d channls, size = %d × %d\n", NC_H, NX_H, NY_H); }
        demand((NC_H == 1) || (NC_H == 2), "hints height map {H} must have 1 or 2 channels");
        demand((NX_H == NX_Z) && (NY_H == NY_Z), "bad hints height map {H}, wrong size");
        pst_vertex_map_clear(H, 1, &(o->clear_min), &(o->clear_max));
      }

    /* Check the reference height map {R}: */
    if (R != NULL)
      { int32_t NC_R, NX_R, NY_R;
        float_image_get_size(R, &NC_R, &NX_R, &NY_R);
        if (o->verbose) { fprintf(stderr, "reference height map {R} has %d channls, size = %d × %d\n", NC_R, NX_R, NY_R); }
        demand((NC_R == 1) || (NC_R == 2), "reference height map {R} must have 1 or 2 channels");
        demand ((NX_R == NX_Z) && (NY_R == NY_Z), "wrong ref heights map {R} size");
        pst_vertex_map_clear(R, 1, &(o->clear_min), &(o->clear_max));
      }

    /* Set {Z} to the specified initial guess: */
    if (strcmp(o->initial_opt, "zero") == 0)
      { fprintf(stderr, "zeroing the initial solution ...\n");
        float_image_fill_channel(Z, 0, 0.0);
      }
    else if (strcmp(o->initial_opt, "reference") == 0)
      { fprintf(stderr, "using the ref height map as the initial solution ...\n");
        demand(R != NULL, "reference map not spacified");
        float_image_assign_channel_rectangle(Z, 0, 0,NX_Z-1, 0,NY_Z-1, R, 0, 0,0);
      }
    else if (strcmp(o->initial_opt, "hints") == 0)
      { fprintf(stderr, "using the hints height map as the initial solution ...\n");
        demand(H != NULL, "hints height map not spacified");
        float_image_assign_channel_rectangle(Z, 0, 0,NX_Z-1, 0,NY_Z-1, H, 0, 0,0);
      }
    else
      { demand(FALSE, "invalid \"-initial\" option"); }
    float_image_fill_channel(Z, 1, 1.0);

    if (o->initial_noise > 0.0)
      { fprintf(stderr, "perturbing the initial guess ...\n");
        double relNoise = o->initial_noise;
        double absNoise = o->initial_noise*fmin(NX_Z, NY_Z)/2;
        pst_height_map_perturb(Z, 1, relNoise, absNoise);
      }

    { fprintf(stderr, "writing out the initial guess ...\n");
      char *fname_ini = jsprintf("%s-ini-Z.fni", o->outPrefix);
      tire_write_fni_file(Z, fname_ini, 0);
      free(fname_ini);
    }

    char *debugDir = jsprintf("%s-iters", o->outPrefix); /* Folder for debug files. */
    mkdir(debugDir, 0755);
    char *debugPrefix = jsprintf("%s/it-", debugDir); /* Prefix for debug file names. */

    auto void reportSys(int32_t level, pst_imgsys_t *S);
      /* This procedure is called by once per level on the way up.
        to report the linear system created from the slope and weight maps. */

    auto void reportHeights
      ( int32_t level,
        int32_t iter,
        double change,
        bool_t final,
        float_image_t *cur_Z,
        float_image_t *cur_R
      );
      /* This procedure is called by {pst_integrate_iterative} to write
        the current height map {cur_Z}. It is called at least once per level on the way up,
        including once with {final=TRUE}.

        The parameter {cur_R} is {NULL} or the reference map {R} scaled
        down to the current level. If {cur_R} is not {NULL}, the
        procedure computes and writes the height error map {cur_E}, and
        its summary */

    auto void reportData
      ( int32_t level,
        float_image_t *cur_G,
        float_image_t *cur_H,
        float_image_t *cur_R
      );
      /* This procedure is called at various times during the recursive
         multiscale integration. See {pst_integrate_recursive}.  */

    /* Call recursive integrator: */
    int32_t level = 0;
    pst_integrate_recursive
      ( level, G, H, o->hints_weight, Z, R,
        o->maxLevel, o->maxIter, o->convTol, o->sortSys,
        o->verbose,
        &reportData,
        &reportSys,
        o->reportStep,
        &reportHeights
      );

    free(debugDir);
    free(debugPrefix);

    return;

    void reportData
      ( int32_t level,
        float_image_t *cur_G,
        float_image_t *cur_H,
        float_image_t *cur_R
      )
      {
        /* This procedure is called once per level, on the way down.*/

        assert((level >= 0) && (level <= o->maxLevel));
        int32_t indent = (level < -1 ? 0 : 2*level+2);

        int32_t NX_cur_G, NY_cur_G;
        float_image_get_size(cur_G, NULL, &NX_cur_G, &NY_cur_G);
        if (o->verbose) { fprintf(stderr, "%*swriting the current slope map {G} %d×%d ...\n", indent, "", NX_cur_G, NY_cur_G); }
        float_image_mscale_write_file(cur_G, o->outPrefix, level, 0, "G");

        if (level == 0) { assert(cur_H == H); }
        if (cur_H != NULL)
          { int32_t NX_cur_H, NY_cur_H;
            float_image_get_size(cur_H, NULL, &NX_cur_H, &NY_cur_H);
            if (o->verbose) { fprintf(stderr, "%*swriting the current hints map {H} %d×%d ...\n", indent, "", NX_cur_H, NY_cur_H); }
            assert((NX_cur_H == NX_cur_G + 1) && (NY_cur_H == NY_cur_G + 1));
            float_image_mscale_write_file(cur_H, o->outPrefix, level, 0, "H");
          }

        if (level == 0) { assert(cur_R == R); }
        if (cur_R != NULL)
          { int32_t NX_cur_R, NY_cur_R;
            float_image_get_size(cur_R, NULL, &NX_cur_R, &NY_cur_R);
            if (o->verbose) { fprintf(stderr, "%*swriting the current reference map {R} %d×%d ...\n", indent, "", NX_cur_R, NY_cur_R); }
            assert((NX_cur_R == NX_cur_G + 1) && (NY_cur_R == NY_cur_G + 1));
            float_image_mscale_write_file(cur_R, o->outPrefix, level, 0, "R");
          }
      }

    void reportSys(int32_t level, pst_imgsys_t *cur_S)
      { assert((level >= 0) && (level <= o->maxLevel));
        int32_t indent = (level < -1 ? 0 : 2*level+2);

        if (cur_S == NULL)
          { if (o->verbose) { fprintf(stderr, "%*sno equation system was generated at this level\n", indent, ""); } }
        else
          { if (o->verbose) { fprintf(stderr, "%*swriting out the equations system %d×%d ...\n", indent, "", cur_S->N, cur_S->N); }
            char *fileName_sys = float_image_mscale_file_name(o->outPrefix, level, 0, "S", "txt");
            FILE* wr_sys = fopen(fileName_sys, "wt");
            pst_imgsys_write(wr_sys, cur_S, "%+10.6f");
            fclose(wr_sys);
            free(fileName_sys);

            if (o->verbose) { fprintf(stderr, "%*swriting out the equations system weight image ...\n", indent, ""); }
            float_image_t *SW = pst_imgsys_make_weight_image(cur_S);
            float_image_mscale_write_file(SW, o->outPrefix, level, 0, "SW");
            float_image_free(SW);
          }
      }

    void reportHeights
      ( int32_t level,
        int32_t iter,
        double change,
        bool_t final,
        float_image_t *cur_Z,
        float_image_t *cur_R
      )
      { assert(level >= 0);
        pst_height_map_analyze_and_write
          ( (final || (iter == 0) ? o->outPrefix : debugPrefix), 
            level, (final ? -1 : iter), change,
            cur_Z, cur_R, o->verbose
          );
      }
  }

float_image_t *tire_read_fni_file(char *fileName, int32_t NX, int32_t NY, int32_t wch, bool_t verbose)
  { demand(fileName != NULL, "file name not given");
    fprintf(stderr, "Reading %s ...\n", fileName);
    FILE *rd = open_read(fileName, FALSE);
    float_image_t *I = float_image_read(rd);
    if (rd != stdin) { fclose(rd); }
    if ((NX != -1) || (NY != -1))
      { int32_t NC_I, NX_I, NY_I;
        float_image_get_size(I, &NC_I, &NX_I, &NY_I);
        if (verbose) { fprintf(stderr, "given map size is %d × %d\n", NX_I, NY_I); }
        if ((NX_I == NX - 1) && (NY_I == NY - 1))
          { /* Seems that a cell map was giveninstead of a vertex map. */
            /* Resize to {NX, NY}: */
            if (verbose) { fprintf(stderr, "Expanding the map to the correct size ...\n"); }
            float_image_t *II = float_image_expand_by_one(I, 1);
            float_image_free(I);
            I = II;
          }
      }
    return I;
  }

void tire_write_fni_file(float_image_t *I, char *fileName, int32_t indent)
  { demand(fileName != NULL, "file name not given");
    FILE* wr = open_write(fileName, FALSE);
    float_image_write(wr, I);
    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
    fprintf(stderr, "%*swrote %s\n", indent, "", fileName);
  }

options_t *tire_parse_options(int32_t argc, char **argv)
  { argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");

    if (argparser_keyword_present(pp, "-slopes"))
      { o->slopes_file = argparser_get_next_non_keyword(pp); }
    else
      { o->slopes_file = NULL; }

    if (argparser_keyword_present(pp, "-normals"))
      { o->normals_file = argparser_get_next_non_keyword(pp); }
    else
      { o->normals_file = NULL; }

    if ((o->slopes_file == NULL) == (o->normals_file == NULL))
      { argparser_error(pp, "must specify exaclty one of \"-slopes\" or \"-normals\""); }

    tire_parse_scale(pp, &(o->slopes_scaleX), &(o->slopes_scaleY));

    if (argparser_keyword_present(pp, "-hints"))
      { o->hints_file = argparser_get_next_non_keyword(pp);
        tire_parse_scale(pp, &(o->hints_scale), NULL);
        o->hints_weight = argparser_get_next_double(pp, 0.0, 1.0e10);
      }
    else
      { o->hints_file = NULL; o->hints_weight = 0.0; }

    if (argparser_keyword_present(pp, "-reference"))
      { o->reference_file = argparser_get_next_non_keyword(pp);
        tire_parse_scale(pp, &(o->reference_scale), NULL);
      }
    else
      { o->reference_file = NULL; }

    argparser_get_keyword(pp, "-initial");
    o->initial_opt = argparser_get_next_non_keyword(pp);
    o->initial_noise = argparser_get_next_double(pp, 0, DBL_MAX);
    
    o->clear_min = i2_vec_new(0);
    o->clear_max = i2_vec_new(0);
    int32_t nclear = 0;
    while (argparser_keyword_present(pp, "-clear"))
      { int32_t clx_min = (int32_t)argparser_get_next_int(pp, 0, INT64_MAX);
        int32_t clx_max = (int32_t)argparser_get_next_int(pp, 0, INT64_MAX);
        int32_t cly_min = (int32_t)argparser_get_next_int(pp, 0, INT64_MAX);
        int32_t cly_max = (int32_t)argparser_get_next_int(pp, 0, INT64_MAX);
        i2_vec_expand(&(o->clear_min), nclear);
        o->clear_min.e[nclear] = (i2_t){{ clx_min, cly_min }};
        i2_vec_expand(&(o->clear_max), nclear);
        o->clear_max.e[nclear] = (i2_t){{ clx_max, cly_max }};
        nclear++;
      }
    i2_vec_trim(&(o->clear_min), (uint32_t)nclear);
    i2_vec_trim(&(o->clear_max), (uint32_t)nclear);
        
    if (argparser_keyword_present(pp, "-maxLevel"))
      { o->maxLevel = (uint32_t)argparser_get_next_int(pp, 0, INT64_MAX); }
    else
      { o->maxLevel = DEFAULT_MAX_ITER; }

    if (argparser_keyword_present(pp, "-maxIter"))
      { o->maxIter = (uint32_t)argparser_get_next_int(pp, 0, INT64_MAX); }
    else
      { o->maxIter = DEFAULT_MAX_ITER; }

    if (argparser_keyword_present(pp, "-convTol"))
      { o->convTol = argparser_get_next_double(pp, 0, DBL_MAX); }
    else
      { o->convTol = DEFAULT_CONV_TOL; }

    if (argparser_keyword_present(pp, "-sortSys"))
      { o->sortSys = argparser_get_next_bool(pp); }
    else
      { o->sortSys = FALSE; }

    if (argparser_keyword_present(pp, "-reportStep"))
      { o->reportStep = (uint32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->reportStep = 0; }

    o->verbose = argparser_keyword_present(pp, "-verbose");

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);

    argparser_finish(pp);

    return o;
  }

void tire_parse_scale(argparser_t *pp, float *scale1_P, float *scale2_P)
  {
    if (argparser_keyword_present_next(pp, "scale"))
      { if (scale1_P != NULL) { (*scale1_P) = (float)argparser_get_next_double(pp, -1.e10, +1.0e10); }
        if (scale2_P != NULL) { (*scale2_P) = (float)argparser_get_next_double(pp, -1.e10, +1.0e10); }
      }
    else
      { if (scale1_P != NULL) { (*scale1_P) = 1.0; }
        if (scale2_P != NULL) { (*scale2_P) = 1.0; }
      }
  }
