#define PROG_NAME "test_build_system"
#define PROG_DESC "test the building of a system of integration equations from a slope map"
#define PROG_VERS "1.0"

#define slope_to_height_C_COPYRIGHT "Copyright © 2024 by the State University of Campinas (UNICAMP)"

/* Last edited on 2025-01-15 06:46:22 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -slopes {IG_FNI_NAME} \\\n" \
  "    -weights {IW_FNI_NAME} \\\n" \
  "    -outPrefix {PREFIX} \\\n" \
  "    [ -full {FULL} ] \\\n" \
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
  "  Created 2025-01-06 by Jorge Stolfi, UNICAMP, based on {slope_to_height.c} by R. Saracchini.\n" \
  "MODIFICATION HISTORY\n" \
  "  2025-01-06 Created.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " slope_to_height_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define stringify(x) strngf(x)
#define strngf(x) #x

#define PROG_INFO_DESC \
  "  The program reads a /slope map/ {IG}, consisting of" \
  " the {X} and {Y} derivatives of some terrain surface {Z(X,Y)}; builds" \
  " the corresponding (overdetermined) system of equations directly fron iy; and" \
  " writes that system out.\n" \
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
  " datum should be ignored."
  
#define PROG_INFO_OUT_FILES \
  "  {PREFIX}-S.sys\n" \
  "    The integration system.\n" \
  "\n" \
  "  {PREFIX}-W.png\n" \
  "    A single-channel image showing the total" \
  " weights {eq[x,y].wtot} of the equation of the" \
  " integration system referring to each pixel {x,y}."

#define PROG_INFO_OPTS \
  "  -slopes {IG_FNI_NAME}\n" \
  "    This mandatory argument specifies the file name" \
  " containing the input slope map {IG}, which must be a two-channel" \
  " image in FNI format" \
  " (see float_image.h).  If {IG_FNI_NAME}" \
  " is \"-\", the slope map is read from standard input.\n" \
  "\n" \
  "  -weights {IW_FNI_NAME}\n" \
  "    This optional argument specifies the file containing the input" \
  " slope weight map {IW}.  It must be a single-channel image with non-negative" \
  " values, preferably between 0 and 1, and with the same size as {IG}.  Slope" \
  " samples with weight 0 are ignored.  If {IW_FNI_NAME}" \
  " is \"-\", the slope weight map is read from standard input.  If this" \
  " argument is omitted, the" \
  " program assumes that all slope samples have weight 1.\n" \
  "\n" \
  "  -full {FULL}\n" \
  "    This optional argument specifies the {full} argument" \
  " to {pst_integrate_build_system} (q.v.).  The {FULL} value" \
  " can be 0, 1, 'F', or 'T'.  If omitted, assumes \"-full F\".\n" \
  "\n" \
  "  -outPrefix {PREFIX}\n" \
  "    Specifies the common prefix for all output" \
  " file names.  Mandatory."

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdlib.h>
#include <values.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <jsstring.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <argparser.h>
#include <float_image.h>
#include <float_image_mscale.h>
#include <float_image_write_gen.h>

#include <pst_slope_map.h>
#include <pst_imgsys.h>
#include <pst_integrate.h>

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  {
    char *slopes_fname;   /* File name of input slope file, or "-" for stdin. */
    char *weights_fname;  /* File name of weight slope file, or "-" for stdin (NULL if none). */
    bool_t full;          /* True to include equations for isolated height pixels. */
    char *outPrefix;      /* Output file name prefix. */
  } options_t;

float_image_t *read_fni_file(char *fname);
  /* Reads a FNI image from file "{fname}" (which should include the extension ".fni").
    If {fname} is "-", reads from standard input. */

void write_system(pst_imgsys_t *S, char *fname);
  /* Writes the system {S} to file "{fname}". If {fname} is "-", writes
    to standard output. */
    
void write_system_weights(pst_imgsys_t *S, char *fname);
  /* Writes the equation weights {S.eq[k].wtot} as a PNG image "{fname}". */

options_t *parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    fprintf(stderr, "Reading the slope map {IG}:\n");
    float_image_t *IG;  /* Input gradient map. */
    IG = read_fni_file(o->slopes_fname);
    
    float_image_t *IW; /* Input (gradient) reliability weight map. */
    if (o->weights_fname == NULL)
      { IW = NULL; }
    else
      { fprintf(stderr, "Reading the slope reliability weight map {IW}:\n");
        IW = read_fni_file(o->weights_fname);
      }
    
    fprintf(stderr, "Building the integration system:\n");
    bool_t full = TRUE;
    pst_imgsys_t *S = pst_integrate_build_system(IG, IW, full);
    
    char *system_fname = jsprintf("%s-S.txt", o->outPrefix);
    write_system(S, system_fname);
    free(system_fname);
    
    char *weights_fname = jsprintf("%s-W.png", o->outPrefix);
    write_system_weights(S, weights_fname);
    free(weights_fname);
    
    float_image_free(IG); IG = NULL;
    if (IW != NULL) { float_image_free(IW); IW = NULL; }
    pst_imgsys_free(S); S = NULL;

    fprintf(stderr, "Done!\n");
    return 0;
  }

float_image_t *read_fni_file(char *fname)
  { demand(fname != NULL, "file name not given");
    FILE *rd = open_read(fname, TRUE);
    float_image_t *I = float_image_read(rd);
    if (rd != stdin) { fclose(rd); }
    fprintf(stderr, "\n");
    return I;
  }

void write_system(pst_imgsys_t *S, char *fname)
  { demand(fname != NULL, "file name not given");
    FILE* wr = open_write(fname, TRUE);
    pst_imgsys_write(wr, S);
    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
    fprintf(stderr, "\n");
  }

void write_system_weights(pst_imgsys_t *S, char *fname)
  { demand(fname != NULL, "file name not given");
    float_image_t *U = float_image_new(1, S->NX, S->NY);
    pst_imgsys_extract_system_eq_tot_weight_image(S, U, 0.0);
    double expoEnc = 1.0;
    double bias = 0.0;
    bool_t verbose = TRUE;
    float_image_write_gen_named(fname, U, image_file_format_PNG, 0.0, 1.0, expoEnc, bias, verbose);
    fprintf(stderr, "\n");
  }

options_t *parse_options(int32_t argc, char **argv)
  { argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = talloc(1, options_t);
    
    if (argparser_keyword_present(pp, "-full"))
      { o->full = argparser_get_next_bool(pp); }
    else
      { o->full = FALSE; }

    argparser_get_keyword(pp, "-slopes");
    o->slopes_fname = argparser_get_next_non_keyword(pp);
    
    if (argparser_keyword_present(pp, "-weights"))
      { o->weights_fname = argparser_get_next(pp); }
    else
      { o->weights_fname = NULL; }
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    
    argparser_finish(pp);
    
    return o;
  }
