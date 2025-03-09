#define PROG_NAME "test_build_system"
#define PROG_DESC "test the building of a system of integration equations from a slope map"
#define PROG_VERS "1.0"

#define slope_to_height_C_COPYRIGHT "Copyright © 2024 by the State University of Campinas (UNICAMP)"

/* Last edited on 2025-03-03 03:47:29 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -slopes {IG_FNI_NAME} \\\n" \
  "    [ -hints {IH_FNI_NAME} ] \\\n" \
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
  "  Created 2025-01-06 by Jorge Stolfi, UNICAMP, based on {slope_to_height.c} by R. Saracchini.\n" \
  "MODIFICATION HISTORY\n" \
  "  2025-01-06 Created.\n" \
  "  2025-02-11 Modified to use {Z} hint map {IH}, option \"-hints\".\n" \
  "  2025-02-11 Replaced \"-weights\" by weights embedded in {IG,IH}.\n" \
  "  2025-02-11 Removed option \"-full\".\n" \
  "  2025-02-11 Removed option to read {stdin}, write {stdout}.\n" \
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
  "  The program reads a /slope map/ image {IG} with" \
  " the {X} and {Y} derivatives of some terrain surface {Z(X,Y)}; and" \
  " a /hints map/ {IH} with external" \
  " estimates of {Z(X,Y)}.  It builds" \
  " the corresponding (overdetermined) system of equations that" \
  " estimates the heights {Z(X,Y)} from {IG} and {IH}; and" \
  " writes that system out.\n" \
  "\n" \
  "  The input slope map {IG} is a three-channel float image.  Each" \
  " pixel {IG[X,Y]} of this image has samples {(xG,yG,wG)}, where " \
  " {(xG,yG)} is the gradient of the {Z} function, averaged" \
  " over the unit square cell (/pixel/) whose lower left corner" \
  " is the point {(X,Y)}; and {wG} is the reliability weight of that gradient.\n" \
  "\n" \
  "  The input hints map {IH} is a two-channel float image.  Each" \
  " pixel {IH[X,Y]} of this image has samples {(zH,wH)}, where " \
  " {zH} is some independent estimate of the {Z} function at the pixel corner" \
  " {(X,Y)}; and {wH} is the reliability weight of that estimate.  The {IH} image" \
  " must therefore have one col and one row more than the {IG} image."
  
#define PROG_INFO_OUT_FILES \
  "  {PREFIX}-S.sys\n" \
  "    The integration system.\n" \
  "\n" \
  "  {PREFIX}-W.png\n" \
  "    A single-channel image showing the total" \
  " weights {eq[x,y].wtot} of the equation of the" \
  " integration system referring to each pixel corner {x,y}."

#define PROG_INFO_OPTS \
  "  -slopes {IG_FNI_NAME}\n" \
  "    This mandatory argument specifies the file name" \
  " containing the input slope map {IG}.  It must be a three-channel" \
  " image in FNI format" \
  " (see float_image.h).\n" \
  "\n" \
  "  -hints {IH_FNI_NAME}\n" \
  "    This optional argument specifies the file name" \
  " containing the hints map {IH}, with the independent {Z}" \
  " estimates. It must be a two-channel" \
  " image in FNI format" \
  " (see float_image.h).\n" \
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
    char *slopes;         /* File name of input slope map. */
    char *hints;          /* File name of input {Z} hints map, or NULL. */
    char *outPrefix;      /* Output file name prefix. */
  } options_t;

float_image_t *read_fni_file(char *fname);
  /* Reads a FNI image from file "{fname}" (which should include the extension ".fni"). */

void write_system(pst_imgsys_t *S, char *fname);
  /* Writes the system {S} to file "{fname}". */

void write_fni_file(char *fname, float_image_t *img);
  /* Writes {img} in FNI format as "{fname}" (which should include the extension ".fni"). */

options_t *parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    fprintf(stderr, "Reading the slope map {IG}:\n");
    float_image_t *IG = read_fni_file(o->slopes);
    
    float_image_t *IH = NULL;
    if (o->hints != NULL)
      { fprintf(stderr, "Reading the {Z} hints map {IH}:\n");
        IH = read_fni_file(o->hints);
      }
    
    fprintf(stderr, "Building the integration system:\n");
    bool_t verbose_build = TRUE;
    double wHMult = 0.1;
    int32_t indent = 6;
    pst_imgsys_t *S = pst_integrate_build_system(IG, IH, wHMult, indent, verbose_build);
    
    char *system_fname = jsprintf("%s-S.txt", o->outPrefix);
    write_system(S, system_fname);
    free(system_fname);
    
    float_image_t *SW = pst_imgsys_make_weight_image(S);

    char *wtot_fname = jsprintf("%s-SW.fni", o->outPrefix);
    write_fni_file(wtot_fname, SW);
    free(wtot_fname);
    
    float_image_free(IG); IG = NULL;
    if (IH != NULL) { float_image_free(IH); IH = NULL; }
    pst_imgsys_free(S); S = NULL;
    float_image_free(SW); SW = NULL;
    
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
    pst_imgsys_write(wr, S, "%+10.6f");
    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
    fprintf(stderr, "\n");
  }

void write_fni_file(char *fname, float_image_t *img)
  { demand(fname != NULL, "file name not given");
    float_image_write_named(fname, img);
  }

options_t *parse_options(int32_t argc, char **argv)
  { argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = talloc(1, options_t);

    argparser_get_keyword(pp, "-slopes");
    o->slopes = argparser_get_next_non_keyword(pp);

    if (argparser_keyword_present(pp, "-hints"))
      { o->hints = argparser_get_next_non_keyword(pp); }
    else
      { o->hints = NULL; }
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    
    argparser_finish(pp);
    
    return o;
  }
