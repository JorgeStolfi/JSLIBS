#define PROG_NAME "make_test_maps"
#define PROG_DESC "generate gradient and height maps for tests of slope_to_height"
#define PROG_VERS "1.0"

/* Copyright © 2005 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2025-03-04 18:33:35 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -function {ZFUNC} \\\n" \
  "    -size {NX} {NY} \\\n" \
  "    -perturbG {PERTURB_G} \\\n" \
  "    -outDir {OUT_DIR}" " \\\n" \
  "    " argparser_help_info_HELP
  
#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Generates a height map {Z(X,Y)} and its gradient (slope)" \
  " map {G(X,Y)}, inclunding" \
  " the corresponding realiability weights, in the format required" \
  " by {slope_to_height}.\n" \
  "\n" \
  "  The height map is defined by the {ZFUNC} parameter, an integer" \
  " that selects one of the program's built-in functions" \
  " ({function_00}, {function_01}, ...). See {pst_proc_map.h} for" \
  " the list of available functions.  The parameters" \
  " {NX} and {NY} are the width and height of the slope" \
  " map, in pixels.\n" \
  "\n" \
  "  The height function is nominally evaluated" \
  " at the points (corners) of a rectangular grid with" \
  " {NX} by {NY} square cells.  The value is the average of the" \
  " heights at several sampling points (sampoints) surrounding" \
  " each corner, weighted by a suitable window weight" \
  " table.  Note that the height map" \
  " has {NX+1} columns and {NY+1} rows.\n" \
  "\n" \
  "  The slope map is nominally evaluated at the center of each cell (pixel)" \
  " of the grid, so it has {NX} columns and {NY} rows.  Two versions" \
  " of the slope map are computed, 'numeric' and 'analytic'. The numeric" \
  " gradient in each grid cell is computed by divided" \
  " differences of the height values computed a collection sampoints" \
  " spanning the cell and part of adjacent cells. The analytic" \
  " gradient is computed by sampling and averaging the analytic" \
  " gradient of the height function at those sampoints.  Note that" \
  " the analytic slope map may be" \
  " very different from the numeric one, and inconsistent" \
  " with the height map, especially if the function" \
  " is discontinuous or highly random inside the pixel.\n" \
  "\n" \
  "   If {RAW_COORDS} is true, the argument{p} passed to the" \
  " height function ranges over the rectangle {[0_NX]×[0_NY]}.  If" \
  " {RAW_COORDS} is false, these coordinates are shifted and scaled" \
  " so that the largest centered square that fits in that" \
  " rectangle is mapped to {[-1 _ +1]^2}.\n" \
  "\n" \
  "  If {PERTURB_G} is positive, the program also randomly perturbs each" \
  " component of the gradient value, to simulate noisy data.  The" \
  " perturbation consists in adding to each gradient component a random number with" \
  " Gaussian distribution, zero mean, and standard" \
  " deviation {PERTURB_G}.  The same random value is added to both gradient" \
  " versions, numeric and anaytic.  The perturbation affects only" \
  " channels 0 and 1 of the gradient map.\n" \
  "\n" \
  "  The program also writes a map with the discrepancy between the" \
  " numeric and analytic gradients at each pixel (before the random noise).\n" \
  "\n" \
  "  The normal map is computed from the (perturbed) numeric slope map.  The Z"\
  " component is always positive.\n" \
  "\n" \
  "  Each map includes a weight channel, whose value at each pixel gives" \
  " the reliability of the corresponding height, gradient, or normal" \
  " in the same pixel.  The weight depends on the map values and the" \
  " assumed implicit noise levels {MAX_GRAD} and {MAX_GDIFF}. See" \
  " {pst_proc_map_make_height_map} and" \
  " {pst_proc_map_make_slope_map} for details.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  The names of output files begin with {OUT_PREFIX} which" \
  " is \"{OUT_DIR}/{ZFUNC}-{FUNC_NAME}-{NOISY}-{XXXX}x{YYYY}\", where" \
  " {NOISY} is 'N' or 'Y' depending on {PERTURB_G}, and {XXXX} and {YYYY} are" \
  " the parameters {NX} and {NY} formatted as '%04d'.\n" \
  "\n" \
  "  The height map is written as channel 0 of the" \
  " two-channel float image file \"{OUT_PREFIX}-Z.fni\".  Note that" \
  " its size is actually {NX+1} by {NY+1}.  The reliability weights are saved in channel 1.\n" \
  "\n" \
  "  The slope maps are written as" \
  " channels 0 (X derivative) and channel 1 (Y derivative) of" \
  " the three-channel float image files \"{OUT_PREFIX}-Gn.fni\" (numeric)" \
  " and \"{OUT_PREFIX}-Ga.fni\" (analytic).  The discrepancy between them" \
  " is saved as \"{OUT_PREFIX}-Gd.fni\" The reliability weights are saved in channel 2.\n" \
  "\n" \
  "  The normal map is written as channels 0, 1, and 2 of the four-channel float image" \
  " file \"{OUT_DIR}-N.fni\". The reliability weights" \
  " are saved in channel 3." \
  "\n" \
  "SEE ALSO\n" \
  "  fni_to_pnm(1), pnm_to_fni(1), slope_to_height(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2005-08-15 by Jorge Stolfi, UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  By J. Stolfi, IC-UNICAMP unless otherwise noted.\n" \
  "\n" \
  "    2010-05-04 Weight is always 1 with \"-numGrad\".\n" \
  "    2025-01-18 Weight is extra channel in each \".fni\" image file.\n" \
  "    2025-01-19 Removed options \"-smoothZ\", \"-smoothG\", \"-noiseW\".\n" \
  "    2025-01-19 Added option \"-maxGrad\".\n" \
  "    2025-01-19 Made \"-numGrad\" take an argument.\n" \
  "    2025-01-21 Removed options \"-maxGrad\", \"-maxGDiff\".\n" \
  "    2025-01-24 Removed option \"-numGrad\".\n" \
  "    2025-01-24 Wrote both versions of the slope map and the difference.\n" \
  "    2025-02-25 Uses {.rawCoords} attribute to func descriptor.\n" \
  "    2025-02-25 Removed the {.rawCoords} attribute from the func descriptor."

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>

#include <float_image.h>
#include <argparser.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <wt_table_hann.h>
#include <r2.h>
#include <r3.h>

#include <pst_proc_map.h>
#include <pst_height_map.h>
#include <pst_slope_map.h>
#include <pst_normal_map.h>

typedef struct sampling_opts_t
  { int N;     /* Number of sampling points per axis. */
    int L;     /* Width of window in pixels. */
  } sampling_opts_t;
  /* A sampling specification. The window width must be an integer
    to ensure the partition-of-unit property. */

typedef struct options_t
  { int function;            /* Integer code of height map. */
    int NX, NY;              /* Image size. */
    double perturbG;         /* Deviation of gaussian gradient noise. */
    char* outDir;            /* Output directory. */
  } options_t;

/* INTERNAL PROTOTYPES */

void write_test_image(char *pref, int32_t fnum, char *fname, bool_t noisy, int32_t NX, int32_t NY, char *tag, float_image_t *I);
  /* Writes the image {I} to file named
    "{pref}-{NN}-{fname}-{S}-{tag}.fni", in FNI format (see
    float_image.h}); where {NN} is {fnum} formatted as "%02d", and {S}
    is {noisy} formatted as 'N' or 'Y'. */

void normalize(double v[]);
  /* Normalizes the three-vector {v[0..2]} to unit Euclidean length. */

void write_fni_image(char *fileName, float_image_t *I);
  /* Writes the image {I} to file "{fileName}" as a FNI image, with
    {float_image_write}. */

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

sampling_opts_t parse_sampling_options(argparser_t *pp, char *key);
  /* Parses an optional smoothing option consisting of the keyword {key} followed
    by the window width {L} and the number of samples per axis {N}. */

int main(int argc, char** argv);

/* IMPLEMENTATIONS */

int main(int argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    pst_proc_map_zfunc_props_t fp = pst_proc_map_function_generic(o->function);

    /* Compute the sample displacements and weights: */
    uint32_t NS = 5;
    double ws[NS];
    uint32_t stride;
    wt_table_hann_fill(NS, 0.0, ws, &stride);
    assert(stride == NS/2 + 1);

    /* Compute height and gradient maps: */
    double minWeight = 1.0e-4;
    
    int32_t xDebug = -1, yDebug = -1;

    int32_t NGMIN = (o->NX < o->NY ? o->NX : o->NY);  /* Smallest dimension of slope map. */
    double pixSize = 2.0/NGMIN; /* Number of grid pixels per {func} slope map domain unit. */
    r2_t org = (r2_t){{ 0.5*o->NX, 0.5*o->NY }};
    
    float_image_t *IZ = pst_proc_map_make_height_map
      ( fp.func, o->NX, o->NY, &org, pixSize, NS, ws, xDebug, yDebug );
      
    float_image_t *IGn = pst_proc_map_make_slope_map
      ( fp.func, o->NX, o->NY, &org, pixSize, NS, ws, TRUE, fp.maxGrad, fp.maxGDiff, minWeight, xDebug, yDebug );
      
    float_image_t *IGa = pst_proc_map_make_slope_map
      ( fp.func, o->NX, o->NY, &org, pixSize, NS, ws, FALSE, fp.maxGrad, fp.maxGDiff, minWeight, xDebug, yDebug );
      
    float_image_t *IGd = float_image_new(3, o->NX, o->NY);
    for (int32_t y = 0; y < o->NY; y++)
      { for (int32_t x = 0; x < o->NX; x++)
          { for (int32_t c = 0; c < 3; c++)
              { float va = float_image_get_sample(IGa, c, x, y);
                float vn = float_image_get_sample(IGn, c, x, y);
                float vd = (float)(c < 2 ? vn - va : 2.0/(1.0/va + 1.0/vn));
                float_image_set_sample(IGd, c, x, y, vd);
              }
          }
      }
    if (o->perturbG != 0.0)
      { uint32_t seed = 271828183;
        pst_slope_map_perturb(IGa, o->perturbG, seed);
        pst_slope_map_perturb(IGn, o->perturbG, seed);
      }
      
    float_image_t *IG = pst_slope_map_merge(IGa, IGn);
      
    float_image_t *IN = pst_normal_map_from_slope_map(IGn);
    
    /* Write images: */
    bool_t noisy = (o->perturbG != 0.0);
    write_test_image(o->outDir, fp.num, fp.name, noisy, o->NX, o->NY, "Z",  IZ);
    write_test_image(o->outDir, fp.num, fp.name, noisy, o->NX, o->NY, "Ga", IGa);
    write_test_image(o->outDir, fp.num, fp.name, noisy, o->NX, o->NY, "Gn", IGn);
    write_test_image(o->outDir, fp.num, fp.name, noisy, o->NX, o->NY, "Gd", IGd);
    write_test_image(o->outDir, fp.num, fp.name, noisy, o->NX, o->NY, "G",  IG);
    write_test_image(o->outDir, fp.num, fp.name, noisy, o->NX, o->NY, "N",  IN);
    
    return 0;
  }

void write_test_image(char *pref, int32_t fnum, char *fname, bool_t noisy, int32_t NX, int32_t NY, char *tag, float_image_t *I)
  { char *fileName = jsprintf("%s/%02d-%s-%c-%04dx%04d-%s.fni", pref, fnum, fname, "NY"[noisy], NX, NY, tag);
    write_fni_image(fileName, I);
    free(fileName);
  }

void normalize(double v[])
  { int c;
    double m2 = 0;
    for (c = 0; c < 3; c++) { double vax = v[c]; m2 += vax*vax; }
    double m = sqrt(m2);
    for (c = 0; c < 3; c++) 
      { if (m > 0)
          { /* Normalize vector to unit length: */
            v[c] /= m;
          }
        else
          { /* Messy spot, give +Z as the normal: */
            v[c] = (c == 2);
          }
      }
  }

void write_fni_image(char *fileName, float_image_t *I)
  { FILE *wr = open_write(fileName, TRUE);
    float_image_write(wr, I);
    fclose(wr);
  }

options_t *parse_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
   
    options_t *o = (options_t *)malloc(sizeof(options_t));

    argparser_get_keyword(pp, "-function"); 
    o->function = (int)argparser_get_next_int
      ( pp, pst_proc_map_MIN_ZFUNC, pst_proc_map_MAX_ZFUNC );
    
    argparser_get_keyword(pp, "-size"); 
    o->NX = (int)argparser_get_next_int(pp, 1, 4095);
    o->NY = (int)argparser_get_next_int(pp, 1, 4095);
   
    argparser_get_keyword(pp, "-perturbG"); 
    o->perturbG = argparser_get_next_double(pp, 0.0, +DBL_MAX);
    
    argparser_get_keyword(pp, "-outDir"); 
    o->outDir = argparser_get_next(pp);
    
    argparser_skip_parsed(pp);

    argparser_finish(pp);
    
    return o;
  }

sampling_opts_t parse_sampling_options(argparser_t *pp, char *key)
  {
    sampling_opts_t smooth;
    if (argparser_keyword_present(pp, key)) 
      { smooth.L = (int)argparser_get_next_int(pp, 0, 2);
        smooth.N = (int)argparser_get_next_int(pp, 1, (smooth.L == 0 ? 1 : 65));
      }
    else
      { smooth.L = 0;
        smooth.N = 1;
      }
    return smooth;
  }

/* COPYRIGHT, AUTHORSHIP, AND WARRANTY NOTICE:
** 
**   Copyright © 2005 by the State University of Campinas (UNICAMP).
**
** Created on 2005-08-15 by Jorge Stolfi, IC-UNICAMP.       
**
** Permission to use, copy, modify, and redistribute this software and
** its documentation for any purpose and without fee is hereby
** granted, provided that: (1) the copyright notice at the top of this
** file and this copyright, authorship, and warranty notice is retained
** in all derived source files and documentation; (2) no executable
** code derived from this file is published or distributed without the
** corresponding source code; and (3) these same rights are granted to
** any recipient of such code, under the same conditions.
** This software is provided "as is", WITHOUT ANY EXPLICIT OR IMPLICIT
** WARRANTIES, not even the implied warranties of merchantibility and
** fitness for a particular purpose. END OF NOTICE.
*/
