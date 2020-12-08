#define PROG_NAME "test_misc_ops"
#define PROG_DESC "test of various unary ops on {float_image_t} images"
#define PROG_VERS "1.0"

/* Last edited on 2017-06-22 18:10:31 by stolfilocal */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_misc_ops_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    { -dilate | -gradSqr | -gradSqrRel } \\\n" \
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
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  convert(1), gimp(1), display(1), ppm(1), pgm(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on 2009-02-15 by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_misc_ops_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  This program writes to {stdout} the gradient-squared or" \
  " relative-gradient-squared image of a PGM or PPM" \
  " image read from {stdin}."

#define PROG_INFO_OPTS \
  "  None."

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include <argparser.h>
#include <affirm.h>
#include <bool.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <sample_conv.h>
#include <float_image_write_pnm.h>
#include <float_image_read_pnm.h>
#include <float_image_mmorph.h>
#include <float_image_gradient.h>
#include <float_image.h>

#define BT_GAMMA (0.450)
#define BT_BIAS (0.0327)
  /* Values of {gamma} and {bias} for {sample_conv_gamma} 
    that approximate the BT.709 encoding. */
    
typedef enum 
  { tgr_op_dilate, 
    tgr_op_gradSqr, 
    tgr_op_gradSqrRel
  } tgr_op_t;
  /* Type of output image desired. */

typedef struct options_t
  { tgr_op_t op; /* Type of output image desired. */
  } options_t;

options_t *tgr_parse_options(int argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

int main(int argc, char **argv);

float_image_t *tgr_make_test_image(int NX, int NY, int m, options_t *o);
  /* Creates a test image using various
    {float_image_gradient.h} tools. */

void tgr_show_mask_statistics(float_image_t *msk, int ic);
  /* Prints statistical proeprties of channel {ic} of {msk}. */

void tgr_write_image(char *name, float_image_t *A);
  /* Writes the image {A} to a file called "{name}.{ext}",
    in the PPM/PGM format, without gamma correction. */

int main(int argc, char **argv)
  {
    /* Parse the command line options: */
    options_t *o = tgr_parse_options(argc, argv);

    /* Read the input image and get its dimensions: */
    bool_t isMask = FALSE; /* Assume smooth pixel distr. */
    float_image_t *fin = float_image_read_pnm_named("-", isMask, 1.000, 0.000, TRUE, TRUE, FALSE);
    int NC, NX, NY;
    float_image_get_size(fin, &NC, &NX, &NY);
    
    /* Compute the ouptut image: */
    float_image_t *fot;
    
    /* Weight table with max value 1.0, for dilation: */
    int hwd = 4;
    double D = 70;
    int nw = 2*hwd+1;
    assert(nw == 9);
    double wtd[9] = { 1/D, 8/D, 28/D, 56/D, 70/D, 56/D, 28/D, 1/D };
    
    int c;
    
    switch(o->op)
      { case tgr_op_dilate:
          fot = float_image_mmorph_dilate(fin, hwd, wtd);
          break;

        case tgr_op_gradSqr:
          fot = float_image_new(NC, NX, NY);
          for (c = 0; c < NC; c++)
            { float_image_gradient_sqr_sobel(fin, c, fot, c); }
          break;

        case tgr_op_gradSqrRel:
          fot = float_image_gradient_sqr_relative(fin, 0.02, FALSE);
          break;

        default:
          assert(FALSE);
      }
      
    /* Find max and min values {vMax,vMin} over all channels: */
    float vMin = +INF, vMax = -INF;
    for (c = 0; c < fot->sz[0]; c++)
      { float_image_update_sample_range(fot, c, &vMin, &vMax);  }
    fprintf(stderr, "output image range = [ %+8.4f _ %+8.4f ]\n", vMin, vMax);
    
    /* Rescale pixels to {[0_1]}: */
    for (c = 0; c < fot->sz[0]; c++)
      { float_image_rescale_samples(fot, c, 0.0, vMax, 0.0, 1.0); }
      
    /* Write output image: */
    float_image_write_pnm_named("-", fot, isMask, 1.000, 0.000, TRUE, TRUE, FALSE);
    
    float_image_free(fin); fin = NULL;
    float_image_free(fot); fot = NULL;
    free(o); o = NULL;

    return 0;
  }

options_t *tgr_parse_options(int argc, char **argv)
  {
    /* INITIALIZATION: */

    /* Start the command line analyzer {pp}: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);

    /* Process "-help" and "-info" options: */
    argparser_process_help_info_options(pp);

    /* Allocate the return record: */
    options_t *o = (options_t *)malloc(sizeof(options_t));

    /* PARSE POSITIONAL ARGUMENTS: */
    if (argparser_keyword_present(pp, "-dilate"))
      { o->op = tgr_op_dilate; }
    else if (argparser_keyword_present(pp, "-gradSqr"))
      { o->op = tgr_op_gradSqr; }
    else if (argparser_keyword_present(pp, "-gradSqrRel"))
      { o->op = tgr_op_gradSqrRel; }
    else
      { argparser_error(pp, "must specify some image op"); }

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }
