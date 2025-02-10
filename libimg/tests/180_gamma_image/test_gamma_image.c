#define PROG_NAME "test_gamma_image"
#define PROG_DESC "test of {float_image_gamma.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-01-30 08:05:49 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_gamma_image_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -kind { BT709 | sRGB | generic }\\\n" \
  "    [ -expo {EXPO} -bias {BIAS} -step {STEP} ] \\\n" \
  "    [ -dots ] [ -vertical ] \\\n" \
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
  "  Created on 2007-08-01 by J. Stolfi, IC-UNICAMP.\n" \
  "MODIFICATION HISTORY\n" \
  "  By J.Stolfi unless noted otherwise.\n" \
  "  2007-08-01 Created.\n" \
  "  2024-12-18 Added sRGB conversion.\n" \
  "  2024-12-20 Split the curve plots off to another prog.\n" \
  "  2024-12-21 Removed the \"interp\" option.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_gamma_image_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program writes a PGM test image that can be used" \
  " to check the light response function (LRF) of a particular" \
  " monitor.  The response function expresses the relationship" \
  " between pixel values, as stored in an image, to the" \
  " brightness of the corresponding pixels on the screen.\n" \
  "\n" \
  "\n" \
  "  The conversion type is specified by the {KIND} parameter, which may be:\n" \
  "\n" \
  "    * \"generic\" a generic modified power-law function with" \
  "       parameters {EXPO} and {BIAS}; or\n" \
  "\n" \
  "    * \"BT707\" the modified power law specified" \
  "      by the ITU-R BT.709 standard; or\n" \
  "\n" \
  "    * \"sRGB\" the modified power law" \
  "      specified by the IEC sRGB standard; or\n" \
  "\n" \
  "  These models are implemented by the" \
  " functions {sample_conv_BT709_encode}, " \
  " {sample_conv_BT709_decode}, {sample_conv_sRGB_encode}, " \
  " {sample_conv_sRGB_decode}, and {sample_conv_gamma}.\n" \
  "\n" \
  "  The output image file is called \"out/test-{KIND}.pgm\", and contains" \
  " one or more /test blocks/. Each test block consists of three" \
  " vertical bands enclosed in a thin black frame.  The" \
  " central band is a gradient of solid gray colors.  The two" \
  " side bands have a varying light/dark texture. To use the image," \
  " put it up on the monitor with \"display -gamma 1.0\".  If the" \
  " LTF specified to the program matches that of the monitor, the three bands" \
  " should have the same average brightness;" \
  " so that, if the image is viewed from such a distance that" \
  " the texture bands becomes smooth gradients, the three bands" \
  " should seem to be a single gradient of gray tones." \
  "\n" \
  "  When {KIND} is \"BT.709\" or \"sRGB\", the\n" \
  " output image has a single test block.\n" \
  "\n" \
  "  When {KIND} is \"generic\", the image consists of a" \
  " matrix of test blocks.  The block in the center" \
  " corresponds to the given values of {EXPO} and {BIAS} the" \
  " other blocks correspond to slightly different values" \
  " of those parameters.\n" \
  "\n" \
  "  If the best match occurs in a block" \
  " that lies to the left of the middle" \
  " column, the {EXPO} parameter should" \
  " be reduced.  If the best match" \
  " occurs in a block that lies below" \
  " the middle row, the {BIAS} must be reduced." \
  " The parameters vary by a factor" \
  " of {STEP} between successive" \
  " rows or colums."

#define PROG_INFO_OPTS \
  "  -kind {KIND}\n" \
  "    This mandatory option specifies the type of transfer function to" \
  " plot. See above for the {KIND} alternatives and meaning.\n" \
  "\n" \
  "  -expo {EXPO}\n" \
  "  -bias {BIAS}\n" \
  "  -step {STEP}\n" \
  "    These options should be present if and only if {KIND}" \
  " is \"generic\", and specify, the exponent {EXPO} and" \
  " bias {BIAS} of {sample_conv_gamma}, as well their" \
  " increment between image blocks. The value" \
  " of {EXPO} must be positive, and {BIAS} must be in" \
  " the range [0 _ 1).  The parameters \n" \
  "\n" \
  "  -dots\n" \
  "    This optional argument specifies that luminosity" \
  " values below 0.25 or above 0.75 should be approximated" \
  " by one-dot-in-four textures, intead of hatched textures.\n" \
  "\n" \
  "  -vertical\n" \
  "    This optional argument specifies that hatched textures" \
  " should use vertical lines, rather than horizontal lines."

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include <argparser.h>
#include <vec.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <bool.h>
#include <sample_conv.h>
#include <sample_conv_gamma.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <float_image_to_uint16_image.h>
#include <float_image.h>

typedef struct options_t
  { /* Encoding/decoding params: */
    char *kind;       /* "BT709", "sRGB", "generic", or "interp". */
    int32_t maxval;   /* Max integer sample value. */
    /* Encoding/decoding params for generic gamma: */
    double expo;      /* Exponent of generic power law. */
    double bias;      /* Bias value for generic power law. */
    double step;      /* INcrement of {expo} and {bias} between steps. */
    /* Texture params: */
    bool_t dots;      /* TRUE allows dot-in-four dither. */
    bool_t vertical;  /* TRUE uses vertical lines. */
  } options_t;

options_t *tgim_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

int32_t main(int32_t argc, char **argv);

float_image_t *tgim_make_gamma_test_image
  ( int32_t mbx,  /* Number of extra block per row, on each side. */
    int32_t mby,  /* Number of extra block per column, on each side. */
    int32_t nx,
    int32_t ny,
    int32_t mrg,
    options_t *o
  );
  /* Creates a test image for testing {sample_conv_gamma},
    {sample_conv_BT709_encode}, or {sample_conv_sRGB_encode}. The image
    consists of an array of blocks with {NBX = 2*mbx+1} columns and {NBY
    = 2*mby+1} rows, where the central block tests the parameters
    specified by {o}, and the other blocks (if any) test of small variations
    thereof. The parameters change by a multiplicative factor {o->step}
    between succesive rows or columns.

    Each block has {nx} by {ny} pixels, including a solid black margin
    {mrg} pixels wide, and is painted with
    {tgim_paint_gamma_test_block}. */

void tgim_paint_gamma_test_block
  ( float_image_t *A,
    int32_t ibx,
    int32_t iby,
    int32_t dx,
    int32_t dy,
    int32_t nx,
    int32_t ny,
    int32_t mrg,
    options_t *o
  );
  /* Paints into {A} a test block for testing {sample_conv_gamma},
    {sample_conv_BT709_encode}, {sample_conv_sRGB_encode}, or {sample_conv_interp} .

    The block is in column {ibx} and row {iby} of the image, where
    {ibx==iby==0} for the central block. It has {nx} by {ny} pixels, and
    its lowest pixel is pixel {(dx,dy)} of {A}. The block consists of
    a central band with a vertical gradient of solid gray colors,
    flanked by bands of a texture whose average brightnes is supposed
    to match that of the central band. The block includes a black
    frame of with {mrg}. Sample values in the whole block are then
    encode as pecified by {o}.

    If {o->dots} is false, the texture consists of alternating light
    and dark single-pixel lines, for all brightness values. If
    {o->dots} is true, uses one-dot-in-four dither for brightness
    values below {0.25} or above {0.75}, and lines otherwise.

    If {o->vertical} is true, the line-based textures consist of
    vertical lines, othwewise they consist of horizontal ones. */

void tgim_write_image(char *name, char *suff, float_image_t *A);
  /* Writes the image {A} to a file called "{name}{suff}",
    in the PPM/PGM format, without gamma correction. */

int32_t main(int32_t argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = tgim_parse_options(argc, argv);

    char *oname = jsprintf("out/test-%s", o->kind); /* Prefix of output file names. */

    /* Test block geometry: */
    int32_t nx;    /* Block width in pixels. */
    int32_t ny;    /* Block height in pixels */
    int32_t mrg;   /* Width of black frame. */
    int32_t mbx;   /* Number of extra {expo} values to try on each side of {o->expo}. */
    int32_t mby;   /* Number of extra {bias} values to try on each side of {o->bias}. */
    if (((strcmp(o->kind, "BT709") == 0)) || ((strcmp(o->kind, "sRGB") == 0)))
      { /* No parameter variation possible: */
        mbx = 0; mby = 0;
        /* We can afford a large block: */
        nx = 92; ny = 400; mrg = 1;
      }
    else if ((strcmp(o->kind, "generic") == 0))
      { /* Gamma variation possible: */
        mbx = 2; 
        /* Bias variation possible if nonzero: */
        mby = (o->bias == 0 ? 0 : 1);
        /* Pick a modest block size: */
        nx = 80; ny = 200; mrg = 1;
      }
    else
      { demand(FALSE, "no image kind?"); }

    /* Generate the test image {TG}: */
    float_image_t *TG = tgim_make_gamma_test_image(mbx, mby, nx, ny, mrg, o);

    /* Write it out: */
    tgim_write_image(oname, ".pgm", TG);

    /* Cleanup: */
    float_image_free(TG); TG = NULL;
    free(o); o = NULL;

    return 0;
  }

float_image_t *tgim_make_gamma_test_image
  ( int32_t mbx,
    int32_t mby,
    int32_t nx,
    int32_t ny,
    int32_t mrg,
    options_t *o
  )
  {
    int32_t NBX = 2*mbx + 1;  /* Blocks per row. */
    int32_t NBY = 2*mby + 1;  /* Blocks per column. */

    int32_t NX = NBX*nx;     /* Image width in pixels. */
    int32_t NY = NBY*ny;     /* Image height in pixels. */

    /* Create image, fill it with  black: */
    float_image_t *A = float_image_new(1, NX, NY);
    float_image_fill(A, 0.0);

    /* Paint the blocks: */
    for (int32_t ibx = 0; ibx < NBX; ibx++)
      { int32_t dx = ibx*nx;
        for (int32_t iby = 0; iby < NBY; iby++)
          { int32_t dy = iby*ny;
            tgim_paint_gamma_test_block(A, ibx-mbx, iby-mby, dx, dy, nx, ny, mrg, o);
          }
      }
    return A;
  }

void tgim_paint_gamma_test_block
  ( float_image_t *A,
    int32_t ibx,
    int32_t iby,
    int32_t dx,
    int32_t dy,
    int32_t nx,
    int32_t ny,
    int32_t mrg,
    options_t *o
  )
  { /* Internal dimensions of block: */
    int32_t txSz = (nx - 2*mrg)/3;          /* Width of textured bands. */
    int32_t gxSz = (nx - 2*mrg - 2*txSz); /* Width of gray band. */
    int32_t gxIni = mrg + txSz;             /* Initial {ix} of gray band. */
    int32_t gxFin = gxIni + gxSz;           /* Final {ix} of gray band. */

    double gBlock = 0, bBlock = 0; /* For generic power law. */
    if ((strcmp(o->kind, "generic") == 0)) 
      { gBlock = o->expo * pow(o->step, ibx);  /* Block's {expo}. */
        bBlock = o->bias * pow(o->step, iby);  /* Block's {bias}. */
        fprintf(stderr, "expo = %6.4f 1/expo = %6.4f  bias = %6.4f\n", gBlock, 1/gBlock, bBlock);
      }
                
    /* Scan block pixels and paint accordingly: */
    for (int32_t iy = 0; iy < ny; iy++)
      {
        /* Check whether pixel row is part of the frame: */
        bool_t y_frame = ((iy < mrg) || (iy >= ny - mrg));

        /* Compute desired brightness {lum} for this line: */
        double lum = (y_frame ? 0 : ((double)iy - mrg + 0.5)/((double)ny - 2*mrg));

        /* Choose texture {tx} and values {loLum,hiLum} of light and dark texture pixels: */
        uint32_t tx; /* Bit {2*t1 + t0} tells pixel value: {loLum} (0) or {hiLum} (1). */
        double hiLum, loLum;
        if (o->dots && (lum <= 0.25))
          { /* Light dots on black field: */
            loLum = 0.0; hiLum = 4*lum; tx = 1;       /* 2_0001 = lower left pixel is light. */
          }
        else if (o->dots && (lum >= 0.75))
          { /* Dark dots on white field: */
            loLum = 4*lum - 3.0; hiLum = 1.0; tx = 7; /* 2_0111 = upper right pixel is dark. */
          }
        else if (lum <= 0.5)
          { /* Light stripes on black field: */
            loLum = 0.0; hiLum = 2*lum; tx = 3;       /* 2_0011 = lower line is light. */
          }
        else
          { /* White stripes on light field: */
            loLum = 2*lum - 1.0 ; hiLum = 1.0; tx = 3; /* 2_0011 = lower line is light. */
          }

        for (int32_t ix = 0; ix < nx; ix++)
          { /* Check whether pixel column is part of the frame: */
            bool_t x_frame = ((ix < mrg) || (ix >= nx - mrg));
            /* Determine pixel value {vRaw}: */
            float vRaw;
            if (x_frame || y_frame)
              { /* Pixel is in frame: */
                vRaw = 0.0;
              }
            else if ((ix >= gxIni) && (ix <= gxFin))
              { /* Pixel is in middle band: */
                vRaw = (float)lum;
              }
            else
              { /* Pixel is in textured bands. */
                /* Determine the distances {bx,by} from the band's low inner corner: */
                int32_t bx = (ix < gxIni ? gxIni - 1 - ix : ix - gxFin - 1);
                int32_t by = iy - mrg;
                /* Determine the texture bit coordinates {t0,t1}: */
                int32_t t0 = (o->vertical ? by : bx) % 2; /* texture coordinate. */
                int32_t t1 = (o->vertical ? bx : by) % 2; /* Second texture coordinate. */
                /* Determine the texture bit selection mask {txsel}: */
                uint32_t txsel = 1 << (2*t1 + t0);
                /* Now pick the color according to the texture: */
                vRaw = (float)((tx & txsel) == 0 ? loLum : hiLum);
              }
            float vCor;
            if ((strcmp(o->kind, "BT709") == 0))
              { vCor = sample_conv_BT709_encode(vRaw); }
            else if ((strcmp(o->kind, "sRGB") == 0))
              { vCor = sample_conv_sRGB_encode(vRaw); }
            else if ((strcmp(o->kind, "generic") == 0))
              { vCor = sample_conv_gamma(vRaw, gBlock, bBlock); }
            else
              { affirm(FALSE, "bug"); }
            float_image_set_sample(A, 0, dx + ix, dy + iy, vCor);
          }
      }
  }

void tgim_write_image(char *name, char *suff, float_image_t *A)
  { char *fname = jsprintf("%s%s", name, suff);
    FILE *wr = open_write(fname, TRUE);
    int32_t chns = (int32_t)A->sz[0];
    bool_t yUp = TRUE, verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    uint16_image_t *pimg = float_image_to_uint16_image(A, isMask, chns, NULL, NULL, NULL, 255, yUp, verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pimg, forceplain, verbose);
    uint16_image_free(pimg);
    fclose(wr);
    free(fname);
  }

options_t *tgim_parse_options(int32_t argc, char **argv)
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

    /* PARSE KEYWORD ARGUMENTS: */

    argparser_get_keyword(pp, "-kind");
    o->kind = argparser_get_next_non_keyword(pp);
    
    /* Defaults even if not used: */
    o->expo = 1.000; o->bias = 0.000; o->step = 1.100;
    
    if (strcmp(o->kind, "generic") == 0)
      { argparser_get_keyword(pp, "-expo");
        o->expo = argparser_get_next_double(pp, 1.0e-100, 1000.0);
        
        argparser_get_keyword(pp, "-bias");
        o->bias = argparser_get_next_double(pp, 0.0, 0.999999);

        argparser_get_keyword(pp, "-step");
        o->step = argparser_get_next_double(pp, 1.001, 10.0);

      }
    else if ((strcmp(o->kind, "BT709") != 0) && (strcmp(o->kind, "sRGB") != 0))
      { argparser_error(pp, "invalid \"-kind\" argument"); }

    /* PARSE POSITIONAL ARGUMENTS: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */
    argparser_finish(pp);

    return o;
  }
