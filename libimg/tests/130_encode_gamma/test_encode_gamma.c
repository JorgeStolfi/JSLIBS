#define PROG_NAME "test_encode_gamma"
#define PROG_DESC "test of {float_image_gamma.h}"
#define PROG_VERS "1.0"

/* Last edited on 2017-06-30 01:08:39 by stolfilocal */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_encode_gamma_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -BT709 | \\\n" \
  "      -generic {GAMMA} {BIAS} [ -step {STEP} ] | \\\n" \
  "      -interp {U[0]} {V[0]}  ...  {U[NP-1]} {V[NP-1]} \\\n" \
  "    ] \\\n" \
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
  "  Criado em 2007-08-01 por J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_encode_gamma_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program writes a PGM test image that can be used" \
  " to check the light response function (LRF) of a particular" \
  " monitor.  The response function expresses the relationship" \
  " between pixel values, as stored in an image, to the" \
  " brightness of the corresponding pixels on the screen.\n" \
  "\n" \
  "  The program assumes that the LRF is either:\n" \
  "    * a generic power-law function with parameters {GAMMA} and {BIAS}; or\n" \
  "    * the modified power law specified by ITU-R BT.709; or\n" \
  "    * a piecewise-affine function defined by" \
  " zero or more user-given input-output pairs, {U[0]->V[0]}, {U[1]->V[1]}, ..." \
  " {U[NP-1]->V[NP-1]}.\n" \
  "\n" \
  "  These models are implemented by the functions {sample_conv_encode_BT709}, " \
  " {sample_conv_decode_BT709}, {sample_conv_gamma}, and {sample_conv_interp}" \
  " in the interface {sample_conv.h}.\n" \
  "\n" \
  "  The output image file is called \"out/test-enc.pgm\", and contains" \
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
  "  The program also writes a file \"out/test-enc.txt\" containing" \
  " a tabulation of the encoding function.  See the script" \
  " {show-sample-map} for details.\n" \
  "\n" \
  "Output for the BT.709 Curve\n" \
  "  For the BT.709 curve, there is a single test block.\n" \
  "\n" \
  "Output for a Generic Power Law\n" \
  "  For the generic power law, the image consists of a" \
  " matrix of test blocks.  The block in the center" \
  " corresponds to the given values of {GAMMA} and {BIAS} the" \
  " other blocks correspond to slightly different values" \
  " of those parameters.\n" \
  "\n" \
  "  If the best match occurs in a block" \
  " that lies to the left of the middle" \
  " column, the {GAMMA} parameter should" \
  " be reduced.  If the best match" \
  " occurs in a block that lies below" \
  " the middle row, the {BIAS} must be reduced." \
  " The parameters vay by a factor" \
  " of {STEP} between successive" \
  " rows/colums.\n" \
  "\n" \
  "Output for a Piecewise Affine Curve\n" \
  "  For a piecewise affine curve, the image consists of a single test block."

#define PROG_INFO_OPTS \
  "  -BT709\n" \
  "    This option, if present, requests the ITU-R BT.709 encoding" \
  " curve.  It is mutually exclusive with \"-generic\".\n" \
  "\n" \
  "  -generic {GAMMA} {BIAS}\n" \
  "    This option, if present, requests a generic power law" \
  " encoding as defined by {sample_conv_gamma} in {sample_conv.h}," \
  " with parameters {gamma=GAMMA} and {bias=BIAS}. The value" \
  " of {GAMMA} must be positive, and {BIAS} must be in" \
  " the range [0 _ 1).\n" \
  "\n" \
  "  -step {STEP}\n" \
  "    This optional argument is allowed only with" \
  " the \"-generic\" option.  It specifies the multiplicative" \
  " factor that is to be applied to the parameters {GAMMA} and" \
  " {BIAS} between successive rows or columns of the image.  The default" \
  " is \"-step 1.10\" (a 10 percent increase).\n" \
  "\n" \
  "  -interp {U[0]} {V[0]}  ...  {U[NP-1]} {V[NP-1]}\n" \
  "    This option, if present, requests a piecewise-affine RTF" \
  " as defined by {sample_conv_interp} in" \
  " {sample_conv.h}.  The function will map 0 to 0," \
  " {U[i]} to {V[i]} for {i} in {0..NP-1}, and 1 to 1; and will" \
  " be anti-symmetric around 0.  In particular, f {NP} is zero" \
  " (i.e. no pairs are given after \"-interp\"), the function" \
  " is the identity.  The values {U[i],V[i]} must" \
  " be in the range (0 _ 1), and each coordinate must" \
  " be strictly increasing with {i}.\n"                \
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

#include <argparser.h>
#include <vec.h>
#include <affirm.h>
#include <jsfile.h>
#include <bool.h>
#include <sample_conv.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <float_image_to_uint16_image.h>
#include <float_image.h>

typedef struct options_t
  { /* BT709 function params: */
    bool_t BT709;     /* TRUE uses BT.709 curve. */
    /* Generic function params: */
    bool_t generic;   /* TRUE uses a generic power law with parameters {gamma,bias}. */
    double gamma;     /* Exponent of generic power law. */
    double bias;      /* Bias value for generic power law. */
    double step;
    /* Interpolated function params: */
    bool_t interp;    /* TRUE uses a piecewise affine map with nodes {(u[i],v[i])}. */
    double_vec_t U;   /* Input data values. */
    double_vec_t V;   /* Output data values. */
    /* Texture params: */
    bool_t dots;      /* TRUE allows dot-in-four dither. */
    bool_t vertical;  /* TRUE uses vertical lines. */
  } options_t;

options_t *fitg_parse_options(int argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} rec%ord. */

void fitg_parse_pair_list(argparser_t *pp, double_vec_t *U, double_vec_t *V);
  /* Parses from the command line parser {pp} a list of {n}
    consecutive number pairs, and stores them into the vectors {U} and
    {V}. The vectors are trimmed to {np} elements each. Checks whether
    each list is strictly increasing and all numbers are in {(0_1)}. */

int main(int argc, char **argv);

void fitg_dump_gamma_table(char *fname, char *suff, int nv, options_t *o);
  /* Writes to file "{name}{suff}" a set of {2*nv+1} quadruples
      "{vRaw} {vCor} {vInf} {mRaw} {mCor} {mInv}",
    one per line, where {vRaw} is a raw sample value, {vCor} is the
    encoded version of {vRaw}, {vInv} is the
    value decoded from {vCor}, and {mRaw,mCor,mInv}
    are the base-10 logarithms of {|vRaw|,|vCor|,|vInv|}, shifted so that
    they are always positive, times their signs.
    The values {vRaw} span {[-1 _ +1]} with {2*nv} equal steps. */

float_image_t *fitg_make_gamma_test_image
  ( int mbx,  /* Number of extra block per row, on each side. */
    int mby,  /* Number of extra block per column, on each side. */
    int nx,
    int ny,
    int mrg,
    options_t *o
  );
  /* Creates a test image for testing {sample_conv_gamma},
    {sample_conv_encode_BT709}, or {sample_conv_interp}. The image
    consists of an array of blocks with {NBX = 2*mbx+1} columns and {NBY
    = 2*mby+1} rows, where the central block tests the parameters
    specified by {o}, and the other blocks (if any) test of small variations
    thereof. The parameters change by a multiplicative factor {o->step}
    between succesive rows or columns.

    Each block has {nx} by {ny} pixels, including a solid black margin
    {mrg} pixels wide, and is painted with
    {fitg_paint_gamma_test_block}. */

void fitg_paint_gamma_test_block
  ( float_image_t *A,
    int ibx,
    int iby,
    int dx,
    int dy,
    int nx,
    int ny,
    int mrg,
    options_t *o
  );
  /* Paints into {A} a test block for testing {sample_conv_gamma},
    {sample_conv_encode_BT709}., or {sample_conv_interp} .

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

void fitg_write_image(char *name, char *suff, float_image_t *A);
  /* Writes the image {A} to a file called "{name}{suff}",
    in the PPM/PGM format, without gamma correction. */

int main(int argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = fitg_parse_options(argc, argv);

    char *oname = "out/test-enc"; /* Prefix of output file names. */

    /* Test block geometry: */
    int nx;   /* Block width in pixels. */
    int ny;   /* Block height in pixels */
    int mrg;  /* Width of black frame. */
    int mbx;   /* Number of extra {gamma} values to try on each side of {o->gamma}. */
    int mby;   /* Number of extra {bias} values to try on each side of {o->bias}. */
    if (o->BT709)
      { /* No parameter variation possible: */
        mbx = 0; mby = 0;
        /* We can afford a large block: */
        nx = 92; ny = 400; mrg = 1;
      }
    else if (o->generic)
      { /* Gamma variation possible: */
        mbx = 2; 
        /* Bias variation possible if nonzero: */
        mby = (o->bias == 0 ? 0 : 1);
        /* Pick a modest block size: */
        nx = 80; ny = 200; mrg = 1;
      }
    else if (o->interp)
      { /* No parameter variation possible: */
        mbx = 0; mby = 0;
        /* We can afford a large block: */
        nx = 92; ny = 400; mrg = 1;
      }
    else
      { demand(FALSE, "no image kind?"); }

    fitg_dump_gamma_table(oname, ".txt", 100, o);

    /* Generate the test image {TG}: */
    float_image_t *TG = fitg_make_gamma_test_image(mbx, mby, nx, ny, mrg, o);

    /* Write it out: */
    fitg_write_image(oname, ".pgm", TG);

    /* Cleanup: */
    float_image_free(TG); TG = NULL;
    free(o); o = NULL;

    return 0;
  }

void fitg_dump_gamma_table(char *name, char *suff, int nv, options_t *o)
  { char *fname = NULL;
    char *fname = jsprintf("%s%s", name, suff);
    FILE *wr = open_write(fname, TRUE);
    float eRaw = (float)(((double)1)/((double)nv));  /* Min positive value of {vRaw} */
    float eCor;  /* Min positive value of {vCor} */
    if (o->BT709)
      { eCor = sample_conv_encode_BT709(eRaw); }
    else if (o->interp)
      { eCor = sample_conv_interp(eRaw, o->U.ne, o->U.e, o->V.e); }
    else if (o->generic)
      { eCor = sample_conv_gamma(eRaw, o->gamma, o->bias); }
    else
      { affirm(FALSE, "bug"); }
    int iv;
    for (iv = -nv; iv <= nv; iv++)
      { float vRaw = (float)(((double)iv)/((double)nv));
        float vCor, vInv;
        if (o->BT709)
          { vCor = sample_conv_encode_BT709(vRaw);
            vInv = sample_conv_decode_BT709(vCor);
          }
        else if (o->interp)
          { vCor = sample_conv_interp(vRaw, o->U.ne, o->U.e, o->V.e);
            vInv = sample_conv_interp(vCor, o->V.ne, o->V.e, o->U.e);
          }
        else if (o->generic)
          { vCor = sample_conv_gamma(vRaw, o->gamma, o->bias);
            vInv = sample_conv_gamma(vCor, 1/o->gamma, o->bias);
          }
        else
          { affirm(FALSE, "bug"); }
        fprintf(wr, " %+7.4f %+7.4f %+7.4f", vRaw, vCor, vInv);
        fprintf(wr, "  ");

        auto double vlog(float v, float e);
        double vlog(float v, float e)
          { return log(fabs(v)/e)/M_LN10 * (v < 0 ? -1 : +1); }

        if (iv != 0)
          {
            double mRaw = vlog(vRaw, eRaw);
            double mCor = vlog(vCor, eCor);
            double mInv = vlog(vInv, eRaw);
            fprintf(wr, "  %+9.4f %+9.4f %+9.4f", mRaw, mCor, mInv);
          }
        fprintf(wr, "\n");
      }
    fclose(wr);
    free(fname);
  }

float_image_t *fitg_make_gamma_test_image
  ( int mbx,
    int mby,
    int nx,
    int ny,
    int mrg,
    options_t *o
  )
  {
    int NBX = 2*mbx + 1;  /* Blocks per row. */
    int NBY = 2*mby + 1;  /* Blocks per column. */

    int NX = NBX*nx;     /* Image width in pixels. */
    int NY = NBY*ny;     /* Image height in pixels. */

    /* Create image, fill it with  black: */
    float_image_t *A = float_image_new(1, NX, NY);
    float_image_fill(A, 0.0);

    /* Paint the blocks: */
    int ibx, iby;
    for (ibx = 0; ibx < NBX; ibx++)
      { int dx = ibx*nx;
        for (iby = 0; iby < NBY; iby++)
          { int dy = iby*ny;
            fitg_paint_gamma_test_block(A, ibx-mbx, iby-mby, dx, dy, nx, ny, mrg, o);
          }
      }
    return A;
  }

void fitg_paint_gamma_test_block
  ( float_image_t *A,
    int ibx,
    int iby,
    int dx,
    int dy,
    int nx,
    int ny,
    int mrg,
    options_t *o
  )
  { /* Internal dimensions of block: */
    int txSz = (nx - 2*mrg)/3;          /* Width of textured bands. */
    int gxSz = (nx - 2*mrg - 2*txSz); /* Width of gray band. */
    int gxIni = mrg + txSz;             /* Initial {ix} of gray band. */
    int gxFin = gxIni + gxSz;           /* Final {ix} of gray band. */

    double gBlock = 0, bBlock = 0; /* For generic power law. */
    if (o->generic) 
      { gBlock = o->gamma * pow(o->step, ibx);  /* Block's {gamma}. */
        bBlock = o->bias * pow(o->step, iby);  /* Block's {bias}. */
        fprintf(stderr, "gamma = %6.4f 1/gamma = %6.4f  bias = %6.4f\n", gBlock, 1/gBlock, bBlock);
      }
                
    /* Scan block pixels and paint accordingly: */
    int ix, iy;
    for (iy = 0; iy < ny; iy++)
      {
        /* Check whether pixel row is part of the frame: */
        bool_t y_frame = ((iy < mrg) || (iy >= ny - mrg));

        /* Compute desired brightness {lum} for this line: */
        double lum = (y_frame ? 0 : ((double)iy - mrg + 0.5)/((double)ny - 2*mrg));

        /* Choose texture {tx} and values {loLum,hiLum} of light and dark texture pixels: */
        unsigned int tx; /* Bit {2*t1 + t0} tells pixel value: {loLum} (0) or {hiLum} (1). */
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

        for (ix = 0; ix < nx; ix++)
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
                int bx = (ix < gxIni ? gxIni - 1 - ix : ix - gxFin - 1);
                int by = iy - mrg;
                /* Determine the texture bit coordinates {t0,t1}: */
                int t0 = (o->vertical ? by : bx) % 2; /* texture coordinate. */
                int t1 = (o->vertical ? bx : by) % 2; /* Second texture coordinate. */
                /* Determine the texture bit selection mask {txsel}: */
                unsigned int txsel = 1 << (2*t1 + t0);
                /* Now pick the color according to the texture: */
                vRaw = (float)((tx & txsel) == 0 ? loLum : hiLum);
              }
            float vCor;
            if (o->BT709)
              { vCor = sample_conv_encode_BT709(vRaw); }
            else if (o->interp)
              { vCor = sample_conv_interp(vRaw, o->U.ne, o->U.e, o->V.e); }
            else if (o->generic)
              { vCor = sample_conv_gamma(vRaw, gBlock, bBlock); }
            else
              { affirm(FALSE, "bug"); }
            float_image_set_sample(A, 0, dx + ix, dy + iy, vCor);
          }
      }
  }

void fitg_write_image(char *name, char *suff, float_image_t *A)
  { char *fname = NULL;
    char *fname = jsprintf("%s%s", name, suff);
    FILE *wr = open_write(fname, TRUE);
    int chns = (int)A->sz[0];
    bool_t yup = TRUE, verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    uint16_image_t *pimg = float_image_to_uint16_image(A, isMask, chns, NULL, NULL, NULL, 255, yup, verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pimg, forceplain, verbose);
    uint16_image_free(pimg);
    fclose(wr);
    free(fname);
  }

options_t *fitg_parse_options(int argc, char **argv)
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

    /* Set all LRF function parameters, just in case. */
    o->BT709 = o->generic = o->interp = FALSE; /* Default */
    o->gamma = 1.000; o->bias = 0.000; o->step = 1.100; /* Just in case. */
    o->U = double_vec_new(0);
    o->V = double_vec_new(0);

    /* Parse the LRF function parameters: */
    if (argparser_keyword_present(pp, "-BT709"))
      { o->BT709 = TRUE; }
    else if (argparser_keyword_present(pp, "-interp"))
      { o->interp = TRUE;
        fitg_parse_pair_list(pp, &(o->U), &(o->V));
        affirm(o->U.ne == o->V.ne, "bug parse pair list"); 
      }
    else if (argparser_keyword_present(pp, "-generic"))
      { o->generic = TRUE;

        o->gamma = argparser_get_next_double(pp, 1.0e-100, 1000.0);
        if (o->gamma <= 0) { argparser_error(pp, "{gamma} must be positive"); }

        o->bias = argparser_get_next_double(pp, 0.0, 0.999999);
        if (o->bias < 0) { argparser_error(pp, "{bias} canot be negative"); }
        if (o->bias >= 1) { argparser_error(pp, "{bias} must be less than 1"); }


        if (argparser_keyword_present(pp, "-step"))
          { o->step = argparser_get_next_double(pp, 0.1, 10.0); }
        else
          { o->step = 1.10; }
      }
    else
      { argparser_error(pp, "must specify either \"-BT709\" or \"-generic\""); }

    o->dots = argparser_keyword_present(pp, "-dots");

    o->vertical = argparser_keyword_present(pp, "-vertical");

    /* PARSE POSITIONAL ARGUMENTS: */

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }

void fitg_parse_pair_list(argparser_t *pp, double_vec_t *U, double_vec_t *V)
  { int np = 0;
    int j;
    while(argparser_next_is_number(pp))
      { int ip = np; np++;
        for (j = 0; j < 2; j++)
          { double_vec_t *w = ( j == 0 ? U : V );
            double_vec_expand(w, ip);
            w->e[ip] = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
            if ((w->e[ip] <= 0) || (w->e[ip] >= 1))
              { argparser_error(pp, "LRF data values must be in (0_1)"); }
            if ((ip > 0) && (w->e[ip] <= w->e[ip-1]))
              { argparser_error(pp, "LRF data values must be increasing"); }
          }
      }
    double_vec_trim(U, np);
    double_vec_trim(V, np);
  }
