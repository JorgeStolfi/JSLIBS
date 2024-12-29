#define PROG_NAME "test_encode_gamma"
#define PROG_DESC "test of gamma enc/dec functions in {sample_conv_gamma.h} etc."
#define PROG_VERS "1.0"

/* Last edited on 2024-12-25 12:48:20 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_encode_gamma_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -kind { BT709 | sRGB | generic |interp }\\\n" \
  "    -maxval {MAXVAL} \\\n" \
  "    [ -expo {EXPO} -bias {BIAS} ] \\\n" \
  "    [ -points {U[0]} {V[0]}  ...  {U[NP-1]} {V[NP-1]} ] \\\n" \
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
  "  Created on 2007-08-01 por J. Stolfi, IC-UNICAMP.\n" \
  "MODIFICATION HISTORY\n" \
  "  By J.Stolfi unless noted otherwise.\n" \
  "  2007-08-01 Created.\n" \
  "  2024-12-18 Added sRGB conversion.\n" \
  "  2024-12-20 Split the gamma gauge image off to another prog.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_encode_gamma_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program plots some non-linear light mapping functions" \
  " implemented in the {libimg} library.  These functions come in pairs, a" \
  " decoding function {D} an an encoding function {E} that express the relationship" \
  " between an integer sample value {z}, as stored in an image in file" \
  " or memory, to the intensity {v = D(z/M)} of the corresponding color channel" \
  " on the screen, also in {[0 _1]}; where {M} is the maximum sample value" \
  " allowed in the image.  The encoding function {E} is supposed to be" \
  " the inverse of {D}.  Both are extended to negative values by" \
  " the identities {D(-u) = -D(u)} and {E(-u) = -E(u)}." \
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
  "    * \"interp\", a piecewise-affine function defined by" \
  "      zero or more user-given input-output pairs, {U[0]->V[0]}, {U[1]->V[1]}, ..." \
  "      {U[NP-1]->V[NP-1]}.\n" \
  "\n" \
  "  These conversions are implemented by the functions {sample_conv_gamma}," \
  " {sample_conv_BT709_encode}, {sample_conv_BT709_decode}, {sample_conv_sRGB_encode}," \
  " {sample_conv_sRGB_decode}, and {sample_conv_interp}.  \n" \
  "\n" \
  "  The program writes a file \"out/test-{KIND}.txt\" containing" \
  " a tabulation of the encoding function.  See the script" \
  " {show_sample_map.sh} for details."

#define PROG_INFO_OPTS \
  "  -maxval {MAXVAL}\n" \
  "    This mandatory option specifies the range or integer sample values" \
  " to assume, namely {-MAXVAL..+MAXVAL}.\n" \
  "\n" \
  "  -kind {KIND}\n" \
  "    This mandatory option specifies the type of transfer function to" \
  " plot. See above for the {KIND} alternatives and meaning.\n" \
  "\n" \
  "  -expo {EXPO}\n" \
  "  -bias {BIAS}\n" \
  "    These options should be present if and only if {KIND} is \"generic\", and specify," \
  " the exponent {EXPO} and bias {BIAS} of {sample_conv_gamma}. The value" \
  " of {EXPO} must be positive, and {BIAS} must be in" \
  " the range [0 _ 1)..\n" \
  "\n" \
  "  -points {U[0]} {V[0]}  ...  {UP[N-1]} {V[NP-1]}\n" \
  "    This option should be present if and only if {KIND} is \"interp\".  It specifies" \
  " the parameters {U} and {V} for {sample_conv_interp}. The values {U[i],V[i]} must" \
  " be in the range (0 _ 1), and each coordinate must" \
  " be strictly increasing with {i}.  The function will map 0 to 0," \
  " {U[i]} to {V[i]} for {i} in {0..NP-1}, and 1 to 1; and will" \
  " be anti-symmetric around 0.  In particular, if {NP} is zero" \
  " (i.e. no pairs are given after \"-points\"), the function" \
  " is the identity. " 

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
    double expo;     /* Exponent of generic power law. */
    double bias;      /* Bias value for generic power law. */
    /* Params for interpolated transfer function: */
    double_vec_t U;   /* Input data values. */
    double_vec_t V;   /* Output data values. */
  } options_t;

options_t *teng_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

void teng_parse_pair_list(argparser_t *pp, double_vec_t *U, double_vec_t *V);
  /* Parses from the command line parser {pp} a list of {n}
    consecutive number pairs, and stores them into the vectors {U} and
    {V}. The vectors are trimmed to {np} elements each. Checks whether
    each list is strictly increasing and all numbers are in {(0_1)}. */

int32_t main(int32_t argc, char **argv);

void teng_dump_gamma_table(char *fname, options_t *o);
  /* Writes to file "{name}{suff}" a set of {2*nv+1} quadruples
      "{vRaw} {vCor} {vInf} {mRaw} {mCor} {mInv}",
    one per line, where {vRaw} is a raw sample value, {vCor} is the
    encoded version of {vRaw}, {vInv} is the
    value decoded from {vCor}, and {mRaw,mCor,mInv}
    are the base-10 logarithms of {|vRaw|,|vCor|,|vInv|}, shifted so that
    they are always positive, times their signs.
    The values {vRaw} span {[-1 _ +1]} with {2*nv} equal steps. */

int32_t main(int32_t argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = teng_parse_options(argc, argv);

    char *oname = jsprintf("out/test-%s", o->kind); /* Prefix of output file names. */

    teng_dump_gamma_table(oname, o);

    /* Cleanup: */
    free(o);
    free(oname);
    return 0;
  }

void teng_dump_gamma_table(char *name, options_t *o)
  { 
    auto double vlog(float v, double e);
      /* Converts {v} to a "signed log" scale assuming that {e} is 
        the smallest positive value.  Returns {NAN} if {v} is zero. */
    
    auto int32_t quant(float v, int32_t maxval);
      /* The value {v*maxval} to the nearest integer, symmetrically about zero. */

    int32_t maxval = o->maxval;
    /* {iRaw}: integer sample value in {0..maxval}, as in file. */ 
    /* {vRaw}: fractional nonlinear sample value converted from {-maxval..+maxval} to {[-1 _ +1]}. */
    /* {vLin}: linearized sample value, {D(vRaw)}. */
    /* {vEnc}: non-linear re-encoded sample value {round(E(vLin)}, rounded to multiple of {1/maxval}. */

    double eRaw = (((double)1)/((double)maxval));  /* Min positive value of {vRaw,vEnc} */
    double eLin;  /* Min positive value of {vLin} */
    if ((strcmp(o->kind, "BT709") == 0))
      { eLin = sample_conv_BT709_decode((float)eRaw); }
    else if ((strcmp(o->kind, "sRGB") == 0))
      { eLin = sample_conv_sRGB_decode((float)eRaw); }
    else if ((strcmp(o->kind, "interp") == 0))
      { assert(o->U.ne == o->V.ne);
        eLin = sample_conv_interp((float)eRaw, (int32_t)o->V.ne, o->V.e, o->U.e);
      }
    else if ((strcmp(o->kind, "generic") == 0))
      { eLin = sample_conv_gamma((float)eRaw, 1/o->expo, o->bias); }
    else
      { affirm(FALSE, "bug"); }

    char *fname = jsprintf("%s%s", name, ".txt");
    FILE *wr = open_write(fname, TRUE);
    for (int32_t iRaw = -maxval; iRaw <= maxval; iRaw++)
      { float vRaw = (float)(((double)iRaw)/((double)maxval));
        assert(iRaw == quant(vRaw, maxval)); /* Check on {quant}. */
        float vLin, vEnc;
        if ((strcmp(o->kind, "BT709") == 0))
          { vLin = sample_conv_BT709_decode(vRaw);
            vEnc = sample_conv_BT709_encode(vLin);
            
          }
        else if ((strcmp(o->kind, "sRGB") == 0))
          { vLin = sample_conv_sRGB_decode(vRaw);
            vEnc = sample_conv_sRGB_encode(vLin);
            
          }
        else if ((strcmp(o->kind, "interp") == 0))
          { vLin = sample_conv_interp(vRaw, (int32_t)o->V.ne, o->V.e, o->U.e);
            vEnc = sample_conv_interp(vLin, (int32_t)o->U.ne, o->U.e, o->V.e);
            
          }
        else if ((strcmp(o->kind, "generic") == 0))
          { vLin = sample_conv_gamma(vRaw, 1/o->expo, o->bias);
            vEnc = sample_conv_gamma(vLin, o->expo, o->bias);
          }
        else
          { affirm(FALSE, "bug"); }
        
        double zErr = (vEnc - vRaw)*(double)maxval; /* Error in quantized value. */
          
        int32_t iEnc = quant(vEnc, maxval);
        
        fprintf(wr, " %6d %15.12f %15.12f %15.12f %15.12f %6d", iRaw, vRaw, vLin, vEnc, zErr, iEnc);
        fprintf(wr, "  ");

        double mRaw = vlog(vRaw, eRaw);
        double mLin = vlog(vLin, eLin);
        double mEnc = vlog(vEnc, eRaw);
        
        fprintf(wr, "   %15.12f %15.12f %15.12f", mRaw, mLin, mEnc);

        fprintf(wr, "\n");
           
      }
    fclose(wr);
    free(fname);
    return;

    double vlog(float v, double e)
      { double x = fabs(v);
        double y;
        if (x == 0)
          { y = NAN; }
        else
          { y = log(x/e)/M_LN10 * (v < 0 ? -1 : +1); }
        return y;
      }

    int32_t quant(float v, int32_t maxval)
      { double x = fabs(v);
        double y = (v < 0 ? -1 : +1)*floor(x*maxval + 0.5);
        return (int32_t)y;
      }
   }

options_t *teng_parse_options(int32_t argc, char **argv)
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

    argparser_get_keyword(pp, "-maxval");
    o->maxval = (int32_t)argparser_get_next_int(pp, 1, 65535);
    
    argparser_get_keyword(pp, "-kind");
    o->kind = argparser_get_next_non_keyword(pp);
    
    /* Defaults even if not used: */
    o->expo = 1.000; o->bias = 0.000;
    o->U = double_vec_new(0);
    o->V = double_vec_new(0);
    
    if (strcmp(o->kind, "generic") == 0)
      { argparser_get_keyword(pp, "-expo");
        o->expo = argparser_get_next_double(pp, 1.0e-100, 1000.0);
        
        argparser_get_keyword(pp, "-bias");
        o->bias = argparser_get_next_double(pp, 0.0, 0.999999);
      }
    else if (strcmp(o->kind, "interp") == 0)
      { argparser_get_keyword(pp, "-points");
        teng_parse_pair_list(pp, &(o->U), &(o->V));
        affirm(o->U.ne == o->V.ne, "bug in {teng_parse_pair_list}"); 
      }
    else if ((strcmp(o->kind, "BT709") != 0) && (strcmp(o->kind, "sRGB") != 0))
      { argparser_error(pp, "invalid \"-kind\" argument"); }

    /* PARSE POSITIONAL ARGUMENTS: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */
    argparser_finish(pp);

    return o;
  }

void teng_parse_pair_list(argparser_t *pp, double_vec_t *U, double_vec_t *V)
  { uint32_t np = 0;
    while(argparser_next_is_number(pp))
      { int32_t ip = (int32_t)np; np++;
        for (int32_t j = 0; j < 2; j++)
          { double_vec_t *w = ( j == 0 ? U : V );
            double_vec_expand(w, (vec_index_t)ip);
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
