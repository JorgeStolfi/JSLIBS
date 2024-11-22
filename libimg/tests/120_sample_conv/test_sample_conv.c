#define PROG_NAME "test_sample_conv"
#define PROG_DESC "checks the {sample_conv.h} routines"
#define PROG_VERS "1.0"

/* Last edited on 2017-06-22 03:27:52 by stolfilocal */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_sample_conv_C_COPYRIGHT \
  "Copyright © 2010  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -maxval {MAXVAL} \\\n" \
  "    -isMask {ISMASK} \\\n" \
  "    [ -range {LO} {HI} ] \\\n" \
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
  "  test_encode_gamma(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2010-08-14 by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_sample_conv_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program checks the {sample_conv_floatize} and {sample_conv_quantize} routines" \
  " in {sample_conv.h}, with given values of {maxval} and given" \
  " setting of the {isMask} flag.  It tabulates the in/out values" \
  " and checks for proper rounding."

#define PROG_INFO_OPTS \
  "  -maxval {MAXVAL}\n" \
  "  -isMask {ISMASK}\n" \
  "  -range {LO} {HI}\n" \
  "    Specify the parameters {maxval,isMask,lo,hi} for {sample_conv_quantize} and" \
  " {sample_conv_floatize}.  The {MAXVAL} must be a positive" \
  " integer, {(!ISMASK)} is a boolean, and {LO,HI} are float values, defaulting to 0 and 1."

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <argparser.h>
#include <vec.h>
#include <affirm.h>
#include <jspnm.h>
#include <jsmath.h>
#include <bool.h>
#include <sample_conv.h>

#define tsc_MAX_SAMPLE (65535u)

typedef struct options_t
  { sample_uint32_t maxval;  /* User-given {maxval}. */
    bool_t isMask;         /* Rounding/flotation of image samples. */
    double lo;             /* Low end of float sample range. */
    double hi;             /* High end of float sample range. */
  } options_t;

options_t *tsc_parse_options(int argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

int main(int argc, char **argv);

void tsc_test_quantize
  ( float fv, 
    sample_uint32_t maxval, 
    bool_t isMask, 
    double lo, 
    double hi,            /* Input value to map to {maxval}. */
    float *vmin,          /* (IN/OUT) Min float input value seen, or NULL. */
    float *vmax,          /* (IN/OUT) Max float input value seen, or NULL. */
    int *clo,             /* (IN/OUT) Count of input values below {lo}, or NULL. */
    int *chi,             /* (IN/OUT) Count of input values above {hi}, or NULL. */
    sample_uint32_t *imin,  /* (IN/OUT) Min output integer value seen, or NULL. */
    sample_uint32_t *imax   /* (IN/OUT) Max output integer value seen, or NULL. */
  );
  /* Tests {sample_conv_quantize} on {fv} with given parameters. */

void tsc_test_floatize
  ( sample_uint32_t iv,      /* Integer sample value to convert. */
    sample_uint32_t maxval, 
    bool_t isMask, 
    double lo, 
    double hi,             /* Input value to map to {maxval}. */
    sample_uint32_t *imin,   /* (IN/OUT) Min integer input value seen, or NULL. */
    sample_uint32_t *imax,   /* (IN/OUT) Max integer input value seen, or NULL. */
    float *vmin,           /* (IN/OUT) Min float output value seen, or NULL. */
    float *vmax            /* (IN/OUT) Max float output value seen, or NULL. */
  );
  /* Tests {sample_conv_floatize} on {iv} with given parameters. */

int main(int argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = tsc_parse_options(argc, argv);
    
    sample_uint32_t imin, imax;
    float vmin, vmax;
    int clo, chi;
    
    fprintf(stderr, "FLT_MAX = %23.16e\n", FLT_MAX);
    fprintf(stderr, "DBL_MAX = %23.16e\n", DBL_MAX);

    /* Denominator for affine mapping: */
    int den = (o->isMask ? o->maxval : o->maxval+1);

    /* Test the quantization functions: */
    fprintf(stderr, "=== testing sample_conv_quantize");
    fprintf(stderr, "  maxval = %d", (int)o->maxval);
    fprintf(stderr, "  isMask = %c", "FT"[o->isMask]);
    fprintf(stderr, "  range = [ %11.8f _ %11.8f ]", o->lo, o->hi);
    fprintf(stderr, " ===\n");
    fprintf(stderr, "\n");
    
    imin = UINT32_MAX; imax = 0;
    vmin = +INF; vmax = -INF;
    clo = 0; chi = 0;
    
    int kv;
    for (kv = -4; kv <= 2*den + 4; kv++)
      { float fvo = (float)(o->lo + (o->hi - o->lo)*((double)kv)/((double)2*den));
        float fvm = nextafterf(fvo,-INF);
        float fvp = nextafterf(fvo,+INF);
        tsc_test_quantize(fvm, o->maxval, o->isMask, o->lo, o->hi, &vmin, &vmax, &clo, &chi, &imin, &imax);
        tsc_test_quantize(fvo, o->maxval, o->isMask, o->lo, o->hi, &vmin, &vmax, &clo, &chi, &imin, &imax);
        tsc_test_quantize(fvp, o->maxval, o->isMask, o->lo, o->hi, &vmin, &vmax, &clo, &chi, &imin, &imax);
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
    sample_conv_print_quantize_stats(0, 0, vmin, vmax, o->lo, o->hi, clo, chi, o->maxval, imin, imax);
    fprintf(stderr, "\n");
    
    /* Test the floatization functions: */
    fprintf(stderr, "=== testing sample_conv_floatize");
    fprintf(stderr, "  maxval = %d", (int)o->maxval);
    fprintf(stderr, "  isMask = %c", "FT"[o->isMask]);
    fprintf(stderr, "  range = [ %11.8f _ %11.8f ]", o->lo, o->hi);
    fprintf(stderr, " ===\n");
    fprintf(stderr, "\n");
    
    imin = UINT32_MAX; imax = 0;
    vmin = +INF; vmax = -INF;
    
    int iv;
    for (iv = 0; iv <= o->maxval; iv++)
      { tsc_test_floatize(iv, o->maxval, o->isMask, o->lo, o->hi, &imin, &imax, &vmin, &vmax); }
    fprintf(stderr, "\n");
    sample_conv_print_floatize_stats(0, 0, imin, imax, o->maxval, o->lo, o->hi, vmin, vmax);
    fprintf(stderr, "\n");

    /* Cleanup: */
    free(o); o = NULL;
    return 0;
  }

void tsc_test_quantize
  ( float fv, 
    sample_uint32_t maxval, 
    bool_t isMask, 
    double lo, 
    double hi,            /* Input value to map to {maxval}. */
    float *vmin,          /* (IN/OUT) Min float input value seen, or NULL. */
    float *vmax,          /* (IN/OUT) Max float input value seen, or NULL. */
    int *clo,             /* (IN/OUT) Count of input values below {lo}, or NULL. */
    int *chi,             /* (IN/OUT) Count of input values above {hi}, or NULL. */
    sample_uint32_t *imin,  /* (IN/OUT) Min output integer value seen, or NULL. */
    sample_uint32_t *imax   /* (IN/OUT) Max output integer value seen, or NULL. */
  )
  {
    demand(lo < hi, "invalid {lo,hi} range"); /* For now. */
    sample_uint32_t mag = (isMask ? maxval : maxval+1); /* Factor of affine mapping. */
    double off = (double)(isMask ? 00.0 : -0.5);      /* Offset of affine mapping before rounding. */
    double av = ((double)mag)*(fv - lo)/(hi - lo) + off;
    sample_uint32_t iv = sample_conv_quantize(fv, maxval, isMask, lo, hi, vmin, vmax, clo, chi, imin, imax);
    float hv = sample_conv_floatize(iv, maxval, isMask, lo, hi, NULL, NULL, NULL, NULL);
    fprintf(stderr, "  fv = %11.8f = %11.8f/%-5d", fv, ((double)fv)*((double)mag), (int)mag);
    fprintf(stderr, "  av = %11.8f", av);
    fprintf(stderr, "  iv = %5d", (int)iv);
    fprintf(stderr, "  hv = %11.8f = %11.8f/%-5d", hv, ((double)hv)*((double)mag), (int)mag);
    fprintf(stderr, "  e = %11.8f", hv - fv);
    fprintf(stderr, "\n");
    demand(iv <= maxval, "invalid quantized result");
    if (fv <= lo)
      { demand(iv == 0, "invalid quantize result for {fv<=lo}"); }
    else if (fv >= hi)
      { demand(iv == maxval, "invalid quantize result for {fv>=hi}"); }
    else
      { /* Check whether rounded value is within 0.5 of raw affine map: */
        demand(fabs(av - (double)iv) <= 0.5, "improper quantize rounding");
        demand(fabs(hv - fv) <= 0.5/((double)mag)*(hi - lo), "excess error in quantize-floatize");
      }
  }

void tsc_test_floatize
  ( sample_uint32_t iv,      /* Integer sample value to convert. */
    sample_uint32_t maxval, 
    bool_t isMask, 
    double lo, 
    double hi,             /* Input value to map to {maxval}. */
    sample_uint32_t *imin,   /* (IN/OUT) Min integer input value seen, or NULL. */
    sample_uint32_t *imax,   /* (IN/OUT) Max integer input value seen, or NULL. */
    float *vmin,           /* (IN/OUT) Min float output value seen, or NULL. */
    float *vmax            /* (IN/OUT) Max float output value seen, or NULL. */
  )
  {
    demand(iv <= maxval, "invalid integer sample");
    demand(lo < hi, "invalid {lo,hi} range"); /* For now. */
    sample_uint32_t mag = (isMask ? maxval : maxval+1); /* Denominator of affine mapping. */
    double ffo = (double)(isMask ? 00.0 : +0.5);      /* Offset of affine mapping. */
    double av = lo + (((double)iv) + ffo)/((double)mag)*(hi - lo);
    float fv = sample_conv_floatize(iv, maxval, isMask, lo, hi, imin, imax, vmin, vmax);
    sample_uint32_t qv = sample_conv_quantize(fv, maxval, isMask, lo, hi, NULL, NULL, NULL, NULL, NULL, NULL);
    fprintf(stderr, "  iv = %5d", (int)iv);
    fprintf(stderr, "  av = %11.8f", av);
    fprintf(stderr, "  fv = %11.8f = %11.8f/%-5d", fv, ((double)fv)*((double)mag), (int)mag);
    fprintf(stderr, "  qv = %5d", (int)qv);
    fprintf(stderr, "  e = %5d", ((int32_t)qv) - ((int32_t)iv));
    fprintf(stderr, "\n");
    /* Check whether floatized value is within 0.5 of raw affine map: */
    demand(fabs(av - fv) <= 0.5/((double)mag)*(hi - lo), "improper floatize");
    demand(qv - iv == 0, "floatize-quantize is not exact");
  }

options_t *tsc_parse_options(int argc, char **argv)
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
    o->maxval = (sample_uint32_t)argparser_get_next_int(pp, 0, tsc_MAX_SAMPLE);
    
    argparser_get_keyword(pp, "-isMask");
    o->isMask = argparser_get_next_bool(pp);

    if (argparser_keyword_present(pp, "-range"))
      { o->lo = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); 
        o->hi = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); 
      }
    else
      { o->lo = 0; o->hi = 1;  }
      
    /* PARSE POSITIONAL ARGUMENTS: */

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }
