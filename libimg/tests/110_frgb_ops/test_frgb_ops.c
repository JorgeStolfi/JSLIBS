#define PROG_NAME "test_frgb_ops"
#define PROG_DESC "test of various ops from {frgb_ops.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-03-07 14:00:31 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_frgb_ops_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
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
  "  " test_frgb_ops_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  This program reads nothing and writes some messages."

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

#include <frgb_ops.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <bool.h>
#include <affirm.h>
#include <argparser.h>

typedef struct options_t
  { bool_t op; /* placeholder. */
  } options_t;
  
typedef void tfop_tr_t(frgb_t *p);
typedef double tfop_fn_t(frgb_t *p);

int main(int argc, char **argv);

options_t *tfop_parse_options(int argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

void tfop_do_tests(options_t *o);
  /* Creates a test image using various
    {float_image_gradient.h} tools. */
  
void tfop_test_transform(tfop_tr_t *fpre, char *sdom, char *simg, tfop_tr_t *fdir, tfop_tr_t *finv);
  /* Tests if a pair of functions {fdir}, {finv} are mutually inverse,
    namely whether {finv(fdir(p)) == p} for a bunch of triples {p}.
    
    The strings {sdom} and {simg} are the names of the color spaces for
    the arguments of {fdir} and {finv}, respectively. If {fpre} is
    {NULL}, the argument of {fdir} is in the RGB color space, hence in
    the unit RGB cube; then {sdom} should be "RGB". If {fpre} is not
    {NULL}, assumes that the argument of {fdir} is in some other color
    space {sdom}, and {fpre} transforms valid RGB values to valid {sdom}
    values. */

void tfop_test_float_fn(char *simg, tfop_tr_t *fdir, int k, char *fn_name, tfop_fn_t *fval);
  /* Tests if a colorspace transformation {fdir} from RGB to some colorspace {simg}
    is consistent with a single-value projection function {fval}; namely, whether {fval(p) = fdir(p).c[k]} for 
    a bunch of non-gray RGB triples {p}. */

bool_t tfop_check_tr(frgb_t *rgb, char *sdom, frgb_t *fp, char *simg, frgb_t *fq, frgb_t *fr);
  /* Compares {*fp} with {*fr}. Complains and returns true if not close enough.
    Assumes that {fq} is the result of mapping {fp} from space {sdom} to space {simg},
    and {fr} is the result of mapping {fq} back to space {sdom}.  The color's RGB coordinates
    are assumed to be {rgb}. */
  
bool_t tfop_check_fn(frgb_t *rgb, char *sdom, frgb_t *fq, int k, char *fn_name, double fvalu);
  /* Compares {fq.c[k]} with {fvalu}. Complains and returns true if not close enough.
    assumes that {fq} is the result of transforming the triple {rgb} from RGB space
    to {simg} space. Also assumes that {fvalu} is the result of applying {fn_name}
    to {rgb}. */

int main(int argc, char **argv)
  {
    /* Parse the command line options: */
    options_t *o = tfop_parse_options(argc, argv);

    tfop_do_tests(o);
    free(o); o = NULL;

    return 0;
  }

void tfop_do_tests(options_t *o)
  {
    /*
    frgb_t frgb_mix(double ca, frgb_t *a, double cb, frgb_t *b);
    frgb_t frgb_scale(double s, frgb_t *a);
    frgb_t frgb_add(frgb_t *a, frgb_t *b);
    frgb_t frgb_sub(frgb_t *a, frgb_t *b);
    frgb_t frgb_mul(frgb_t *a, frgb_t *b);
    bool_t frgb_is_all_zeros(frgb_t *a);
    bool_t frgb_is_all_ones(frgb_t *a);
    bool_t frgb_eq(frgb_t *a, frgb_t *b);
    
    double frgb_gamma_encoding_gray(double y, double gamma, double bias);
    double frgb_gamma_decoding_gray(double y, double gamma, double bias);
    
    typedef frgb_t frgb_adjuster_t(frgb_t *p, int col, int row);
    frgb_t frgb_correct_arg(frgb_t *p, frgb_t *inGamma, int gray);
    double frgb_log_scale_gray(double x);
    void frgb_log_scale(frgb_t *p, int chns);
    double frgb_clip_gray(double p);
    void frgb_clip_rgb(frgb_t *p);
    void frgb_clip_rgb_towards(frgb_t *p, frgb_t *q);
    void frgb_clip_rgb_towards_grey(frgb_t *p);
    double frgb_apply_kappa_gray(double y, double kappa);
    void frgb_apply_glob_kappa_sat_clip(frgb_t *p, double kap, double satf);
    int frgb_dequal(double *a, double *b, int chns);
    int frgb_fequal(float *a, float *b, int chns);
    double frgb_floatize(int ival, int maxval, double zero, double scale);
    int frgb_quantize(double fval, double zero, double scale, int maxval);
    frgb_t frgb_parse(argparser_t *pp, double lo, double hi);
    frgb_t frgb_read(FILE *rd, double lo, double hi);
    frgb_t frgb_parse_color(argparser_t *pp);
    frgb_t frgb_read_color(FILE *rd);

    double frgb_Y_pbm(frgb_t *p);
    
    void frgb_YUV_to_yuv(frgb_t *p, double ybias);
    void frgb_YUV_from_yuv(frgb_t *p, double ybias);
    
    void frgb_YUV_to_Yuv(frgb_t *p, double ybias);
    void frgb_YUV_from_Yuv(frgb_t *p, double ybias);

    void frgb_print(FILE *f, char *pref, frgb_t *p, int chns, char *fmt, char *suff);
    void frgb_print_int_pixel(FILE *f, char *pref, int *p, int chns, char *suff);
    void frgb_debug(char *label, int col, int row, frgb_t *p, int chns, char *tail);
    void frgb_debug_int_pixel(char *label, int col, int row, int *p, int chns, char *tail);
    */
    
    
    tfop_test_transform(NULL, "RGB", "CIE_XYZrec601_1", &frgb_to_CIE_XYZrec601_1, &frgb_from_CIE_XYZrec601_1);
    tfop_test_float_fn("CIE_XYZrec601_1", &frgb_to_CIE_XYZrec601_1, 1, "luminance_CIE_XYZrec601_1", &frgb_luminance_CIE_XYZrec601_1);
    
    tfop_test_transform(NULL, "RGB", "CIE_XYZccir709", &frgb_to_CIE_XYZccir709, &frgb_from_CIE_XYZccir709);
    tfop_test_float_fn("CIE_XYZccir709", &frgb_to_CIE_XYZccir709, 1, "luminance_CIE_XYZccir709", &frgb_luminance_CIE_XYZccir709);
    
    tfop_test_transform(NULL, "RGB", "CIE_XYZitu_D65", &frgb_to_CIE_XYZitu_D65, &frgb_from_CIE_XYZitu_D65);
    tfop_test_float_fn("CIE_XYZitu_D65", &frgb_to_CIE_XYZitu_D65, 1, "luminance_CIE_XYZitu_D65", frgb_luminance_CIE_XYZitu_D65);
    
    tfop_test_transform(NULL, "RGB", "YUV", &frgb_to_YUV, &frgb_from_YUV);
    tfop_test_float_fn("YUV", &frgb_to_YUV, 0, "get_Y", &frgb_get_Y);
    
    tfop_test_transform(&frgb_to_YUV, "YUV", "YHS", &frgb_YUV_to_YHS, &frgb_YHS_to_YUV);
    
    tfop_test_transform(NULL, "RGB", "HSV_CG", &frgb_to_HSV_CG, &frgb_from_HSV_CG);
    
    tfop_test_transform(NULL, "RGB", "HTY_UV", &frgb_to_HTY, &frgb_from_HTY);
    tfop_test_float_fn("HTY_UV", &frgb_to_HTY, 0, "get_H", frgb_get_H);
    
    tfop_test_transform(NULL, "RGB", "YIQ", &frgb_to_YIQ, &frgb_from_YIQ);
    tfop_test_float_fn("YIQ", &frgb_to_YIQ, 0, "get_Y", &frgb_get_Y);
    
    tfop_test_transform(NULL, "RGB", "YCbCr_601_1", &frgb_to_YCbCr_601_1, &frgb_from_YCbCr_601_1);
    tfop_test_float_fn("YCbCr_601_1", &frgb_to_YCbCr_601_1, 0, "get_Y", &frgb_get_Y);
    
    tfop_test_transform(NULL, "RGB", "YUV_a", &frgb_to_YUV_a, &frgb_from_YUV_a);
    tfop_test_float_fn("YUV_a", &frgb_to_YUV_a, 0, "get_Y", &frgb_get_Y);
    
    tfop_test_transform(NULL, "RGB", "YUV_b", &frgb_to_YUV_b, &frgb_from_YUV_b);
    tfop_test_float_fn("YUV_b", &frgb_to_YUV_b, 0, "get_Y", &frgb_get_Y);
  }

void tfop_test_transform(tfop_tr_t *fpre, char *sdom, char *simg, tfop_tr_t *fdir, tfop_tr_t *finv)
  {
    fprintf(stderr, "--- tfop_test_transform: frgb");
    if (fpre != NULL) { fprintf(stderr, "_%s", sdom); }
    fprintf(stderr, "_{to,from}_%s ---\n", simg); 
    
    int N = 11; /* Better be odd? */
    int bugs = 0;
    int tics = 0;
    int ir, ig, ib;
    frgb_t rgb, fp, fq, fr;
    for (ir = 0; ir < N; ir++)
      { for (ig = 0; ig < N; ig++) 
          { for (ib = 0; ib < N; ib++)
              { rgb.c[0] = (float)(0.0001 + 0.9998*((double)ir)/(N - 1 + 1.0e-200));
                rgb.c[1] = (float)(0.0001 + 0.9998*((double)ig)/(N - 1 + 1.0e-200));
                rgb.c[2] = (float)(0.0001 + 0.9998*((double)ib)/(N - 1 + 1.0e-200));
                fp = rgb; /* Argument of the direct function. */
                if (fpre != NULL) { fpre(&fp); }
                fq = fp; fdir(&fq); /* Result of the direct function. */
                fr = fq; finv(&fr); /* Result of the inverse function. */
                bool_t res_tr = tfop_check_tr(&rgb, sdom, &fp, simg, &fq, &fr); 
                if (! res_tr) { bugs++; if (bugs > 100) { exit(1); } }
                tics++;
              }
          }
      }
    fprintf(stderr, "tested %d points in the RGB cube, %d errors\n", tics, bugs); 
  }
  
bool_t tfop_check_tr(frgb_t *rgb, char *sdom, frgb_t *fp, char *simg, frgb_t *fq, frgb_t *fr)
  {
    for (int c = 0; c < 3; c++)
      {if (fabs(fp->c[c] - fr->c[c]) > 0.0001) 
          { char *sa = (strcmp(sdom, "RGB") == 0  ? "" : "_");
            char *sb = (strcmp(sdom, "RGB") == 0  ? "" : sdom);
            fprintf(stderr, "** inconsistency for frgb%s%s_{to,from}_%s:\n", sa, sb, simg);
            frgb_print(stderr, "rgb = ", rgb, 3, "%7.4f", "\n");
            frgb_print(stderr, "fp =  ", fp, 3, "%7.4f", "\n");
            frgb_print(stderr, "fq =  ", fq, 3, "%7.4f", "\n");
            frgb_print(stderr, "fr =  ", fr, 3, "%7.4f", "\n");
            return FALSE;
          }
      }
    return TRUE;
  }

void tfop_test_float_fn(char *simg, tfop_tr_t *fdir, int k, char *fn_name, tfop_fn_t *fval)
  {
    fprintf(stderr, "--- tfop_test_float_fn: frgb_to_%s.c[%d] and frgb_%s ---\n", simg, k, fn_name);

    int N = 11; /* Better be odd? */
    int bugs = 0;
    int tics = 0;
    int ir, ig, ib;
    frgb_t rgb, fq;
    for (ir = 0; ir < N; ir++)
      { for (ig = 0; ig < N; ig++) 
          { for (ib = 0; ib < N; ib++)
              { rgb.c[0] = (float)(0.0001 + 0.9998*((double)ir)/(N - 1 + 1.0e-200));
                rgb.c[1] = (float)(0.0001 + 0.9998*((double)ig)/(N - 1 + 1.0e-200));
                rgb.c[2] = (float)(0.0001 + 0.9998*((double)ib)/(N - 1 + 1.0e-200));
                fq = rgb; /* Argument of the transform. */
                fdir(&fq); /* Result of the direct function. */
                bool_t res_fn = TRUE;
                if ((ir != ig) || (ig != ib))
                  { double fvalu = fval(&rgb);
                    res_fn = tfop_check_fn(&rgb, simg, &fq, k, fn_name, fvalu);
                  }
                if (! res_fn) { bugs++; if (bugs > 100) { exit(1); } }
                tics++;
              }
          }
      }
    fprintf(stderr, "tested %d points in the RGB cube, %d errors\n", tics, bugs); 
  }

bool_t tfop_check_fn(frgb_t *rgb, char *simg, frgb_t *fq, int k, char *fn_name, double fvalu)
  {
    if (fabs(fq->c[k] - fvalu) > 0.0001) 
      { fprintf(stderr, "** inconsistency between frgb_to_%s.c[%d] and frgb_%s:\n", simg, k, fn_name);
        frgb_print(stderr, "rgb = ", rgb, 3, "%7.4f", "\n");
        frgb_print(stderr, "fq = ", fq, 3, "%7.4f", "\n");
        fprintf(stderr, "coord fq.c[%d] = %12.9f\n", k, fq->c[k]);
        fprintf(stderr, "function value = %12.9f\n", fvalu);
        return FALSE;
      }
    else
      { return TRUE; }
  }

options_t *tfop_parse_options(int argc, char **argv)
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

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }
