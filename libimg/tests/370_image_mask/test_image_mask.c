#define PROG_NAME "test_image_mask"
#define PROG_DESC "test of {float_image_mask.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-01-21 18:31:21 by stolfi */

#define test_image_mask_C_COPYRIGHT \
  "Copyright � 2007  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    {NX} {NY} \\\n" \
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
  "  2012-01-15 J. Stolfi: Added size options.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_image_mask_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program writes PGM masks with {NX} columns and {NY} rows, generated by calling" \
  " {float_image_mask.h} procedures with various parameters."

#define PROG_INFO_OPTS \
  "  None."

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
#include <jsprintf.h>
#include <jsrandom.h>
#include <sample_conv.h>
#include <sample_conv_gamma.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <float_image_to_uint16_image.h>
#include <float_image_mask.h>
#include <float_image.h>

#define BT_ENC_EXPO sample_conv_gamma_BT709_ENC_EXPO
#define BT_ENC_BIAS sample_conv_gamma_BT709_BIAS
  /* Values of {expo} and {bias} for {sample_conv_gamma} 
    that approximate the BT.709 encoding. */

typedef struct options_t
  { int NX;
    int NY;
  } options_t;

options_t *tim_parse_options(int argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

int main(int argc, char **argv);

float_image_t *tim_make_test_image(int NX, int NY, int m, options_t *o);
  /* Creates a test image using various
    {float_image_mask.h} tools. */

void tim_show_mask_statistics(float_image_t *msk, int ic);
  /* Prints statistical proeprties of channel {ic} of {msk}. */

void tim_write_image(char *name, float_image_t *A);
  /* Writes the image {A} to a file called "{name}.{ext}",
    in the PPM/PGM format, without gamma correction. */

int main(int argc, char **argv)
  {
    /* Parse the command line options: */
    options_t *o = tim_parse_options(argc, argv);

    char *prefix = "out/msk"; /* Prefix of output file names. */

    int NC = 1;
    int NX = o->NX;
    int NY = o->NY;
    /* Allocate the mask: */
    float_image_t *msk = float_image_new(NC, NX, NY);
    
    int ord;    /* Continuity order (before blurring). */
    int ic = 0; /* Channel to use. */
    for (ord = 0; ord <= 2; ord++)
      { int iround; /* Support: 0 = rectangular, 1 = oval. */
        for (iround = 0; iround <= 1; iround++)
          { int imodf; /* Modification: 0 = none, 1 = Gaussian, 2 = InvPower. */
            for (imodf = 0; imodf <= 2; imodf++)
              { fprintf(stderr, "--------------------------------------------------\n");
                fprintf
                  ( stderr, "basic order = %d %s\n",
                    ord, (char *[3]){"Flat", "Para", "Hann"}[ord]
                  );
                fprintf
                  ( stderr, "basic shape = %d %s\n", 
                    iround, (char *[2]){"Rect", "Oval"}[iround]
                  );
                fprintf
                  ( stderr, "mul factor = %d %s\n",
                    imodf, (char *[3]){"Unifm", "Gauss", "Power"}[imodf]
                  );
                /* Store into {msk} the basic pattern for {ord,iround}: */
                float_image_mask_window(msk, ic, ord, (iround != 0));
                /* Multiply by the modifier function: */
                switch (imodf)
                  { case 0: 
                      break;
                    case 1:
                      { double sx = NX/7.0;
                        double sy = NY/5.0;
                        fprintf(stderr, "params = ( %8.3f %8.3f )\n", sx, sy);
                        float_image_mask_mul_gauss(msk, ic, sx, sy);
                      }
                      break;
                    case 2:
                      { double sx = NX/7.0;
                        double sy = NY/5.0;
                        double pwr = 2.0;
                        fprintf
                          ( stderr, "params = ( %8.3f %8.3f ) power = %8.3f\n",
                            sx, sy, pwr
                          );
                        float_image_mask_mul_power(msk, ic, sx, sy, pwr);
                      }
                      break;
                    default:
                      assert(FALSE);
                  }
                /* Compute parameters: */
                tim_show_mask_statistics(msk, ic);
                /* Write mask to disk: */
                int igamma; /* 0 = linear encoding, 1 = BT.709 */
                for (igamma = 0; igamma <= 1; igamma++)
                  { if (igamma != 0)
                      { float_image_apply_gamma(msk, ic, BT_ENC_EXPO, BT_ENC_BIAS); }
                    char *msk_name = jsprintf("%s-%04dx%04d-o%01d-r%01d-m%01d-%c",
                             prefix, NX, NY, ord, iround, imodf, "LG"[igamma]
                    );
                    tim_write_image(msk_name, msk);
                    free(msk_name);
                  }
                fprintf(stderr, "--------------------------------------------------\n");
                fprintf(stderr, "\n");
              }
          }
      }
      
    float_image_free(msk); msk = NULL;
    free(o); o = NULL;

    return 0;
  }

void tim_show_mask_statistics(float_image_t *msk, int ic)
  {
    float_image_mask_stats_t S = float_image_mask_stats_get(msk, ic);
    fprintf(stderr, "\n");
    fprintf(stderr, "mask properties:\n");
    float_image_mask_stats_print(stderr, &S);
    fprintf(stderr, "\n");
  }

void tim_write_image(char *img_name, float_image_t *img)
  { char *suff = (img->sz[0] == 1 ? ".pgm" : ".ppm"); 
    char *fname = jsprintf("%s%s", img_name, suff);
    FILE *wr = open_write(fname, TRUE);
    int chns = (int)img->sz[0];
    bool_t yup = TRUE, verbose = TRUE;
    bool_t isMask = TRUE; /* Assume 0 and 1 are important in encoding/decoding. */
    uint16_image_t *pimg = float_image_to_uint16_image(img, isMask, chns, NULL, NULL, NULL, 255, yup, verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pimg, forceplain, verbose);
    uint16_image_free(pimg);
    fclose(wr);
    free(fname);
  }

options_t *tim_parse_options(int argc, char **argv)
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
    o->NX = (int)argparser_get_next_int(pp, 10, 65535);
    o->NY = (int)argparser_get_next_int(pp, 10, 65535);

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }
