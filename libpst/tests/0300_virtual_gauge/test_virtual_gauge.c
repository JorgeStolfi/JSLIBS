#define PROG_NAME "compare_gauge_virtual"
#define PROG_DESC "compares a simple lighting model to a spherical light gauge image"

/* Last edited on 2025-01-30 05:01:13 by stolfi */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <jsfile.h>
#include <jsprintf.h>
#include <math.h>
#include <r3.h>
#include <rn.h>
#include <box.h>
#include <rmxn.h>
#include <ellipse_crs.h>
#include <ellipse_crs_args.h>
#include <affirm.h>
#include <sample_conv_gamma.h>
#include <float_image.h>
#include <float_image_read_pnm.h>
#include <float_image_write_gen.h>
#include <argparser.h>
#include <jsrandom.h>

#include <pst_virtual_gauge.h>

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "   -inFile  {IMAGE_FILE_IN} \\\n" \
  "   -gauge \\\n" \
  "       " pst_virtual_gauge_args_parse_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION" \
  "  ???...\n" \
  "\n" \
  "    The \"view\" option, if present, must immediately follow the other gauge parameters." \
  "OPTIONS\n" \
  "  " pst_virtual_gauge_args_parse_HELP_INFO "\n" \
  "\n" \
  "  " pst_light_spec_HELP_INFO ""

typedef struct options_t
  { pst_virtual_gauge_data_t *gd;
    pst_light_t *lht;
    char *file_in;          /* Filename of image containing the gauge, or {NULL}. */
    char *outPrefix;
    double gammaDec;
  } options_t;

options_t* parse_options(int32_t argc, char **argv);

float_image_t* read_gauge_image(char *inFile, double gammaDec);

void tvga_write_image(char *prefix, char *tag, float_image_t *img, float vmin, float vmax);
  /* Writes the image {img} in FNI format as "{prefix}-{tag}.fni" and in PPM format
    as "{prefix}-{tag}.ppm".  In the latter, the sample range {[vmin _ vmax]} is 
    mapped to {0..maxval}.  */

#define sRGB_ENC_EXPO sample_conv_gamma_sRGB_ENC_EXPO
#define sRGB_DEC_EXPO sample_conv_gamma_sRGB_DEC_EXPO
#define sRGB_BIAS sample_conv_gamma_sRGB_BIAS

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = parse_options(argc,argv);

    int32_t NX, NY;
    int32_t NCS = 4; /* Channels in synthetic image. */
    float_image_t *img_in;
    if (o->file_in != NULL)
      { fprintf(stderr,"reading gauge image...\n");
        double expoDec = sRGB_ENC_EXPO;
        img_in = read_gauge_image(o->file_in, expoDec);
        int32_t NCI, NXI, NYI;
        float_image_get_size(img_in, &NCI, &NXI, &NYI);
        demand(NCI == 3, "input image must be RGB (3 channels)");
        NX = NXI; NY = NYI;
      }
    else
      { img_in = NULL; NX = 600; NY = 400; }

    auto frgb_t shading(r3_t *nrm);
      /* Shading function for the virtual gauge. */

    fprintf(stderr,"making virtual gauge image...\n");
    frgb_t bg = frgb_Black;
    float_image_t *img_sy = float_image_new(NCS, NX, NY); /* Synthetic image with opacity channel. */
    float_image_t *img_nr = float_image_new(3, NX, NY);     /* Normal map. */
    for (int32_t c = 0; c < 3; c++) { float_image_fill_channel(img_sy, c, bg.c[c]); }
    pst_virtual_gauge_paint(o->gd, shading, img_sy, img_nr);
    tvga_write_image(o->outPrefix, "sy", img_sy, 0.0, 1.0);
    tvga_write_image(o->outPrefix, "nr", img_nr, -1.0, +1.0);

    if (img_in != NULL)
      { fprintf(stderr,"computing error image...");
        float_image_t *img_nr = NULL; /* For now */
        float_image_t *img_er = pst_shading_difference_image(img_in, img_sy, NULL, img_nr);
        tvga_write_image(o->outPrefix, "er", img_er, -1.0, +1.0);
      }

    return 0;

    frgb_t shading(r3_t *nrm)
      {
        return pst_light_shading(o->lht, nrm);
      }
  }

float_image_t* read_gauge_image(char *inFile, double expoDec)
  { char *ext = strrchr(inFile, '.');
    demand(ext != NULL, "image name has no extension");
    float_image_t *img;
    if ((strcmp(ext, ".ppm") == 0) || (strcmp(ext, ".pgm") == 0))
      { bool_t isMask = FALSE;
        double bias = 0.0; /* For now */
        bool_t yup = TRUE;
        bool_t warn = FALSE;
        bool_t verbose = TRUE;
        img = float_image_read_pnm_named
          ( inFile, isMask, expoDec, bias, yup, warn, verbose );
      }
    else if (strcmp(ext, ".fni") == 0)
      { FILE *rd = open_read(inFile, TRUE);
        img = float_image_read(rd);
        if (rd != stdin) { fclose(rd); }
      }
    else
     { demand(FALSE, "unknown image name extension"); }
    return img;
  }

void tvga_write_image(char *prefix, char *tag, float_image_t *img, float vmin, float vmax)
  { int32_t NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    char *fni_name = jsprintf("%s-%s.fni", prefix, tag);
    float_image_write_named(fni_name, img);
    free(fni_name);

    char *png_name = jsprintf("%s-%s.png", prefix, tag);
    bool_t yUp = FALSE;
    double expoEnc = sRGB_ENC_EXPO;
    double bias = sRGB_BIAS;
    bool_t verbose = TRUE;
    float_image_write_gen_named(png_name, img, image_file_format_PNG, yUp, vmin, vmax, expoEnc, bias, verbose);
    free(png_name);
  }

options_t* parse_options(int32_t argc, char **argv)
  {
    options_t *o = (options_t*)malloc(sizeof(options_t));
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    if (argparser_keyword_present(pp, "-file"))
      { o->file_in = argparser_get_next_non_keyword(pp); }
    else
      { o->file_in = NULL; }

    argparser_get_keyword(pp, "-gauge");
    o->gd = talloc(1, pst_virtual_gauge_data_t);
    pst_virtual_gauge_args_parse(pp, o->gd);

    argparser_get_keyword(pp, "-light");
    o->lht = pst_light_spec_parse(pp, FALSE);

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);

    if (argparser_keyword_present(pp, "-gammaDec"))
        { o->gammaDec = argparser_get_next_double(pp, 0.0, 10.0); }
      else
        { o->gammaDec = 1.0; }

    argparser_finish(pp);
    return o;
  }
