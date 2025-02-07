#define PROG_NAME "test_paint"
#define PROG_DESC "test of {float_image_paint.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-20 18:18:42 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_paint_C_COPYRIGHT \
  "Copyright � 2007  by the State University of Campinas (UNICAMP)"

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
  "  Criado em 2007-08-01 por J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_paint_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program writes a PPM test image showing the" \
  " result of various {float_image_paint.h} procedures."

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
#include <jsmath.h>
#include <jsprintf.h>
#include <jsrandom.h>
#include <ellipse_ouv.h>
#include <ellipse_crs.h>
#include <sample_conv.h>
#include <sample_conv_gamma.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <float_image_to_uint16_image.h>
#include <float_image_paint.h>
#include <float_image.h>

#define BT_ENC_EXPO sample_conv_gamma_BT709_ENC_EXPO
#define BT_ENC_BIAS sample_conv_gamma_BT709_BIAS
  /* Values of {expo} and {bias} for {sample_conv_gamma} 
    that approximate the BT.709 encoding. */

typedef struct options_t
  { 
  } options_t;

options_t *fitp_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

int32_t main(int32_t argc, char **argv);

float_image_t *fitp_make_test_image(int32_t nx, int32_t ny, uint32_t m, options_t *o);
  /* Creates a test image using various
    {float_image_paint.h} tools. */

void fitp_write_image(char *name, float_image_t *A);
  /* Writes the image {A} to a file called "{name}.{ext}",
    in the PPM/PGM format, without gamma correction. */

void test_paint_dot
  ( float_image_t *A,
    int32_t channel, 
    double xctr,
    double yctr, 
    double rad,
    double hwd, 
    float vfill,
    float vdraw,
    uint32_t m
  );

void test_paint_ellipse
  ( float_image_t *A,
    int32_t channel, 
    double xctr,
    double yctr, 
    double rad,
    double hwd, 
    float vfill,
    float vdraw,
    uint32_t m
  );

void test_paint_cross
  ( float_image_t *A,
    int32_t channel, 
    double xctr,
    double yctr, 
    double rad,
    double hwd, 
    float vfill,
    float vdraw,
    uint32_t m
  );

void test_paint_smudge
  ( float_image_t *A,
    int32_t channel, 
    double xctr,
    double yctr, 
    double rad,
    float vfill,
    uint32_t m
  );

void test_paint_rectangle
  ( float_image_t *A,
    int32_t channel,  
    double xctr,
    double yctr, 
    double rad, 
    double hwd,
    float vfill,
    float vdraw,
    uint32_t m
  );

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = fitp_parse_options(argc, argv);

    char *outPrefix = "out/img"; /* Prefix of output file names. */

    int32_t nx = 640;
    int32_t ny = 480;
    
    uint32_t m; /* Subsampling order. */
    for (m = 0; m <= 3; m++)
      { /* Generate the test image {A} with antialiasing {m}: */
        float_image_t *A = fitp_make_test_image(nx, ny, m, o);
        char *filename = jsprintf("%s-%02d", outPrefix, m);
        int32_t c;
        for (c = 0; c < 3; c++) 
          { float_image_apply_gamma(A, c, BT_ENC_EXPO, BT_ENC_BIAS); }
        fitp_write_image(filename, A);
        free(filename);
        float_image_free(A); A = NULL;
      }
      
     free(o); o = NULL;

    return 0;
  }

float_image_t *fitp_make_test_image(int32_t nx, int32_t ny, uint32_t m, options_t *o)
  {
    /* Create image, fill it with  black: */
    float_image_t *A = float_image_new(3, nx, ny);
    float_image_fill_channel(A, 0, 0.300f);
    float_image_fill_channel(A, 1, 0.280f);
    float_image_fill_channel(A, 2, 0.150f);
    srandom(4615 + 418*747);

    /* Paint various shapes: */
    double rad_max = 15;
    int32_t step = (int32_t)ceil(2.5*rad_max);
    
    uint32_t ntypes = 5; /* See below. */
    
    uint32_t trial = 0;
    for (int32_t tx = step/2; tx < nx; tx += step)
      { for (int32_t ty = step/2; ty < ny; ty += step)
          { bool_t sync = (trial < 3*ntypes); /* If TRUE, uses only integer or half-integer centers. */
            fprintf(stderr, "\n");

            /* Choose the channel where to paint: */
            int32_t channel = int32_abrandom(0,2);

            /* Choose the center of the symbol: */
            double xctr = tx + 2.0*drandom() - 1;
            double yctr = ty + 2.0*drandom() - 1;
            if (sync)
              { if (drandom() < 0.5)
                  { /* Sync to integer: */
                    xctr = (int32_t)floor(xctr + 0.5);
                    yctr = (int32_t)floor(yctr + 0.5);
                  }
                else
                  { /* Sync to half-integer: */
                    xctr = (int32_t)floor(xctr) + 0.5;
                    yctr = (int32_t)floor(yctr) + 0.5;
                  }
              }

            /* Choose some parameters used by more than one shape: */
            double rad = drandom()*rad_max;
            double hwd = drandom()*4;
            float vdraw = (float)(drandom() < 0.250 ? NAN : 0.600 + 0.400*drandom());
            float vfill = (float)(drandom() < 0.250 ? NAN : 0.400 + 0.600*drandom());

            uint32_t mtype = trial % ntypes;

            fprintf(stderr, " channel = %d", channel);
            fprintf(stderr, " type = %d", mtype);

            switch (mtype)
              { case 0:
                  test_paint_ellipse(A, channel, xctr, yctr, rad, hwd, vfill, vdraw, m);
                  break;
                case 1:
                  test_paint_cross(A, channel, xctr, yctr, rad, hwd, vfill, vdraw, m);
                  break;
                case 2:
                  test_paint_dot(A, channel, xctr, yctr, rad, hwd, vfill, vdraw, m);
                  break;
                case 3:
                  test_paint_smudge(A, channel, xctr, yctr, rad, vfill, m);
                  break;
                case 4:
                  test_paint_rectangle(A, channel, xctr, yctr, rad, hwd, vfill, vdraw, m);
                  break;
                default:
                  assert(FALSE);
              }
              
            trial++;
          }
      }
    return A;
  }

void test_paint_dot
  ( float_image_t *A,
    int32_t channel, 
    double xctr,
    double yctr, 
    double rad,
    double hwd, 
    float vfill,
    float vdraw,
    uint32_t m
  )
  { /* Paint a dot: */
    fprintf(stderr, " [dot]\n");

    bool_t round = (drandom() < 0.333);
    bool_t diagonal = (drandom() < 0.500);

    fprintf(stderr, "  ");
    fprintf(stderr, " ctr = ( %9.5f %9.5f )", xctr, yctr);
    fprintf(stderr, " rad = %9.5f", rad);
    fprintf(stderr, " round = %c", "FT"[round]);
    fprintf(stderr, " diag = %c", "FT"[diagonal]);
    fprintf(stderr, "  ");
    fprintf(stderr, " hwd = %9.5f", hwd);
    fprintf(stderr, " vfill=%5.3f", vfill);
    fprintf(stderr, " vdraw = %5.3f", vdraw);

    double w_tot = float_image_paint_dot
      ( A, channel, xctr, yctr, rad, hwd, round, diagonal, vfill, vdraw, m );
    
    fprintf(stderr, " w_tot = %8.4f", w_tot);
    fprintf(stderr, "\n");
  }

void test_paint_ellipse
  ( float_image_t *A,
    int32_t channel, 
    double xctr,
    double yctr, 
    double rad,
    double hwd, 
    float vfill,
    float vdraw,
    uint32_t m
  )
  { /* Paint an ellipse: */
    fprintf(stderr, " [ellipse]\n");

    ellipse_crs_t E;
    E.rad = 2*rad;
    E.ctr = (r2_t) {{ xctr, yctr }};
    E.str.c[0] = (2*drandom()-1)*20;
    E.str.c[1] = (2*drandom()-1)*20;

    fprintf(stderr, "  ");
    fprintf(stderr, " ctr = ( %9.5f %9.5f )", xctr, yctr);
    fprintf(stderr, " rad = %9.5f", E.rad);
    fprintf(stderr, "  ");
    fprintf(stderr, " hwd = %5.3f", hwd);
    fprintf(stderr, " vfill = %5.3f", vfill);
    fprintf(stderr, " vdraw = %5.3f", vdraw);

    /* Paint ellipse: */
    double w_tot = float_image_paint_ellipse_crs
      ( A, channel, &E, hwd, vfill, vdraw, m );

    fprintf(stderr, " w_tot = %8.4f", w_tot);
    fprintf(stderr, "\n");
  }

void test_paint_cross
  ( float_image_t *A,
    int32_t channel, 
    double xctr,
    double yctr, 
    double rad,
    double hwd, 
    float vfill,
    float vdraw,
    uint32_t m
  )
  { /* Paint a cross: */
    fprintf(stderr, " [cross]\n");

    bool_t empty = (drandom() < 0.500);
    bool_t diagonal = (drandom() < 0.500);

    fprintf(stderr, "  ");
    fprintf(stderr, " ctr = ( %9.5f %9.5f )", xctr, yctr);
    fprintf(stderr, " rad = %9.5f", rad);
    fprintf(stderr, " empty = %c", "FT"[empty]);
    fprintf(stderr, " diagonal = %c", "FT"[diagonal]);
    fprintf(stderr, "  ");
    fprintf(stderr, " hwd = %9.5f", hwd);
    fprintf(stderr, " vdraw = %5.3f", vdraw);

    double w_tot = float_image_paint_cross(A, channel, xctr, yctr, rad, empty, hwd, diagonal, vdraw, m);

    fprintf(stderr, " w_tot = %8.4f", w_tot);
    fprintf(stderr, "\n");
  }

void test_paint_smudge
  ( float_image_t *A,
    int32_t channel, 
    double xctr,
    double yctr, 
    double rad, 
    float vfill,
    uint32_t m
  )
  { /* Paint a smudge: */
    fprintf(stderr, " [smudge]\n");

    double stretch = exp(2*drandom() - 1);
    double xdev = rad/stretch;
    double ydev = rad*stretch;

    fprintf(stderr, "  ");
    fprintf(stderr, " ctr = ( %9.5f %9.5f )", xctr, yctr);
    fprintf(stderr, " xdev = %9.5f", xdev);
    fprintf(stderr, " ydev = %9.5f", ydev);
    fprintf(stderr, "  ");
    fprintf(stderr, " vfill=%5.3f", vfill);

    double w_tot = float_image_paint_smudge(A, channel, xctr, yctr, xdev, ydev, vfill, m);
    
    fprintf(stderr, " w_tot = %8.4f", w_tot);
    fprintf(stderr, "\n");
  }

void test_paint_rectangle
  ( float_image_t *A,
    int32_t channel, 
    double xctr,
    double yctr, 
    double rad,
    double hwd, 
    float vfill,
    float vdraw,
    uint32_t m
  )
  { /* Paint a dot: */
    fprintf(stderr, " [rectangle]\n");

    double xrad = rad*drandom();
    double yrad = rad*drandom();

    double xmin = xctr - xrad;
    double xmax = xctr + xrad;
    double ymin = yctr - yrad; 
    double ymax = yctr + yrad; 

    fprintf(stderr, "  ");
    fprintf(stderr, " xmin = %9.5f", xmin);
    fprintf(stderr, " xmax = %9.5f", xmax);
    fprintf(stderr, " ymin = %9.5f", ymin);
    fprintf(stderr, " ymax = %9.5f", ymax);
    fprintf(stderr, "  ");
    fprintf(stderr, " hwd = %9.5f", hwd);
    fprintf(stderr, " vfill=%5.3f", vfill);
    fprintf(stderr, " vdraw = %5.3f", vdraw);

    double w_tot = float_image_paint_rectangle
      ( A, channel, xmin, xmax, ymin, ymax, hwd, vfill, vdraw, m );

    fprintf(stderr, " w_tot = %8.4f", w_tot);
    fprintf(stderr, "\n");
  }


void fitp_write_image(char *name, float_image_t *A)
  { char *suff = (A->sz[0] == 1 ? ".pgm" : ".ppm"); 
    char *fname = jsprintf("%s%s", name, suff);
    FILE *wr = open_write(fname, TRUE);
    int32_t chns = (int32_t)A->sz[0];
    bool_t yup = TRUE, verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    uint16_image_t *pimg = float_image_to_uint16_image(A, isMask, chns, NULL, NULL, NULL, 255, yup,verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pimg, forceplain, verbose);
    uint16_image_free(pimg);
    fclose(wr);
    free(fname);
  }

options_t *fitp_parse_options(int32_t argc, char **argv)
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

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }
