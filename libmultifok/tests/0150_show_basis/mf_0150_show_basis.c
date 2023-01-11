#define PROG_NAME "mf_0150_show_basis"
#define PROG_DESC "test of basis functions from {multifok_focus_op.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-01-08 04:31:47 by stolfi */ 
/* Created on 2012-01-25 by J. Stolfi, UNICAMP */

#define mf_0150_show_basis_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>
#include <bool.h>
#include <float_image_write_pnm.h>
#include <float_image.h>

#include <multifok_focus_op.h>

typedef struct mfsb_options_t 
  { char *outPrefix;   /* Output image file name. */
    int32_t winSize;   /* Size of window. */
  } mfsb_options_t;
  /* Command line parameters. */

int32_t main(int32_t argn, char **argv);

mfsb_options_t *mfsb_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

void mfsb_show_sample_weights
  ( int32_t NW,           /* Window width and height. */
    double ws[],          /* Sample window weights. */
    char *outPrefix       /* Prefix for outupt file name. */
  );
  /* Writes the sample weights {ws[0..NS-1]} as an image file "{outPrefix}-n{NNN}-ws.ppm" where 
    {NNN} is {NW} zero-padded to 3 digits. */

void mfsb_show_single_basis
  ( int32_t NW,                            /* Window width and height. */
    double ws[],                           /* Sample window weights. */
    multifok_focus_op_basis_type_t bType,  /* Basis type. */
    bool_t ortho,                          /* True to ortho normalize basis. */
    char *outPrefix                        /* Prefix for outupt file names. */
  );
  /* The input {timg} and the output {rimg} must be images with the same
    number of channels {NC}. For each channel index {ic} in {0..NC-1},
    sets each pixel of channel {ic} {rimg} to
    {multifok_focus_op_score_simple} applied to a window centered at the
    corresponding pixel of channel {ic} of {timg}, defined by the window
    sample weights {ws[0..NS-1]}.
    
    Uses the basis {phi[0..NB-1]} and coeff weights {wc[0..NB-1]} to compute
    the focus operator. 
    
    If {NB=0} uses a simple focus indicator, ignoring {phi} and {wc} (which may be null)
    
    Writes each basis element {phi[k]} as an image file "{outPrefix}-n{NNN}-b{bType}-o{ortho}-{KKK}.ppm" where 
    {NNN} and {KKK} are {NW} and {k} zero-padded to 3 digits.  Each image pixel is a basis value {phi[kb][ks]}
    times {sqrt(ws[ks])}. */

float_image_t *mfsb_read_image(char *fname);
  /* Reads an image from file "{fname}". */ 

void mfsb_write_image(char *fname, float_image_t *img);
  /* Writes image {img} to file "{fname}". */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfsb_options_t *o = mfsb_parse_options(argc, argv);
    
    int32_t NW = o->winSize;
    demand((NW % 2 == 1) && (NW >= 3), "invalid window size");

    double *ws = multifok_focus_op_sample_weights(NW);
    mfsb_show_sample_weights(NW, ws, o->outPrefix);

    bool_t ortho;
    for (int32_t bt = 0; bt < 3; bt++)
      { multifok_focus_op_basis_type_t bType = (multifok_focus_op_basis_type_t)bt; 
        if ((NW == 3) || (bType != multifok_focus_op_basis_type_DIFF))
          { ortho = FALSE;
            mfsb_show_single_basis(NW, ws, bType, ortho, o->outPrefix);
            ortho = TRUE;
            mfsb_show_single_basis(NW, ws, bType, ortho, o->outPrefix);
          }
      }

    return 0;
  }

#define mfsb_image_gamma 1.000
  /* Assumed encoding gamma of input and output images. */

#define mfsb_image_bias 0.0327
  /* Assumed encoding bias of input and output images. */

void mfsb_write_image(char *fname, float_image_t *img)
  {
    assert(img->sz[0] == 1);

    float vMin = +INF, vMax = -INF;
    float_image_update_sample_range(img, 0, &vMin, &vMax);
    double vR = fmax(fabs(vMin), fabs(vMax));
    vMin = (float)(vMin < 0.0 ? -vR : 0.0);
    vMax = (float)(vMax > 0.0 ? +vR : 0.0);
    float_image_rescale_samples(img, 0, vMin, vMax, 0.0, 1.0);
    
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    double gamma = mfsb_image_gamma;
    double bias = mfsb_image_bias;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = TRUE;
    float_image_write_pnm_named(fname, img, isMask, gamma, bias, yup, warn, verbose);
  }

void mfsb_show_sample_weights
  ( int32_t NW,           /* Window width and height. */
    double ws[],          /* Sample window weights. */
    char *outPrefix       /* Prefix for outupt file names. */
  )
  {
    fprintf(stderr, "--- testing sample weights NW = %d -------------------\n", NW);

    int32_t NC = 1;
    int32_t NX = NW;
    int32_t NY = NW;
    float_image_t *img = float_image_new(NC, NX, NY);
    
    /* Write weights as an image: */
    for (int32_t ix = 0; ix < NX; ix++)
      { for (int32_t iy = 0; iy < NY; iy++) 
          { int32_t ks = iy*NW + ix;
            fprintf(stderr, "%3d %3d %3d %16.8f\n", ks, ix, iy, ws[ks]);
            float_image_set_sample(img, 0, ix, iy, (float)ws[ks]);
          }
      }
    char *fname = NULL;
    asprintf(&fname, "%s-n%03d-ws.ppm", outPrefix, NW);
    mfsb_write_image(fname, img);
    free(fname);

    float_image_free(img);

    fprintf(stderr, "---------------------------------------------------------------------\n");
  }  

void mfsb_show_single_basis
  ( int32_t NW,           /* Window width and height. */
    double ws[],          /* Sample window weights. */
    multifok_focus_op_basis_type_t bType,          /* Basis type. */
    bool_t ortho,         /* True to orthonormalize the basis. */
    char *outPrefix       /* Prefix for outupt file names. */
  )
  {
    char bTypeX = "NDH"[bType]; 
    fprintf(stderr, "--- showing basis NW = %d bType = %c ortho = %c -------------------\n", NW, bTypeX, "FT"[ortho]);

    int32_t NC = 1;
    int32_t NX = NW;
    int32_t NY = NW;
    float_image_t *img = float_image_new(NC, NX, NY);
    
    /* Get and print the basis for the focus indicator: */
    int32_t NB;
    double **phi = NULL;
    double *wc = NULL;
    char **name = NULL;
    multifok_focus_op_basis_make(NW, ws, bType, ortho, &NB, &phi, &wc, &name);

    multifok_focus_op_basis_print(stderr, NW, NB, phi, wc, name);
    multifok_focus_op_basis_ortho_check(stderr, NW, ws, NB, phi);

    /* Write each basis element as an image: */
    for (int32_t kb = 0; kb < NB; kb++)
      { for (int32_t ix = 0; ix < NX; ix++)
          { for (int32_t iy = 0; iy < NY; iy++) 
              { int32_t ks = iy*NW + ix;
                double wsk = ws[ks];
                double phk = phi[kb][ks];
                double vsk = phk*sqrt(wsk);
                float_image_set_sample(img, 0, ix, iy, (float)vsk);
              }
          }
        char *fname = NULL;
        asprintf(&fname, "%s-n%03d-b%c-o%c-%03d.ppm", outPrefix, NW, bTypeX, "FT"[ortho], kb);
        mfsb_write_image(fname, img);
        free(fname);
      }

    float_image_free(img);
    multifok_focus_op_basis_free(NB, phi, wc, name);

    fprintf(stderr, "---------------------------------------------------------------------\n");
  }  

mfsb_options_t *mfsb_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, "HELP!");
    argparser_set_info(pp, "KNOW THYSELF");
    argparser_process_help_info_options(pp);
    
    mfsb_options_t *o = notnull(malloc(sizeof(mfsb_options_t)), "no mem");

    argparser_get_keyword(pp, "-winSize");
    o->winSize = (int32_t)argparser_get_next_int(pp, 3, 99);  
    if (o->winSize % 2 != 1) { argparser_error(pp, "window size must be odd"); }

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
