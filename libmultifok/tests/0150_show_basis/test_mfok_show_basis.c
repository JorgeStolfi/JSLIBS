#define PROG_NAME "test_mfok_show_basis"
#define PROG_DESC "test an analysis of basis functions from {multifok_focus_op.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-04-28 17:30:11 by stolfi */ 
/* Created on 2023-01-05 by J. Stolfi, UNICAMP */

#define test_mfok_show_basis_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "Duh?"

#define PROG_INFO \
  "SYNOPSIS" \
  "  Writes images showing the values of basis elements in the local {NW} by {NW} window\n" \
  "OUTPUTS\n" \
  "\n" \
  "  WINDOW WEIGHT IMAGE FILE\n" \
  "    The program writes a grayscale image" \
  " showing the {NW} by {NW} values of basis element {kb}.\n" \
  "\n" \
  "    The image is written to file \"{outPrefix}-n{NNN}-ws.pgm\" where" \
  " {NNN} is the window size {NW}, zero-padded to 3 digits.\n" \
  "\n" \
  "  BASIS IMAGE FILES\n" \
  "    For each basis and each basis element index {kb}, the program writes grayscale image" \
  " showing the {NW} by {NW} values of basis element {kb}.\n" \
  "\n" \
  "    The image is written to file \"{outPrefix}-n{NNN}-b{bType}-o{ortho}-{KKK}.pgm\" where" \
  " {NNN} is the window size {NW}, {bType} is the basis type (\"LAPL\", \"DIFF\", etc.) {ortho}" \
  " is 0 for the raw version of the basis and 1 for the orthonormalized one, and {KKK} is" \
  " the element index {kb}. Both {NNN} and {KKK} are zero-padded to 3 digits.\n" \
  "\n" \
  "  BASIS NAMES FILE\n" \
  "    Also writes to \"{outPrefix}-n{NNN}-b{bType}-belnames.txt\" the names of " \
  " the basis elements, like \"DXDY\", one per line."

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

#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_term.h>
#include <multifok_score.h>
#include <multifok_test.h>

typedef struct mfsb_options_t 
  { int32_t winSize;   /* Size of window. */
    char *outPrefix;   /* Output image file name. */
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
    multifok_basis_type_t bType,  /* Basis type. */
    bool_t ortho,                          /* True to ortho normalize basis. */
    char *outPrefix                        /* Prefix for outupt file names. */
  );
  /* Writes each basis element {bas[k]} as an image file "{outPrefix}-n{NNN}-b{bType}-o{ortho}-{KKK}.ppm" where 
    {NNN} and {KKK} are {NW} and {k} zero-padded to 3 digits.  Each image pixel is a basis value {bas[kb][ks]}
    times {sqrt(ws[ks])}. */

/* IMPLEMENTATIONS */

#define CANC multifok_basis_type_CANC  
#define LAPL multifok_basis_type_LAPL  
#define DIFF multifok_basis_type_DIFF  
#define HART multifok_basis_type_HART  
#define FIRST multifok_basis_type_FIRST
#define LAST multifok_basis_type_LAST

int32_t main (int32_t argc, char **argv)
  {
    mfsb_options_t *o = mfsb_parse_options(argc, argv);
    
    int32_t NW = o->winSize;
    demand((NW % 2 == 1) && (NW >= 3), "invalid window size");

    double *ws = multifok_window_sample_weights(NW);
    mfsb_show_sample_weights(NW, ws, o->outPrefix);

    bool_t ortho;
    for (multifok_basis_type_t bt = FIRST; bt <= LAST; bt++)
      { if ((NW == 3) || (bt != DIFF))
          { ortho = FALSE;
            mfsb_show_single_basis(NW, ws, bt, ortho, o->outPrefix);
            ortho = TRUE;
            mfsb_show_single_basis(NW, ws, bt, ortho, o->outPrefix);
          }
      }

    return 0;
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
    float_image_t *wsimg = float_image_new(NC, NX, NY);
    
    /* Write weights as an image: */
    for (int32_t ix = 0; ix < NW; ix++)
      { for (int32_t iy = 0; iy < NW; iy++) 
          { int32_t ks = iy*NW + ix;
            float_image_set_sample(wsimg, 0, ix, iy, (float)ws[ks]);
            double ws_check = (double)float_image_get_sample(wsimg, 0, ix, iy);
            fprintf(stderr, "%3d %3d %3d %16.12f = %12.7f\n", ks, ix, iy, ws[ks], ws_check);
          }
      }
      
    char *tag = NULL;
    asprintf(&tag, "-n%03d", NW);
    multifok_test_write_sample_weights_image(wsimg, outPrefix, tag);
    free(tag);
  
    float_image_free(wsimg);

    fprintf(stderr, "---------------------------------------------------------------------\n");
  }  

void mfsb_show_single_basis
  ( int32_t NW,           /* Window width and height. */
    double ws[],          /* Sample window weights. */
    multifok_basis_type_t bType,   /* Basis type. */
    bool_t ortho,         /* True to orthonormalize the basis. */
    char *outPrefix       /* Prefix for outupt file names. */
  )
  {
    char *bTypeX = multifok_basis_type_to_text(bType); 
    fprintf(stderr, "--- showing basis NW = %d bType = %s ortho = %c -------------------\n", NW, bTypeX, "FT"[ortho]);

    int32_t NC = 1;
    int32_t NX = NW;
    int32_t NY = NW;
    float_image_t *bcimg = float_image_new(NC, NX, NY);
    
    /* Get and print the basis for the focus indicator: */
    int32_t NB;
    double **bas = NULL;
    char **belName = NULL;
    multifok_basis_make(bType, NW, ws, ortho, &NB, &bas, &belName);

    multifok_basis_print(stderr, NW, NB, bas, belName);
    multifok_basis_ortho_check(stderr, NW, NB, bas);

    /* Write each basis element as an image: */
    for (int32_t kb = 0; kb < NB; kb++)
      { for (int32_t ix = 0; ix < NX; ix++)
          { for (int32_t iy = 0; iy < NY; iy++) 
              { int32_t ks = iy*NW + ix;
                double phk = bas[kb][ks];
                float_image_set_sample(bcimg, 0, ix, iy, (float)phk);
              }
          }
        char *tag = NULL;
        asprintf(&tag, "-n%03d-%s-o%c-%03d", NW, bTypeX, "FT"[ortho], kb);
        multifok_test_write_basis_elem_image(bcimg, outPrefix, tag);
        free(tag);
      }
 
    float_image_free(bcimg);
    multifok_basis_free(NB, bas, belName);

    fprintf(stderr, "---------------------------------------------------------------------\n");
  }  

mfsb_options_t *mfsb_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
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
