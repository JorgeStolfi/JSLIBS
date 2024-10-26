#define PROG_NAME "test_mfok_show_basis"
#define PROG_DESC "test an analysis of basis functions from {multifok_focus_op.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-10-22 13:08:38 by stolfi */ 
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
  "    For each window weight distribution type {wType} (\"BIN\", \"GLD\", etc.), the" \
  " program writes an {NW} by {NW} grayscale image showing the weight values.  The" \
  " file name will be \"{outDir}/weights-nw{NNN}-wt{wType}.png\", where" \
  " {NNN} is the window size {NW}, zero-padded to 3 digits.\n" \
  "\n" \
  "  BASIS IMAGE FILES\n" \
  "     For each each basis type {bType} (\"LAPL\", \"DIFF\", etc.), each window weight" \
  " distribution type {wType},  with and without orthogonalization, and each" \
  " basis element {kb}, the program writes  agrayscale image showing the {NW} by {NW} values" \
  " of that element. The file name is \"{outDir}/basis-nw{NNN}-bt{bType}-wt{wType}-or{ortho}-{KKK}.pgm\" where" \
  " {NNN} is the window size {NW}, {ortho}" \
  " is \"F\" for the raw version of the basis and \"T\" for the orthonormalized" \
  " version, and {KKK} is the element index {kb}. Both {NNN} and {KKK} are" \
  " zero-padded to 3 digits.\n" \
  "\n" \
  "  BASIS NAMES FILE\n" \
  "    Also writes to \"{outDir}/bas-nw{NNN}-bt{bType}-belnames.txt\" the names of " \
  " the basis elements, like \"FXY\", one per line."

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
#include <multifok_image.h>
#include <multifok_basis.h>
#include <multifok_term.h>
#include <multifok_score.h>
#include <multifok_test.h>

typedef struct mfsb_options_t 
  { int32_t winSize;   /* Size of window. */
    char *outDir;   /* Output image file name. */
  } mfsb_options_t;
  /* Command line parameters. */

int32_t main(int32_t argn, char **argv);

mfsb_options_t *mfsb_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

void mfsb_show_sample_weights
  ( multifok_window_type_t wType,  /* Window weights distribution type. */
    int32_t NW,                    /* Window width and height. */
    double ws[],                   /* Sample window weights. */
    char *outDir                   /* Output directory. */
  );
  /* Writes the sample weights {ws[0..NS-1]} as an image file 
    "{outDir}/weights.png" where 
    {NNN} is {NW} zero-padded to 3 digits. */

void mfsb_show_single_basis
  ( int32_t NW,                    /* Window width and height. */
    double ws[],                   /* Sample window weights. */
    multifok_basis_type_t bType,   /* Basis type. */
    multifok_window_type_t wType,  /* Window weights distribution type. */
    bool_t ortho,                  /* True to ortho normalize basis. */
    char *outDir                /* Prefix for outupt file names. */
  );
  /* Writes each basis element {bas[k]} as an image file "{outDir}/basis-nw{NNN}-bt{bType}-wt{wType}-or{ortho}-{KKK}.ppm" where 
    {NNN} and {KKK} are {NW} and {k} zero-padded to 3 digits.  Each image pixel is a basis value {bas[kb][ks]}
    times {sqrt(ws[ks])}. */

/* IMPLEMENTATIONS */

#define CANC multifok_basis_type_CANC  
#define LAPL multifok_basis_type_LAPL  
#define DIFF multifok_basis_type_DIFF  
#define HART multifok_basis_type_HART  

#define BT_FIRST multifok_basis_type_FIRST
#define BT_LAST multifok_basis_type_LAST

#define WT_FIRST multifok_window_type_FIRST
#define WT_LAST multifok_window_type_LAST

int32_t main (int32_t argc, char **argv)
  {
    mfsb_options_t *o = mfsb_parse_options(argc, argv);
    
    int32_t NW = o->winSize;
    demand((NW % 2 == 1) && (NW >= 3), "invalid window size");

    for (multifok_window_type_t wt = WT_FIRST; wt <= WT_LAST; wt++)
      { double *ws = multifok_window_weights(NW, wt);
        mfsb_show_sample_weights(wt, NW, ws, o->outDir);
        for (multifok_basis_type_t bt = BT_FIRST; bt <= BT_LAST; bt++)
          { if ((NW == 3) || (bt != DIFF))
              { for (int32_t io = 0; io <= 1; io++)
                  { bool_t ortho = (io == 1);
                    char *basDir = NULL;
                    asprintf(&basDir, "%s/basis-nw%03d-bt%s-wt%s-or%c", o->outDir, NW, bTypeX, wTypeX, "FT"[ortho]);
                    mfsb_show_single_basis(NW, ws, bt, wt, ortho, basDir);
                    free(basDir);
                  }
              }
          }
      }

    return 0;
  }

void mfsb_show_sample_weights
  ( multifok_window_type_t wType,  /* Window weights distribution type. */
    int32_t NW,                    /* Window width and height. */
    double ws[],                   /* Sample window weights. */
    char *outDir                /* Prefix for outupt file names. */
  )       
  {
    char *wTypeX = multifok_window_type_to_text(wType); 
    fprintf(stderr, "--- testing sample weights wType = %s NW = %d -------------------\n", wTypeX, NW);

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
      
    char *fileName = NULL;
    asprintf(&fileName, "%s/weights-wt%s-nw%03d", outDir, wTypeX, NW);
    multifok_image_sample_weights_write(wsimg, fileName);
    free(fileName);
  
    float_image_free(wsimg);

    fprintf(stderr, "---------------------------------------------------------------------\n");
  }  

void mfsb_show_single_basis
  ( int32_t NW,                    /* Window width and height. */
    double ws[],                   /* Sample window weights. */
    multifok_basis_type_t bType,   /* Basis type. */
    multifok_window_type_t wType,  /* Basis type. */
    bool_t ortho,                  /* True to orthonormalize the basis. */
    char *basDir                   /* Prefix for outupt file names. */
  )
  {
    char *bTypeX = multifok_basis_type_to_text(bType); 
    char *wTypeX = multifok_window_type_to_text(wType); 
    fprintf(stderr, "--- showing basis NW = %d bType = %s wType = %s ortho = %c -------------------\n", NW, bTypeX, wTypeX, "FT"[ortho]);

    int32_t NC = 1;
    int32_t NX = NW;
    int32_t NY = NW;
    float_image_t *bcimg[NB];
    
    /* Get and print the basis for the focus indicator: */
    multifok_basis_t *basis = multifok_basis_make(bType, NW, ws, ortho);

    multifok_basis_print(stderr, basis);
    multifok_basis_ortho_check(stderr, basis);

    /* Write each basis element as an image: */
    for (int32_t kb = 0; kb < basis->NB; kb++)
      { bcimg[kb] = float_image_new(NC, NX, NY);
        for (int32_t ix = 0; ix < NX; ix++)
          { for (int32_t iy = 0; iy < NY; iy++) 
              { int32_t ks = iy*NW + ix;
                double phk = basis->bas[kb][ks];
                float_image_set_sample(bcimg[kb], 0, ix, iy, (float)phk);
              }
          }
      }
    multifok_image_basis_kernels_write(NB, bcimg, basDir);

    float_image_free(bcimg);
    multifok_basis_free(basis);

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

    argparser_get_keyword(pp, "-outDir");
    o->outDir = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
