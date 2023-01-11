#define PROG_NAME "mf_0200_op_correl"
#define PROG_DESC "Generates data for regression between actual focus and basis of {multifok_focus_op.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-01-09 14:05:52 by stolfi */ 
/* Created on 2012-01-25 by J. Stolfi, UNICAMP */

#define mf_0200_op_correl_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "Duh?"

#define PROG_INFO \
  "SYNOPSIS" \
  "  Writes to \"out/{outPrefix}{tag}-b{basisType}-vals.txt\" a file with columns\n" \
  "\n" \
  "  {ix} {iy} {rFoc} {rFoc^2} {b[0]} .. {b[NB-1]} \n" \
  "\n" \
  " where {ix} and {iy} are the column and row of the pixel, {rFoc} is the nominal" \
  " focal blur radius at that pixel (as read from the focal blur radius image), and" \
  " {b[0..NB-1]} as the coefficients of the window samples in the basis elements." \
  " The basis coefficients are normalized so as to be independent of local image" \
  " brightness and contrast.\n" \
  "\n" \
  "  Also writes to \"out/{outPrefix}{tag}-b{basisType}-names.txt\" the names of " \
  " the basis elements."

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
#include <float_image_read_pnm.h>
#include <float_image_write_pnm.h>
#include <float_image_map_channels.h>
#include <float_image.h>

#include <multifok_focus_op.h>

typedef struct mfoc_options_t 
  { char *inPrefix;    /* Prefix for input image file names. */
    double zFocus;     /* Focus plane {Z} coordinate, for input file names. */
    int32_t winSize;   /* Size of window. */
    multifok_focus_op_basis_type_t basisType;   /* Basis type. */
    double noise;      /* Assumed noise level in image samples. */
    char *outPrefix;   /* Prefix for output filenames. */
  } mfoc_options_t;
  /* Command line parameters. */

int32_t main(int32_t argc, char **argv);

mfoc_options_t *mfoc_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

float_image_t *mfoc_read_color_image(char *inPrefix, char *tag);
  /* Reads an image from file "{inPrefix}{tag}-c.ppm". */ 

float_image_t *mfoc_read_radius_image(char *inPrefix, char *tag);
  /* Reads an image from file "{inPrefix}{tag}-r.pgm". */ 

void mfoc_write_gray_image(float_image_t *img, char *outPrefix, char *tag, char bTypeX, char *tag2, int32_t kb);
  /* Writes the blur radius image {rimg} to file "{outPrefix}{tag}-b{bTypeX}-{tag2}{KKK}.pgm"
    where {KKK} is {kb} zero-padded to 3 digits. If {kb} is negative, the "{KKK}" is suppressed */

void mfoc_write_basis_elem_names(char *outPrefix, char *tag, char bTypeX, int32_t NB, char *name[]);
  /* Writes the basis element names {name[0..NB-1]} to a file file "{outPrefix}{tag}-b{bTypeX}-names.txt"
    one per line. */

FILE *mfoc_open_regression_data_file(char *outPrefix, char *tag, char bTypeX);
  /* Returns the open handle of a file called "{outPrefix}{tag}-b{bTypeX}-vals.txt" for writing. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfoc_options_t *o = mfoc_parse_options(argc, argv);
    
    /* Read the color images: */
    int32_t NC, NX, NY;
    char *tag = NULL;
    asprintf(&tag, "-z%08.4f", o->zFocus);

    float_image_t *cimg = mfoc_read_color_image(o->inPrefix, tag);
    float_image_get_size(cimg, &NC, &NX, &NY);
    assert(NC == 3);
    
    /* Convert to grayscale: */
    float_image_t *gimg = float_image_new(1, NX, NY);
    float_image_map_channels_RGB_to_YUV(cimg, gimg);

    /* Read the focal blur radius image: */
    float_image_t *rimg = mfoc_read_radius_image(o->inPrefix, tag);
    float_image_check_size(rimg, 1, NX, NY);
    
    /* Get the window dimnsions: */
    int32_t NW = o->winSize;
    demand((NW % 2 == 1) && (NW >= 3), "invalid window size");
    int32_t NS = NW*NW; /* Number of samples in window. */

    /* Generate the window sample weights {ws}: */
    double *ws = multifok_focus_op_sample_weights(NW);

    /* Generate the focus op basis {phi[0..NB-1]} and coeff weights {wc[0..NB-1]}: */
    char bTypeX = "NDH"[o->basisType];
    int32_t NB;
    double **phi = NULL;
    double *wc = NULL;
    char **name = NULL;
    bool_t ortho = TRUE;
    multifok_focus_op_basis_make(NW, ws, o->basisType, ortho, &NB, &phi, &wc, &name);
    fprintf(stderr, "obtained NB = %d independent basis elements\n", NB);
    multifok_focus_op_basis_print(stderr, NW, NB, phi, wc, name);
    
    /* Allocate the images with basis coeffs: */
    fprintf(stderr, "allocating images for basis coefs...\n");
    float_image_t *bimg[NB];
    for (int32_t kb = 0; kb < NB; kb++) 
      { bimg[kb] = float_image_new(1, NX, NY); }
      
    /* File with basis element names: */
    mfoc_write_basis_elem_names(o->outPrefix, tag, bTypeX, NB, name);

    /* Allocate composite score image: */
    float_image_t *simg = float_image_new(1, NX, NY);
      
    /* File with nominal blur radius and basis element coeffs for regression: */
    FILE *wr = mfoc_open_regression_data_file(o->outPrefix, tag, bTypeX);

    /* Enumerate all windows in the image and compute the terms for regression: */
    int32_t HW = NW/2;
    float fsmp[NS];
    double dsmp[NS];
    double coef[NB]; /* Coefficients of normalized window in be basis. */
    for (int32_t ix = HW; ix < NX-HW; ix++)
      { for (int32_t iy = HW; iy < NY-HW; iy++) 
          { /* Get the nominal focal blur radius at this pixel: */
            float rFoc = float_image_get_sample(rimg, 0, ix, iy);
            fprintf(wr, "%5d %5d %12.6f %12.6f ", ix, iy, rFoc, rFoc*rFoc);
            
            /* Get the samples in the window and normalize them for brightness and contrast: */
            float_image_get_window_samples(gimg, 0,ix,iy, NW, NW, FALSE, fsmp);
            for (int32_t ks = 0; ks < NS; ks++) { dsmp[ks] = fsmp[ks]; }
            multifok_focus_op_normalize_samples(NW, dsmp, ws, o->noise);
            
            /* Compute and store the coefficients of the basis: */
            for (int32_t kb = 0; kb < NB; kb++)
              { coef[kb] = multifok_focus_op_prod(NW, dsmp, phi[kb], ws);
                fprintf(wr, " %16.12f", coef[kb]);
                float_image_set_sample(bimg[kb], 0, ix, iy, (float)coef[kb]);
              }
              
            /* Compute and store the composite score: */
            for (int32_t ks = 0; ks < NS; ks++) { dsmp[ks] = fsmp[ks]; }
            double score = multifok_focus_op_score_from_basis(NW, dsmp, o->noise, ws, NB, phi, wc);
            float_image_set_sample(simg, 0, ix, iy, (float)score);

            fprintf(wr, "\n");
          }
      }
      fclose(wr);
      
      /* Write out the basis images: */
      for (int32_t kb = 0; kb < NB; kb++)
        { mfoc_write_gray_image(bimg[kb], o->outPrefix, tag, bTypeX, "c", kb); }
      
      /* Write the composite score image: */
       mfoc_write_gray_image(simg, o->outPrefix, tag, bTypeX, "s", -1);
        
      return 0;
  }

#define mfoc_image_gamma 1.000
  /* Assumed encoding gamma of input and output images. */

#define mfoc_image_bias 0.0327
  /* Assumed encoding bisa of input and output images. */

float_image_t *mfoc_read_color_image(char *inPrefix, char* tag)
  {
    char *fname = NULL;
    asprintf(&fname, "%s%s-c.ppm", inPrefix, tag);
    bool_t isMask = FALSE;
    double gamma = 1.000;
    double bias = 0.000;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = FALSE;
    float_image_t *img = float_image_read_pnm_named(fname, isMask, gamma, bias, yup, warn, verbose);
    free(fname);
    return img;
  }

float_image_t *mfoc_read_radius_image(char *inPrefix, char* tag)
  {
    char *fname = NULL;
    asprintf(&fname, "%s%s-r.pgm", inPrefix, tag);
    bool_t isMask = FALSE;
    double gamma = 1.000;
    double bias = 0.000;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = FALSE;
    float_image_t *img = float_image_read_pnm_named(fname, isMask, gamma, bias, yup, warn, verbose);
    free(fname);
    return img;
  }

void mfoc_write_gray_image(float_image_t *img, char *outPrefix, char *tag, char bTypeX, char *tag2, int32_t kb)
  {
    assert(img->sz[0] == 1);

    /* Remap image values to {[0_1]}: */
    float vMin = +INF;
    float vMax = -INF;
    float_image_update_sample_range(img, 0, &vMin, &vMax); /* Ignores infinities. */
    if (vMin > vMax) { vMin = 0.0; vMax = 1.0; }
    float vR = (float)fmax(fabs(vMin), fabs(vMax));
    vMax = vR;
    vMin = (vMin < 0.0f ? -vR : 0.0f);
    float_image_rescale_samples(img, 0, vMin, vMax, 0.0, 1.0);

    char *fname = NULL;
    if (kb >= 0)
      { asprintf(&fname, "%s%s-b%c-%s%03d-r.pgm", outPrefix, tag, bTypeX, tag2, kb); }
    else
      { asprintf(&fname, "%s%s-b%c-%s-r.pgm", outPrefix, tag, bTypeX, tag2); }
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    double gamma = mfoc_image_gamma;
    double bias = mfoc_image_bias;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = TRUE;
    float_image_write_pnm_named(fname, img, isMask, gamma, bias, yup, warn, verbose);
    free(fname);
  }
    
void mfoc_write_basis_elem_names(char *outPrefix, char *tag, char bTypeX, int32_t NB, char *name[])  
  { char *fname = NULL;
    asprintf(&fname, "%s%s-b%c-names.txt", outPrefix, tag, bTypeX);
    FILE *wr = open_write(fname, TRUE);
    for (int32_t kb = 0; kb < NB; kb++)
      { fprintf(wr, "%s\n", name[kb]); }
    fclose(wr);
  }
    
FILE *mfoc_open_regression_data_file(char *outPrefix, char *tag, char bTypeX)
  { char *fname = NULL;
    asprintf(&fname, "%s%s-b%c-vals.txt", outPrefix, tag, bTypeX);
    FILE *wr = open_write(fname, TRUE);
    return wr;
  }

mfoc_options_t *mfoc_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    mfoc_options_t *o = notnull(malloc(sizeof(mfoc_options_t)), "no mem");

    argparser_get_keyword(pp, "-inPrefix");
    o->inPrefix = argparser_get_next(pp);

    argparser_get_keyword(pp, "-winSize");
    o->winSize = (int32_t)argparser_get_next_int(pp, 3, 99);  
    if (o->winSize % 2 != 1) { argparser_error(pp, "window size must be odd"); }

    argparser_get_keyword(pp, "-basisType");
    char *bType = argparser_get_next_non_keyword(pp);
    if (strcmp(bType, "D") == 0) 
      { o->basisType = multifok_focus_op_basis_type_DIFF; }
    else if (strcmp(bType, "H") == 0)
      { o->basisType = multifok_focus_op_basis_type_HART; }
    else if (strcmp(bType, "N") == 0)
      { o->basisType = multifok_focus_op_basis_type_NONE; }
    else 
      { argparser_error(pp, "invalid basis type"); }

    argparser_get_keyword(pp, "-noise");
    o->noise = argparser_get_next_double(pp, 0.0, 100.0);  

    argparser_get_keyword(pp, "-zFocus");
    o->zFocus = argparser_get_next_double(pp, 0.0, 1.0e200);  

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
