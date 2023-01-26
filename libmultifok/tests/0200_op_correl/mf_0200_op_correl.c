#define PROG_NAME "mf_0200_op_correl"
#define PROG_DESC "Generates data for regression between actual sharpness and basis of {multifok_focus_op.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-01-25 05:32:12 by stolfi */ 
/* Created on 2023-01-07 by J. Stolfi, UNICAMP */

#define mf_0200_op_correl_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "Duh?"

#define PROG_INFO \
  "SYNOPSIS" \
  "\n" \
  "  The terms {term[0..NT-1]} are derived from" \
  " the coefficients of the window samples projected on a specified local operator basis, from a small" \
  " window surrounding of the pixel, ." \
  " The basis coefficients are computed from normalized window samples so as to" \
  " be independent of local image brightness and contrast.\n" \
  "INPUTS" \
  "  Reads a collection of image quadruples {cimg[ki]}, {zimg[ki]}, {dimg[ki]}, {simg[ki]} where\n" \
  "\n" \
  "   {cimg[ki]} is sythetic image of a 3D scene  with simulated focus blurring.\n" \
  "   {zimg[ki]} specifies the average {Z}-coord of the scene within each pixel.\n" \
  "   {dimg[ki]} specifies the deviation of that {Z} coordinate within each pixel.\n" \
  "   {simg[ki]} specifies the sharpness of {cimg[k]} at each pixel.\n" \
  "\n" \
  "  The {cimg} images are in color, the other images are greyscale with linear" \
  " encoding (gamma = 1).  They may have different sizes and focus plane" \
  " positions {zFoc}, but must all have the same depth of focus {zDep}.\n" \
  "\n" \
  "OUTPUTS\n" \
  "\n" \
  "  REGRESSION DATA FILE\n" \
  "    Writes to \"{outPrefix}-rdata.txt\" a file with columns\n" \
  "\n" \
  "    P{ki}.{ix}.{iy} {sharp^exp} {term[0]} .. {term[NT-1]} \n" \
  "\n" \
  " where {ki} is the image index, {ix} and {iy} are the column and row" \
  " of the pixel, {sharp} is the \"actual\" sharpness at that" \
  " pixel (as read from the {simg[ki]} image), {wp} is a pixel importance weight for correlation, and" \
  " {term[0..NT-1]} are terms that are to be used to estimate the sharpness, as" \
  " per {multifok_focus_op_score_from_basis_coeffs}.  The exponent {exp} is either" \
  " 1 or 2, depending on the \"-squared\" option.\n" \
  "\n" \
  "  SCORE COMPARISON FILE\n" \
  "\n" \
  "    Also writes to \"{outPrefix}-fcomp.txt\" a file with columns\n" \
  "\n" \
  "    P{ki}.{ix}.{iy} {sharp} {score} {zave} {zdev} \n" \
  "\n" \
  " where {ki,ix,iy,sharp} are as above, {zave} is the average scene height in the" \
  " pixel (as read from the {zimg[ki]} image), compute {score}, {score} is" \
  " the sharpness estimated by {multifok_focus_op_score_from_basis_coeffs}\n" \
  "\n" \
  "    Optionally, the fields {sharp} and {score} are squared before being written to the regression file.\n" \
  "\n" \
  "  BASIS NAMES FILE\n" \
  "\n" \
  "    Also writes to \"{outPrefix}-belnames.txt\" the names of " \
  " the basis elements, like \"DXDY\", one per line.\n" \
  "\n" \
  "  TERM NAMES FILE\n" \
  "\n" \
  "    Also writes to \"{outPrefix}-terms.txt\" the names and weights of " \
  " the terms used to compute the score, like \"DX*DY +0.531\", one per line."

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <i2.h>
#include <jsfile.h>
#include <argparser.h>
#include <bool.h>
#include <float_image_read_pnm.h>
#include <float_image_write_pnm.h>
#include <float_image_map_channels.h>
#include <float_image.h>

#include <multifok_focus_op.h>
#include <multifok_test.h>
#include <multifok_scene.h>

typedef struct mfoc_options_t 
  { string_vec_t inPrefix;                     /* Prefixes for input image file names. */
    int32_t winSize;                           /* Size of window. */
    multifok_focus_op_basis_type_t basisType;  /* Basis type. */
    char *termsFile;                           /* Name of file with term names and weights. */
    double noise;                              /* Assumed noise level in image samples. */
    bool_t squared;                            /* If true does regression on the square of sharpness. */ 
    char *outPrefix;                           /* Prefix for output filenames. */
  } mfoc_options_t;
  /* Command line parameters. */

int32_t main(int32_t argc, char **argv);

void mfoc_compute_basis
  ( int32_t NW,
    double ws[],
    multifok_focus_op_basis_type_t bType,
    int32_t *NB_P,
    double ***bas_P,
    char ***belName_P,
    char *outPrefix
  );
  /* Creates an operator basis of the specified {basisType} for an {NW} by {NW} window with
    with window sample weights {ws[0..NS-1]}, where {NS=NW*NW}.  See {multifok_focus_op_basis_make}. 
    The basis will be orrthonormal.  The number {NB} of basis elements is returned in {*NB_P},
    the basis elements {bas[0..NB-1][0..NS-1]}  are returned in {*bas_P}, and their names
    {belName[0..NB-1]} are returned in {*belName_P}.
    
    Also prints the basis and writes its names out to file with
    {mfoc_write_basis_elem_names(outPrefix,NP,belName)}. */

void mfoc_process_image_set
  ( FILE *wr_reg,
    FILE *wr_cmp,
    int32_t ki, 
    char *inPrefix,
    int32_t NW,
    double ws[],
    double noise,
    int32_t NB,
    double **bas,
    int32_t NT,
    double wt[],
    i2_t tix[],
    bool_t squared,
    char *outPrefix
  );
  /* Processes the image set with index {ki} which is read from files "{inPrefix}-{tail}" where
    {tail} is "-c.ppm" (simulated view), "-z.pgm" ({Z}-coordinates), "-d.pgm" ({Z}-deviation),
    and "-s.pgm" ("actual" sharpness {sharp}). 
    
    Writes images "{outPrefix}-{KKK}-b.pgm" with the basis 
    coeff values at each pixel, where {KKK} is the basis element index 
    in {0..NB-1} zero-padded to 3 digits.
    
    Also writes images {fimg} = "{outPrefix}-f.pgm" with the sharpness {score} computed
    from the basis coeffs, and {eimg} = "{outPrefix}-e.pgm" with the error {score-sharp}.
    See {multifok_focus_op_score_from_basis_coeffs} for how the {score} is computed.
    
    Writes out to {wr_reg} and {wr_cmp} one line for each pixel considered
    useful for the regression, as explained in {PROG_INFO}.
    
    If {squared} is true, assumes that the combination of the terms is
    an estimate of the squared sharpness, so that the {score} is the
    square root of that combination. See
    {multifok_focus_op_score_from_basis_coeffs}. In that case, writes {sharp^2}
    instead of {sharp} to {wr_reg}. */

mfoc_options_t *mfoc_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

FILE *mfoc_open_text_file(char *outPrefix, char *tag);
  /* Returns the open handle of a file called "{outPrefix}{tag}.txt" for writing. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfoc_options_t *o = mfoc_parse_options(argc, argv);
    
    /* Read the images: */
    int32_t NI = o->inPrefix.ne; /* Number of input images. */
   
    /* Get the window dimnsions: */
    int32_t NW = o->winSize;
    demand((NW % 2 == 1) && (NW >= 3), "invalid window size");

    /* Generate the window sample weights {ws}: */
    double *ws = multifok_focus_op_sample_weights(NW);

    int32_t NB;
    double **bas = NULL;
    char **belName = NULL;
    mfoc_compute_basis(NW, ws, o->basisType, &NB, &bas, &belName, o->outPrefix);

    /* Get term indices and coefficients for computed sharpness score, if any: */
    int32_t NT;   /* Number of terms in the score formula. */
    char **tname; /* Term names. */
    double *wt;   /* Term weights. */
    i2_t *tix;    /* Basis element indices of terms. */
    multifok_test_read_term_names_and_weights(o->termsFile, NB, belName, &NT, &tname, &wt, &tix, TRUE);
    multifok_test_write_term_names_and_weights(o->outPrefix, NT, tname, wt);

    /* File with nominal blur radius and basis element coeffs for regression: */
    FILE *wr_reg = mfoc_open_text_file(o->outPrefix, "-rdata");
    FILE *wr_cmp = mfoc_open_text_file(o->outPrefix, "-cdata");

    for (int32_t ki = 0; ki < NI; ki++)
      { mfoc_process_image_set
          ( wr_reg, wr_cmp, 
            ki, o->inPrefix.e[ki],
            NW, ws, o->noise, 
            NB, bas, 
            NT, wt, tix, o->squared,
            o->outPrefix
          );
      } 
      
    fclose(wr_reg);
    fclose(wr_cmp);
      
    return 0;
  }
    
void mfoc_compute_basis
  ( int32_t NW,
    double ws[],
    multifok_focus_op_basis_type_t bType,
    int32_t *NB_P,
    double ***bas_P,
    char ***belName_P,
    char *outPrefix
  )
  {
    /* Generate the focus op basis {bas[0..NB-1]}: */
    int32_t NB;
    double **bas = NULL;
    char **belName = NULL;
    bool_t ortho = TRUE;
    multifok_focus_op_basis_make(NW, ws, bType, ortho, &NB, &bas, &belName);
    fprintf(stderr, "obtained NB = %d independent basis elements\n", NB);
    multifok_focus_op_basis_print(stderr, NW, NB, bas, belName);
    multifok_focus_op_basis_ortho_check(stderr, NW, ws, NB, bas);
      
    /* File with basis element names: */
    multifok_test_write_basis_elem_names(outPrefix, NB, belName);
      
    (*NB_P) = NB;
    (*bas_P) = bas;
    (*belName_P) = belName;
  }

void mfoc_process_image_set
  ( FILE *wr_reg,
    FILE *wr_cmp,
    int32_t ki, 
    char *inPrefix,
    int32_t NW,
    double ws[],
    double noise,
    int32_t NB,
    double **bas,
    int32_t NT,
    double wt[],
    i2_t tix[],
    bool_t squared,
    char *outPrefix
  )
  {
    int32_t NS = NW*NW; /* Number of samples in window. */

    /* Read the color image: */
    int32_t NC, NX, NY;
    float_image_t *cimg = multifok_test_read_scene_color_image(inPrefix, "");
    float_image_get_size(cimg, &NC, &NX, &NY);
    assert(NC == 3);
    
    /* Convert to grayscale: */
    float_image_t *gimg = float_image_new(1, NX, NY);
    float_image_map_channels_RGB_to_YUV(cimg, gimg);

    /* Read the actual sharpness image: */
    float_image_t *simg = multifok_test_read_sharpness_image(inPrefix, "");
    float_image_check_size(simg, 1, NX, NY);
    
    /* Read the scene {Z} average image: */
    float_image_t *zimg = multifok_test_read_zave_image(inPrefix, "");
    float_image_check_size(zimg, 1, NX, NY);
     
    /* Read the scene {Z} deviation image: */
    float_image_t *dimg = multifok_test_read_zdev_image(inPrefix, "");
    float_image_check_size(dimg, 1, NX, NY);
    
    /* Allocate the images with basis coeffs: */
    fprintf(stderr, "allocating images for basis coefs...\n");
    float_image_t *bimg[NB];
    for (int32_t kb = 0; kb < NB; kb++) 
      { bimg[kb] = float_image_new(1, NX, NY); }

    /* Allocate composite score image: */
    float_image_t *fimg = float_image_new(1, NX, NY);
      
    /* Allocate composite score error image: */
    float_image_t *eimg = float_image_new(1, NX, NY);

    /* Enumerate all windows in the image and compute the terms for regression: */
    int32_t HW = NW/2;
    float fsmp[NS];
    double dsmp[NS];
    double coeff[NB]; /* Coefficients of normalized window in be basis. */
    double term[NT]; /* Score terms derived from {coeff[0..NB-1]}. */
    int32_t NP_scan = 0;    /* Number of pixels processed. */
    int32_t NP_edge = 0;    /* Number of pixels discarded because of high {zdev}. */
    int32_t NP_badsco = 0;  /* Number of pixels discarded because of {score} out of {[0_1]}. */
    int32_t NP_reg = 0;     /* Number of pixels considered in regression. */
    for (int32_t ix = HW; ix < NX-HW; ix++)
      { for (int32_t iy = HW; iy < NY-HW; iy++) 
          { /* Get the samples in the window and normalize them for brightness and contrast: */
            float_image_get_window_samples(gimg, 0,ix,iy, NW, NW, FALSE, fsmp);
            for (int32_t ks = 0; ks < NS; ks++) { dsmp[ks] = fsmp[ks]; }
            double ave, dev;
            multifok_focus_op_normalize_window_samples(NW, dsmp, ws, noise, &avg, &dev); 
            
            /* Compute coefficients of window in basis and estimated sharpness score: */
            multifok_focus_op_compute_basis_coeffs(NW, dsmp, ws, NB, bas, coeff);
            double score = multifok_focus_op_score_from_basis_coeffs(NW, NB, coeff, NT, tix, wt, term, squared);
            
            /* Save coeffs in basis coeffs images: */
            for (int32_t kb = 0; kb < NB; kb++)
              { float_image_set_sample(bimg[kb], 0, ix, iy, (float)coeff[kb]); }
              
            /* Save the computed score in its image: */
            float_image_set_sample(fimg, 0, ix, iy, (float)score);
              
            /* Get the "true" sharpness {sharp} at this pixel. */
            double sharp = float_image_get_sample(simg, 0, ix, iy);
            assert((sharp >= 0.0) && (sharp <= 1.0));
            
            /* Get the scene {Z} height average {zave} at this pixel: */
            float zave = float_image_get_sample(zimg, 0, ix, iy);
            
            /* Get the scene {Z} height deviation {zdev} at this pixel: */
            float zdev = float_image_get_sample(dimg, 0, ix, iy);
            
            /* Save the computed score error in its image: */
            double ds = score - sharp;
            float_image_set_sample(eimg, 0, ix, iy, (float)ds);
            
            /* Decide whether the pixel is useful for the regression: */
            double zdev_max = 0.05*multifok_scene_ZMAX; /* Ignore pixels which more than this variance. */
            if (zdev > zdev_max)
              { NP_edge++; }
            else if ((score <= -0.500) || (score >= 1.500))
              { NP_badsco++; }
            else
              { double wp = sharp*sharp;
                /* double wp = 1.0; */

                /* Write the data for regression: */
                fprintf(wr_reg, "P%d.%d.%d ", ki, ix, iy);
                fprintf(wr_reg, " %14.10f ", (squared ? sharp*sharp : sharp)); /* "Actual" sharpness. */
                fprintf(wr_reg, " %12.6f ", wp); /* Pixel weight. */
                for (int32_t kt = 0; kt < NT; kt++) { fprintf(wr_reg, " %16.12f", term[kt]); }
                fprintf(wr_reg, "\n");

                /* Write the data for comparisons: */
                fprintf(wr_cmp, "P%d.%d.%d ", ki, ix, iy);
                fprintf(wr_cmp, " %12.6f ", wp); /* Pixel weight. */
                fprintf(wr_cmp, " %14.10f ", sharp); /* "Actual" sharpness. */
                fprintf(wr_cmp, " %14.10f ", score); /* Computed sharpness score. */
                fprintf(wr_cmp, " %+12.6f ", zave); /* Avg scene height rel to focus plane in pixel. */
                fprintf(wr_cmp, " %12.6f ", zdev); /* Deviation of scene height in pixel. */
                fprintf(wr_cmp, "\n");
                
                NP_reg++;
              }
            NP_scan++;
          }
      }
    fprintf(stderr, "scanned %d pixels\n", NP_scan); 
    fprintf(stderr, "ignored %d pixels for high {zdev}\n", NP_edge); 
    fprintf(stderr, "ignored %d pixels for {score} out of {[0_1]}\n", NP_badsco); 
    fprintf(stderr, "using %d pixels for regression\n", NP_reg); 

    demand(NP_reg > 0, "ABORTED");

    /* Write out the basis images: */
    for (int32_t kb = 0; kb < NB; kb++)
      { char *tag = NULL;
        asprintf(&tag, "-%03d", kb);
        multifok_test_write_basis_coeff_image(bimg[kb], outPrefix, tag);
        free(tag);
      }

    /* Write the computed score image: */
     multifok_test_write_score_image(fimg, outPrefix, "");

    /* Write the computed score error image: */
     multifok_test_write_score_error_image(eimg, outPrefix, "");
  }

FILE *mfoc_open_text_file(char *outPrefix, char *tag)
  { char *fname = NULL;
    asprintf(&fname, "%s%s.txt", outPrefix, tag);
    FILE *wr = open_write(fname, TRUE);
    return wr;
    free(fname);
  }

mfoc_options_t *mfoc_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    mfoc_options_t *o = notnull(malloc(sizeof(mfoc_options_t)), "no mem");
    
    o->inPrefix = string_vec_new(20);
    int32_t NI = 0; /* Number of input image sets. */
    while (argparser_keyword_present(pp, "-inPrefix"))
      { string_vec_expand(&(o->inPrefix), NI);
        o->inPrefix.e[NI] = argparser_get_next(pp);
        NI++;
      }
    string_vec_trim(&(o->inPrefix), NI);
    if (NI == 0) { argparser_error(pp, "must specify at least one \"-inPrefix\""); }

    argparser_get_keyword(pp, "-winSize");
    o->winSize = (int32_t)argparser_get_next_int(pp, 3, 99);  
    if (o->winSize % 2 != 1) { argparser_error(pp, "window size must be odd"); }

    argparser_get_keyword(pp, "-basisName");
    char *bName = argparser_get_next_non_keyword(pp);
    o->basisType = multifok_focus_op_basis_type_from_text(bName, FALSE)
    if (bName < 0) { argparser_error(pp, "invalid basis type"); }
      
    argparser_get_keyword(pp, "-termsFile");
    o->termsFile = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-squared");
    o->squared = argparser_get_next_bool(pp);  

    argparser_get_keyword(pp, "-noise");
    o->noise = argparser_get_next_double(pp, 0.0, 100.0);  

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
