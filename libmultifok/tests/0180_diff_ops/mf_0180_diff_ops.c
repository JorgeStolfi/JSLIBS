#define PROG_NAME "mf_0180_diff_ops"
#define PROG_DESC "Amalyzes dependency of differential operators on sharpness.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-01-25 23:53:52 by stolfi */ 
/* Created on 2023-01-24 by J. Stolfi, UNICAMP */

#define mf_0180_diff_ops_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "Duh?"

#define PROG_INFO \
  "SYNOPSIS" \
  "  Reads {NI} images with focus blur. Computes the coefficients of the {NB} elements" \
  " of a specified local operator basis at every pixel. Writes images of those coefficients.  Writes" \
  " this data as a text file for plot.\n" \
  "\n" \
  "  The basis coefficients are computed from normalized window samples so as to" \
  " be independent of local image brightness and contrast.\n" \
  "INPUTS" \
  "  Reads a collection of image quadruples {cimg[ki]}, {zimg[ki]}, {dimg[ki]}, {simg[ki]} for {ki} in {0..NI-1}, where\n" \
  "\n" \
  "   {cimg[ki]} is sythetic image of a 3D scene  with simulated focus blurring.\n" \
  "   {zimg[ki]} specifies the average {Z}-coord of the scene within each pixel.\n" \
  "   {dimg[ki]} specifies the deviation of that {Z} coordinate within each pixel.\n" \
  "   {simg[ki]} specifies the sharpness of {cimg[k]} at each pixel.\n" \
  "\n" \
  "  The {cimg} images should be in color, but are converted internally to a" \
  " grayscale image {gimg[ki]}. The other images are greyscale with linear" \
  " encoding (gamma = 1).  They may have different sizes and focus plane" \
  " positions {zFoc}, but must all have the same depth of focus {zDep}.\n" \
  "\n" \
  "OUTPUTS\n" \
  "\n" \
  "  BASIS IMAGE FILES\n" \
  "    For each input image {cimg[ki]} and each basis element {kb} in {0..NB-1}, writes" \
  " to \"{outPrefix}-{NNN}-bq.pgm\" a grayscale" \
  " image showing the SQUARED coefficient of that basis element computed on normalized window" \
  " samples; where {NNN} is the element index {kb} zero-padded to 3 digits." \
  "\n" \
  "  LOCAL AVERAGE AND DEVIATION IMAGE FILES\n" \
  "    For each input image {cimg[ki]} the program computes images {avimg[ki]} and {dvimg[ki]} with" \
  " the window sample average and deviation around each pixel of the grayscale image {gimg[ki]}, and writes" \
  " them to \"{outPrefix}{tail}\" where {tail} is \"-av.pgm\" and \"-dv.pgm\", respectively.\n" \
  "\n" \
  "  LOCALLY NORMALIZED IMAGE FILES\n" \
  "    For each input image {cimg[ki]} writes" \
  " to \"{outPrefix}-n.txt\" a grayscale image showing the central window pixel after window" \
  " normalization, mapped from {[-1_+1]} to {[0_1]}.\n" \
  "\n" \
  "  PLOT DATA FILE\n" \
  "    Writes to \"{outPrefix}-odata.txt\" a file with columns\n" \
  "\n" \
  "    P{ki}.{ix}.{iy} {wp} {vave} {vdev} {sharp} {coeff[0]} .. {coeff[NB-1]} \n" \
  "\n" \
  " where {ki} is the image index, {ix} and {iy} are the column and row" \
  " of the pixel, {vave} is the window sample average and deviation, {wp} is a pixel" \
  " relevance weight, {sharp} is the \"actual\" sharpness at that" \
  " pixel (as read from the {simg[ki]} image), and" \
  " {coeff[0..NB-1]} are the coefficients of the normalized window samples" \
  " in the basis elements {bas[0..NB-1]}.\n" \
  "\n" \
  "  The values of {vave} and {vdev} are computed from the window of {gimg[hi]} centered" \
  " at the pixel, before they are normalized, taking the window sample weights" \
  " into account.  The local differential operators are computed after the" \
  " window sanmples have been normalized for brightnes and contrast.\n" \
  "\n" \
  "  BASIS NAMES FILE\n" \
  "    Also writes to \"{outPrefix}-belnames.txt\" the names of " \
  " the basis elements, like \"DXDY\", one per line."

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
#include <float_image.h>
#include <float_image_read_pnm.h>
#include <float_image_map_channels.h>

#include <multifok_focus_op.h>
#include <multifok_test.h>

typedef struct mfdo_options_t 
  { char *inPrefix;                            /* Prefix for input image filenames. */
    string_vec_t image;                        /* Prefixes for input image file names. */
    multifok_focus_op_basis_type_t basisType;  /* Basis type. */
    double noise;                              /* Noise level to assume when normalizing window samples. */
    char *outPrefix;                           /* Prefix for output filenames. */
  } mfdo_options_t;
  /* Command line parameters. */

int32_t main(int32_t argc, char **argv);

void mfdo_create_basis
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
    {mfdo_write_basis_elem_names(outPrefix,NP,belName)}. */

void mfdo_process_image_set
  ( FILE *wr_ops,
    int32_t ki, 
    char *inPrefix,
    char *imageName,
    int32_t NW,
    double ws[],
    double noise,
    int32_t NB,
    double **bas,
    char *outPrefix
  );
  /* Processes the image set with index {ki} which is read from files "{inPrefix}{imageName}-{tail}" where
    {tail} is "-c.ppm" (simulated view), "-z.pgm" ({Z}-coordinates), "-d.pgm" ({Z}-deviation),
    and "-s.pgm" ("actual" sharpness {sharp}). 
    
    Writes images "{outPrefix}-{KKK}-bq.pgm" with the basis 
    coeff values SQUARED at each pixel, where {KKK} is the basis element index 
    in {0..NB-1} zero-padded to 3 digits.
    
    Also writes the local average image {avimg[ki]}, the local deviation image {dvimg[ki]}, and
    the normalized image {nimg[ki]} to files \"{outPrefix}{tail}\" where {tail} is \"-av.pgm\",
    \"-dv.pgm\", and \"-n.pgm\", respectively.
    
    Writes out to {wr_ops} one line for each pixel considered
    useful for plotting and analysis, as explained in {PROG_INFO}. */

mfdo_options_t *mfdo_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

FILE *mfdo_open_text_file(char *outPrefix, char *tag);
  /* Returns the open handle of a file called "{outPrefix}{tag}.txt" for writing. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfdo_options_t *o = mfdo_parse_options(argc, argv);
    
    /* Read the images: */
    int32_t NI = o->image.ne; /* Number of input images. */
   
    /* Get the window dimnsions: */
    int32_t NW = 3;
    demand((NW % 2 == 1) && (NW >= 3), "invalid window size");

    /* Generate the window sample weights {ws}: */
    double *ws = multifok_focus_op_sample_weights(NW);

    int32_t NB;
    double **bas = NULL;
    char **belName = NULL;
    mfdo_create_basis(NW, ws, o->basisType, &NB, &bas, &belName, o->outPrefix);

    /* File with nominal blur radius and basis element coeffs for regression: */
    FILE *wr_ops = mfdo_open_text_file(o->outPrefix, "-odata");

    for (int32_t ki = 0; ki < NI; ki++)
      { mfdo_process_image_set
          ( wr_ops, 
            ki, o->inPrefix, o->image.e[ki],
            NW, ws, o->noise, 
            NB, bas, 
            o->outPrefix
          );
      } 
      
    fclose(wr_ops);
      
    return 0;
  }
    
void mfdo_create_basis
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

void mfdo_process_image_set
  ( FILE *wr_ops,
    int32_t ki, 
    char *inPrefix,
    char *imageName,
    int32_t NW,
    double ws[],
    double noise,
    int32_t NB,
    double **bas,
    char *outPrefix
  )
  {
    int32_t NS = NW*NW; /* Number of samples in window. */

    /* Read the color image: */
    int32_t NC, NX, NY;
    float_image_t *cimg = multifok_test_read_scene_color_image(inPrefix, imageName);
    float_image_get_size(cimg, &NC, &NX, &NY);
    assert(NC == 3);
    
    /* Convert to grayscale: */
    float_image_t *gimg = float_image_new(1, NX, NY);
    float_image_map_channels_RGB_to_YUV(cimg, gimg);

    /* Read the actual sharpness image: */
    float_image_t *simg = multifok_test_read_sharpness_image(inPrefix, imageName);
    float_image_check_size(simg, 1, NX, NY);
    
    /* Read the scene {Z} average image: */
    float_image_t *zimg = multifok_test_read_zave_image(inPrefix, imageName);
    float_image_check_size(zimg, 1, NX, NY);
     
    /* Read the scene {Z} deviation image: */
    float_image_t *dimg = multifok_test_read_zdev_image(inPrefix, imageName);
    float_image_check_size(dimg, 1, NX, NY);
    
    fprintf(stderr, "allocating images for basis coefs squared...\n");
    float_image_t *bqimg[NB];
    for (int32_t kb = 0; kb < NB; kb++) 
      { bqimg[kb] = float_image_new(1, NX, NY); }

    fprintf(stderr, "allocating window average and deviation images...\n");
    float_image_t *avimg = float_image_new(1, NX, NY);
    float_image_t *dvimg = float_image_new(1, NX, NY);
      
    fprintf(stderr, "allocating the locally normalized image...\n");
    float_image_t *nimg = float_image_new(1, NX, NY);

    /* Enumerate and process all windows in the image: */
    int32_t HW = NW/2;
    float fsmp[NS];
    double dsmp[NS];
    double coeff[NB];       /* Coefficients of normalized window in be basis. */
    int32_t NP_scan = 0;    /* Number of pixels processed. */
    int32_t NP_edge = 0;    /* Number of pixels discarded because of high {zdev}. */
    int32_t NP_ops = 0;     /* Number of pixels written to {wr_ops}. */
    for (int32_t ix = HW; ix < NX-HW; ix++)
      { for (int32_t iy = HW; iy < NY-HW; iy++) 
          { /* Get the samples in the window and normalize them for brightness and contrast: */
            float_image_get_window_samples(gimg, 0,ix,iy, NW, NW, FALSE, fsmp);
            for (int32_t ks = 0; ks < NS; ks++) { dsmp[ks] = fsmp[ks]; }
            double vave, vdev; /* Weighted window average and deviation */
            multifok_focus_op_normalize_window_samples(NW, dsmp, ws, noise, &vave, &vdev); 
            
            /* Save the average, deviation, and normalized value in respective images: */
            float_image_set_sample(avimg, 0, ix, iy, (float)vave);
            float_image_set_sample(dvimg, 0, ix, iy, (float)vdev);
            float_image_set_sample(nimg, 0, ix, iy, (float)dsmp[HW*NW + HW]);
            
            /* Compute coefficients of window in basis and estimated sharpness score: */
            multifok_focus_op_compute_basis_coeffs(NW, dsmp, ws, NB, bas, coeff);
           
            /* Save coeffs squared in basis coeffs images: */
            for (int32_t kb = 0; kb < NB; kb++)
              { float_image_set_sample(bqimg[kb], 0, ix, iy, (float)(coeff[kb]*coeff[kb])); }
              
            /* Get the "true" sharpness {sharp} at this pixel. */
            double sharp = float_image_get_sample(simg, 0, ix, iy);
            assert((sharp >= 0.0) && (sharp <= 1.0));
            
            /* Get the scene {Z} height average {zave} at this pixel: */
            float zave = float_image_get_sample(zimg, 0, ix, iy);
            
            /* Get the scene {Z} height deviation {zdev} at this pixel: */
            float zdev = float_image_get_sample(dimg, 0, ix, iy);
            
            /* Decide whether the pixel is useful for analysis: */
            double zdev_max = 0.05*multifok_scene_ZMAX; /* Ignore pixels which more than this variance. */
            if (zdev > zdev_max)
              { NP_edge++; }
            else
              { double wp = sharp*sharp;
                /* double wp = 1.0; */

                /* Write the data for analysis: */
                fprintf(wr_ops, "P%d.%d.%d ", ki, ix, iy);
                fprintf(wr_ops, " %12.6f ", vave); /* Window average. */
                fprintf(wr_ops, " %12.6f ", vdev); /* Window deviation. */
                fprintf(wr_ops, " %12.6f ", wp); /* Pixel weight. */
                fprintf(wr_ops, " %14.10f ", sharp); /* "Actual" sharpness. */
                fprintf(wr_ops, " %+12.6f ", zave); /* Avg scene height rel to focus plane in pixel. */
                fprintf(wr_ops, " %12.6f ", zdev); /* Deviation of scene height in pixel. */
                for (int32_t kb = 0; kb < NB; kb++) { fprintf(wr_ops, " %16.12f", coeff[kb]); }
                fprintf(wr_ops, "\n");
                
                NP_ops++;
              }
            NP_scan++;
          }
      }
    fprintf(stderr, "scanned %d pixels\n", NP_scan); 
    fprintf(stderr, "ignored %d pixels for high {zdev}\n", NP_edge); 
    fprintf(stderr, "written %d pixels for analysis\n", NP_ops); 

    /* Write out the basis images: */
    for (int32_t kb = 0; kb < NB; kb++)
      { char *tag = NULL;
        asprintf(&tag, "-%s-%03d", imageName, kb);
        multifok_test_write_basis_coeff_squared_image(bqimg[kb], outPrefix, tag);
        free(tag);
      }

     /* Write the average, deviation, and normalized images: */
     char *tag = NULL;
     asprintf(&tag, "-%s", imageName);
     multifok_test_write_window_average_image(avimg, outPrefix, tag);
     multifok_test_write_window_deviation_image(dvimg, outPrefix, tag);
     multifok_test_write_normalized_image(nimg, outPrefix, tag);
     free(tag);
  }

FILE *mfdo_open_text_file(char *outPrefix, char *tag)
  { char *fname = NULL;
    asprintf(&fname, "%s%s.txt", outPrefix, tag);
    FILE *wr = open_write(fname, TRUE);
    return wr;
    free(fname);
  }

mfdo_options_t *mfdo_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    mfdo_options_t *o = notnull(malloc(sizeof(mfdo_options_t)), "no mem");

    argparser_get_keyword(pp, "-inPrefix");
    o->inPrefix = argparser_get_next(pp);
    
    o->image = string_vec_new(20);
    int32_t NI = 0; /* Number of input image sets. */
    while (argparser_keyword_present(pp, "-image"))
      { string_vec_expand(&(o->image), NI);
        o->image.e[NI] = argparser_get_next(pp);
        NI++;
      }
    string_vec_trim(&(o->image), NI);
    if (NI == 0) { argparser_error(pp, "must specify at least one \"-image\""); }

    argparser_get_keyword(pp, "-basisName");
    char *bName = argparser_get_next_non_keyword(pp);
    o->basisType = multifok_focus_op_basis_type_from_text(bName, FALSE);
    if (bName < 0) { argparser_error(pp, "invalid basis type"); }

    argparser_get_keyword(pp, "-noise");
    o->noise = argparser_get_next_double(pp, 0.0, 100.0);  

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
