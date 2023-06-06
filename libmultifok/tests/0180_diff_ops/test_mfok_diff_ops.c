#define PROG_NAME "test_mfok_diff_ops"
#define PROG_DESC "Amalyzes dependency of differential operators on sharpness.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-04-28 17:30:53 by stolfi */ 
/* Created on 2023-01-24 by J. Stolfi, UNICAMP */

#define test_mfok_diff_ops_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "Duh?"

#define PROG_INFO \
  "SYNOPSIS\n" \
  "  Tries to indentify an optimum local focus sharpness formula from" \
  " images where that information is known.\n" \
  "\n" \
  "  Reads {NI} images with focus blur.  For each image, enumerates" \
  " the {NW} by {NW} windows centered at each pixel.  For each window" \
  " position, maps the {NS = NW*NW} image samples in the window to" \
  " coefficients {coeff[0..NB-1]} of a  specified local operator" \
  " basis with {NB} elements.  Writes {NI*NB} images showing those" \
  " coefficients.  Then computes {NT} specified quadratic" \
  " terms {term[0..NT-1]} from those coefficients.  Writes" \
  " this data as a text file for plotting and regression.\n" \
  "\n" \
  "  Before computing the basis coefficients, he samples in each" \
  " window are normalized to have zero average and gradient and unit deviation.  Thus" \
  " the basis coefficients and the quadratic terms are independent of local" \
  " image brightness, gradient, and contrast.\n" \
  "INPUTS\n" \
  "  INPUT IMAGES\n" \
  "    Reads a collection of image quadruples {csimg[ki]}, {azimg[ki]}," \
  " {dzimg[ki]}, {shimg[ki]} for {ki} in {0..NI-1}, where\n" \
  "\n" \
  "      {csimg[ki]} is sythetic image of a 3D scene  with simulated focus blurring.\n" \
  "      {azimg[ki]} specifies the average {Z}-coord of the scene within each pixel.\n" \
  "      {dzimg[ki]} specifies the deviation of that {Z} coordinate within each pixel.\n" \
  "      {shimg[ki]} specifies the actual sharpness of {csimg[k]} at each pixel, from the ray tracing.\n" \
  "\n" \
  "    The {csimg} images should be in color, but are converted internally to a" \
  " grayscale image {grimg[ki]}. The other images are greyscale with linear" \
  " encoding (gamma = 1).  They may have different sizes and focus plane" \
  " positions {zFoc}, but must all have the same depth of focus {zDep}.\n" \
  "\n" \
  "  TERMS TABLE\n" \
  "    The formula terms {term[0..NT-1]} are" \
  " computed from the basis coeffs as described" \
  " under {multifok_focus_op_terms_from_basis_coeffs}, using" \
  " the terms descriptions read from the specified {termSetFile}. See" \
  " {multifok_term_indices_from_names} and {multifok_term_read_term_table} for the file format.\n" \
  "\n" \
  "OUTPUTS\n" \
  "\n" \
  "  BASIS IMAGE FILES\n" \
  "    For each input image {csimg[ki]} and each basis element {kb} in {0..NB-1}, writes" \
  " to \"{outPrefix}-{NNN}-bq.pgm\" a grayscale" \
  " image showing the SQUARED coefficient of that basis element computed on normalized window" \
  " samples; where {NNN} is the element index {kb} zero-padded to 3 digits." \
  "\n" \
  "  TERM IMAGE FILES\n" \
  "    For each input image {csimg[ki]} and each quadratic term {kt} in {0..NT-1}, writes" \
  " to \"{outPrefix}-{NNN}-tm.pgm\" a grayscale" \
  " image showing the value of that quadratic term computed from the basis" \
  " coefficients above; where {NNN} is the term index {kt} zero-padded to 3 digits." \
  "\n" \
  "  LOCAL AVERAGE AND DEVIATION IMAGE FILES\n" \
  "    For each input image {csimg[ki]} the program computes images {avimg[ki]} and {dzimg[ki]} with" \
  " the window sample average and deviation around each pixel of the grayscale image {grimg[ki]}, and writes" \
  " them to \"{outPrefix}{tail}\" where {tail} is \"-av.pgm\" and \"-dv.pgm\", respectively.\n" \
  "\n" \
  "  LOCALLY NORMALIZED IMAGE FILES\n" \
  "    For each input image {csimg[ki]} writes" \
  " to \"{outPrefix}-n.txt\" a grayscale image showing the central window pixel after window" \
  " normalization, mapped from {[-1_+1]} to {[0_1]}.\n" \
  "\n" \
  "  REGRESSION DATA FILE\n" \
  "    Writes to \"{outPrefix}-odata.txt\" a file with columns\n" \
  "\n" \
  "    P{ki}.{ix}.{iy} {vave} {vgrd} {vdev} {sharp} {zrav} {zdev} {coeff[0]} .. {coeff[NB-1]} {term[0]} .. {term[NT-1]} \n" \
  "\n" \
  " where \n" \
  "\n" \
  "    {ki} is the image index.\n" \
  "    {ix} and {iy} are the column and row of the pixel.\n" \
  "    {vave} is the average of the samples in the window centered at that pixel.\n" \
  "    {vgrd} is the gradient modulus within that window.\n" \
  "    {vdev} is the rms of window samples values after removing average and gradient.\n" \
  "    {sharp} is the \"actual\" sharpness at that pixel (as read from the {shimg[ki]} image).\n" \
  "    {zrav} is average of the {Z}-coord of the scene within that pixel, rel to focused plane {Z}.\n" \
  "    {zdev} is the deviation of the {Z}-coord of the scene within that pixel.\n" \
  "    {coeff[0..NB-1]} are the basis coefficients computed from the normalized window samples.\n" \
  "    {term[0..NT-1]} are the quadratic terms computed from those coefficients.\n" \
  "\n" \
  "  The values of {vave} and {vdev} are computed from the window of {grimg[hi]} centered" \
  " at the pixel, before they are normalized, taking the window sample weights" \
  " into account.  The basis coefficients and quadratic terms are computed after the" \
  " window samples have been normalized for brightness, gradient, and contrast.\n" \
  "\n" \
  "  PLOT MASK IMAGE FILE\n" \
  "    Writes to \"{outPrefix}-pixmask.pgm\" a binary image file that is 1 for pixels" \
  " that were written out to the \"-odata\" file, 0 otherwise.\n" \
  "\n" \
  "  BASIS NAMES, PRODUCTS, AND TERM FILE\n" \
  "    Also writes, for documentation, the following files:\n" \
  "\n" \
  "      \"{outPrefix}-bnames.txt\"  names of the {NB} basis elements, like \"DXDY\";\n" \
  "      \"{outPrefix}-prix.txt\"    table mapping products of pairs of coefficients to terms;\n" \
  "      \"{outPrefix}-tnames.txt\"  list of quadratic term formulas, like \"DX*DX+DY*DY\"." \
  "\n" \
  "    See {multifok_term_read_index_table} for the format and semantics of the \"-prix\" file."

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
#include <float_image.h>
#include <float_image_read_pnm.h>
#include <float_image_map_channels.h>

#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_term.h>
#include <multifok_test.h>
#include <multifok_scene.h>

typedef struct mfdo_image_spec_t
  { char *sceneType;  /* Scene type shown in image. */
    int32_t NX;       /* Width of image.  */
    int32_t NY;       /* Height of image.  */
    char *pattern;    /* Name of pattern used to texture image. */
    double zDep;      /* Depth of focus of image. */
    double zFoc;      /* Focus plane position of image. */
  } mfdo_image_spec_t;
  /* Data that determines the name of each input image file. */
  
vec_typedef(mfdo_image_spec_vec_t,mfdo_image_spec_vec,mfdo_image_spec_t);

#define ZMAX multifok_scene_ZMAX

typedef struct mfdo_options_t 
  { /* Input images: */
    char *inDir;                               /* Directory where input image files live. */
    mfdo_image_spec_vec_t image;               /* Specifications of the input image files. */
    multifok_basis_type_t basisType;           /* Local operator basis type. */
    char *termSetFile;                         /* Name of file with formulas of quadratic terms. */
    double noise;                              /* Noise level to assume when normalizing window samples. */
    char *outPrefix;                           /* Prefix for output filenames. */
  } mfdo_options_t;
  /* Command line parameters. */

int32_t main(int32_t argc, char **argv);

void mfdo_create_basis
  ( int32_t NW,
    double ws[],
    multifok_basis_type_t bType,
    int32_t *NB_P,
    double ***bas_P,
    char ***belName_P,
    char *outPrefix
  );
  /* Creates an operator basis of the specified {basisType} for an {NW} by {NW} window with
    with window sample weights {ws[0..NS-1]}, where {NS=NW*NW}.  See {multifok_basis_make}. 
    The basis will be orthonormal.  The number {NB} of basis elements is returned in {*NB_P},
    the basis elements {bas[0..NB-1][0..NS-1]}  are returned in {*bas_P}, and their names
    {belName[0..NB-1]} are returned in {*belName_P}.
    
    Also prints the basis and writes its names out to file with
    {mfdo_write_basis_elem_names(outPrefix,NP,belName)}. */

void mfdo_read_term_names(char *fname, int32_t *NT_P, char ***termName_P);
  /* Reads term names (formulas) from file {fname}.
    See {multifok_term_read_names} for the file format. */
    
void mfdo_process_image_set
  ( FILE *wr_dat,
    int32_t ki, 
    char *inDir,
    mfdo_image_spec_t *imo,
    int32_t NW,
    double ws[],
    double noise,
    int32_t NB,
    double **bas,
    int32_t NP,
    multifok_term_prod_t prix[],
    int32_t NT,
    char *outPrefix
  );
  /* Processes the image set with index {ki} which is read from files "{inPrefix}{imageName}-{tail}" where
    {tail} is "-cs.ppm" (simulated view), "-az.pgm" ({Z}-coordinates), "-dz.pgm" ({Z}-deviation),
    and "-sh.pgm" ("actual" sharpness {sharp}). 
    
    For each {kb} in {0..NB-1}, writes an image "{outPrefix}-{KKK}-bq.pgm" with the coefficient SQUARED
    of basis element {kb} at each pixel, where {KKK} is {kb} zero-padded to 3 digits.
    
    Also, for each {kt} in {0..NT-1}, writes an image "{outPrefix}-{KKK}-tm.pgm" with the value of quadratic 
    term {kt} at each pixel, where {KKK} is {kt} zero-padded to 3 digits.
    
    Also writes the local average image {avimg[ki]}, the local deviation image {dzimg[ki]}, and
    the normalized image {nrimg[ki]} to files "{outPrefix}{tail}" where {tail} is "-av.pgm",
    "-dv.pgm", and "-nr.pgm", respectively.
    
    Writes out to {wr_dat} one line for each pixel considered
    useful for plotting and analysis, as explained in {PROG_INFO}.
    
    Also writes a mask image file "{outPrefix}-mk.pgm" that shows which pixels 
    had their data written out. */

void mfdo_read_term_index_table
  ( char *fname, 
    int32_t NB,
    char *belName[],
    int32_t *NP_P,
    multifok_term_prod_t **prix_P,
    int32_t *NT_P,
    char ***termName_P
  );
  /* Calls {multifok_term_read_index_table} on the file "{fname}". */

bool_t mfdo_pixel_is_useful
  ( double vave, 
    double vgrd, 
    double vdev, 
    double zave, 
    double zdev, 
    double zFoc, 
    double zDep,
    int32_t *nd_vdev_P,
    int32_t *nd_zdev_P,
    int32_t *nd_zave_P,
    int32_t *nd_zrav_P
  );
  /* Decides if a pixel is useful for analysis or regression, based on
    the average {vave}, gradient modulus {vgrd}, anddeviation {vdev} of
    window samples, the average {zave} and deviation {zdev} of absolute
    scene {Z} within the pixel, the {Z} position of the focused plane
    {zFoc}, and the focus depth {zDep}.
    
    Increments the counters on these conditions:
      
      {*nr_vdev_P} Sample deviation {vdev} in window too high.
      {*nr_zdev_P} Scene {Z} deviation in sampled area too high.
      {*nr_zave_P} Scene {Z} average in sampled area too low (fat floor).
      {*nr_zrav_P} Scene {Z} average too far from focused plane position {zFoc}. 
  */

mfdo_options_t *mfdo_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

/* IMPLEMENTATIONS */

vec_typeimpl(mfdo_image_spec_vec_t,mfdo_image_spec_vec,mfdo_image_spec_t);

int32_t main (int32_t argc, char **argv)
  {
    mfdo_options_t *o = mfdo_parse_options(argc, argv);
    
    /* Read the images: */
    int32_t NI = o->image.ne; /* Number of input images. */
   
    /* Get the window size: */
    int32_t NW = 3;
    demand((NW % 2 == 1) && (NW >= 3), "invalid window size");

    /* Generate the window sample weights {ws}: */
    double *ws = multifok_window_sample_weights(NW);

    /* Create the basis elements:*/
    int32_t NB;
    double **bas = NULL;
    char **belName = NULL;
    mfdo_create_basis(NW, ws, o->basisType, &NB, &bas, &belName, o->outPrefix);

    /* Get term descriptions: */
    int32_t NT;   /* Number of quadratic terms to analyze. */
    char **termName; /* Term names. */
    int32_t NP;   /* Number of coeff products in the quadratic terms. */
    multifok_term_prod_t *prix;    /* Mapping of basis element index pairs to terms. */
    multifok_test_read_term_weights_names_get_indices
      ( o->termSetFile,
        NB, belName,
        &NT, NULL, &termName,
        &NP, &prix,
        TRUE
      ); 
    multifok_test_write_term_names(o->outPrefix, NT, termName);
    multifok_test_write_term_index_table(o->outPrefix, NP, prix, TRUE);

    /* File with basis coeffs and quadratic terms for regression: */
    FILE *wr_dat = multifok_test_open_text_file(o->outPrefix, "-odata");

    for (int32_t ki = 0; ki < NI; ki++)
      { mfdo_process_image_set
          ( wr_dat, 
            ki, o->inDir, &(o->image.e[ki]),
            NW, ws, o->noise, 
            NB, bas,
            NP, prix, 
            NT,
            o->outPrefix
          );
      } 
      
    fclose(wr_dat);
      
    return 0;
  }
    
void mfdo_create_basis
  ( int32_t NW,
    double ws[],
    multifok_basis_type_t bType,
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
    multifok_basis_make(bType, NW, ws, ortho, &NB, &bas, &belName);
    fprintf(stderr, "obtained NB = %d independent basis elements\n", NB);
    multifok_basis_print(stderr, NW, NB, bas, belName);
    multifok_basis_ortho_check(stderr, NW, NB, bas);
      
    /* File with basis element names: */
    multifok_test_write_basis_elem_names(outPrefix, NB, belName);
      
    (*NB_P) = NB;
    (*bas_P) = bas;
    (*belName_P) = belName;
  }

void mfdo_process_image_set
  ( FILE *wr_dat,
    int32_t ki, 
    char *inDir,
    mfdo_image_spec_t *imo,
    int32_t NW,
    double ws[],
    double noise,
    int32_t NB,
    double **bas,
    int32_t NP,
    multifok_term_prod_t prix[],
    int32_t NT,
    char *outPrefix
  )
  {
    int32_t NS = NW*NW; /* Number of samples in window. */
    
    int32_t NX = imo->NX;
    int32_t NY = imo->NY;
    
    /* Assemble the image name: */
    char *inPrefix = NULL;
    asprintf(&inPrefix, "%s/st%s-%04dx%04d-%s/frame", inDir, imo->sceneType, NX, NY, imo->pattern);
    char *frameTag = NULL;
    asprintf(&frameTag, "-fd%05.2f-zf%05.2f", imo->zDep, imo->zFoc);

    /* Read the color image: */
    int32_t NC_csimg = 3;
    float_image_t *csimg = multifok_test_read_scene_color_image(inPrefix, frameTag);
    float_image_check_size(csimg, NC_csimg, NX, NY);
    
    /* Convert to grayscale: */
    float_image_t *grimg = float_image_new(1, NX, NY);
    float_image_map_channels_RGB_to_YUV(csimg, grimg); /* Discard {U,V} channels. */

    /* Read the actual sharpness image: */
    float_image_t *shimg = multifok_test_read_sharpness_image(inPrefix, frameTag);
    float_image_check_size(shimg, 1, NX, NY);
    
    /* Read the scene {Z} average image: */
    float_image_t *azimg = multifok_test_read_zave_image(inPrefix, frameTag);
    float_image_check_size(azimg, 1, NX, NY);
     
    /* Read the scene {Z} deviation image: */
    float_image_t *dzimg = multifok_test_read_zdev_image(inPrefix, frameTag);
    float_image_check_size(dzimg, 1, NX, NY);
    
    fprintf(stderr, "allocating images for basis coefs squared...\n");
    float_image_t *bqimg[NB];
    for (int32_t kb = 0; kb < NB; kb++) 
      { bqimg[kb] = float_image_new(1, NX, NY); }

    fprintf(stderr, "allocating images for quadratic term values...\n");
    float_image_t *tmimg[NT];
    for (int32_t kt = 0; kt < NT; kt++) 
      { tmimg[kt] = float_image_new(1, NX, NY); }

    fprintf(stderr, "allocating window average, gradient, and deviation images...\n");
    float_image_t *avimg = float_image_new(1, NX, NY);
    float_image_t *gvimg = float_image_new(1, NX, NY);
    float_image_t *dvimg = float_image_new(1, NX, NY);
      
    fprintf(stderr, "allocating the locally normalized image...\n");
    float_image_t *nrimg = float_image_new(1, NX, NY);

    fprintf(stderr, "allocating the pixel selection mask...\n");
    float_image_t *mkimg = float_image_new(1, NX, NY);

    /* Enumerate and process all windows in the image: */
    int32_t HW = NW/2;
    float fsmp[NS];
    double dsmp[NS];
    double coeff[NB];       /* Coefficients of normalized window in be basis. */
    double term[NT];        /* Quadratic terms. */
    int32_t NL_scan = 0;    /* Number of pixels processed. */
    int32_t nd_vdev_P = 0;  /* Number of pixels rejected for high {vdev}. */ 
    int32_t nd_zdev_P = 0;  /* Number of pixels rejected for high {zdev}. */ 
    int32_t nd_zave_P = 0;  /* Number of pixels rejected for low {zave}. */ 
    int32_t nd_zrav_P = 0;  /* Number of pixels rejected for {zave} too far from {zFoc}. */ 
    int32_t NL_dat = 0;     /* Number of pixels written to {wr_dat}. */
    for (int32_t ix = HW; ix < NX-HW; ix++)
      { for (int32_t iy = HW; iy < NY-HW; iy++) 
          { /* Get the samples in the window and normalize them for brightness and contrast: */
            float_image_get_window_samples(grimg, 0,ix,iy, NW, NW, FALSE, fsmp);
            for (int32_t ks = 0; ks < NS; ks++) { dsmp[ks] = fsmp[ks]; }
            double vave, vgrd, vdev; /* Weighted window average and deviation */
            multifok_window_normalize_samples(NW, dsmp, ws, noise, &vave, &vgrd, &vdev); 
            
            /* Save the average, deviation, and normalized value in respective images: */
            float_image_set_sample(avimg, 0, ix, iy, (float)vave);
            float_image_set_sample(gvimg, 0, ix, iy, (float)vgrd);
            float_image_set_sample(dvimg, 0, ix, iy, (float)vdev);
            float_image_set_sample(nrimg, 0, ix, iy, (float)dsmp[HW*NW + HW]);
            
            /* Compute coefficients of window in basis: */
            multifok_basis_compute_coeffs(NW, dsmp, NB, bas, coeff);
           
            /* Save coeffs squared in basis coeffs images: */
            for (int32_t kb = 0; kb < NB; kb++)
              { float_image_set_sample(bqimg[kb], 0, ix, iy, (float)(coeff[kb]*coeff[kb])); }
              
            /* Compute the quadratic terms: */
            multifok_term_values_from_basis_coeffs(NB, coeff, NP, prix, NT, term);
           
            /* Save terms in term images: */
            for (int32_t kt = 0; kt < NT; kt++)
              { float_image_set_sample(tmimg[kt], 0, ix, iy, (float)(term[kt])); }
              
            /* Get the "true" sharpness {sharp} at this pixel. */
            double sharp = float_image_get_sample(shimg, 0, ix, iy);
            assert((sharp >= 0.0) && (sharp <= 1.0));
            
            /* Get the absolute scene {Z} height average {zave} at this pixel: */
            float zave = float_image_get_sample(azimg, 0, ix, iy);
            double zrav = zave - imo->zFoc; /* Average scene {Z} rel to focused plane. */
            
            /* Get the scene {Z} height deviation {zdev} at this pixel: */
            float zdev = float_image_get_sample(dzimg, 0, ix, iy);
            
            /* Decide whether the pixel is useful for analysis: */
            bool_t pix_ok = mfdo_pixel_is_useful
              ( vave, vgrd, vdev, zave, zdev, 
                imo->zFoc, imo->zDep,
                &nd_vdev_P, &nd_zdev_P, &nd_zave_P, &nd_zrav_P 
              );
            float mkval;  /* Value to write in pixel mask image. */
            if (pix_ok)
              { mkval = 0.0;
              }
            else
              { /* Write the data for analysis: */
                fprintf(wr_dat, "P%d.%d.%d ", ki, ix, iy);
                fprintf(wr_dat, " %12.6f ", vave);   /* Window average. */
                fprintf(wr_dat, " %12.6f ", vgrd);   /* Window gradient modulus. */
                fprintf(wr_dat, " %12.6f ", vdev);   /* Window deviation after removing average and gradient. */
                fprintf(wr_dat, " %14.10f ", sharp); /* "Actual" sharpness. */
                fprintf(wr_dat, " %+12.6f ", zrav);  /* Avg scene height rel to focus plane in pixel. */
                fprintf(wr_dat, " %12.6f  ", zdev);  /* Deviation of scene height in pixel. */
                fprintf(wr_dat, "   ");
                for (int32_t kb = 0; kb < NB; kb++) { fprintf(wr_dat, " %16.12f", coeff[kb]); }
                fprintf(wr_dat, "   ");
                for (int32_t kt = 0; kt < NT; kt++) { fprintf(wr_dat, " %16.12f", term[kt]); }
                fprintf(wr_dat, "\n");
                
                NL_dat++;
                mkval = 1.0;
              }
            float_image_set_sample(mkimg, 0, ix, iy, mkval);
            NL_scan++;
          }
      }
    fprintf(stderr, "scanned %d pixels\n", NL_scan); 
    fprintf(stderr, "rejected %d pixels for high {vdev}\n", nd_vdev_P);  
    fprintf(stderr, "rejected %d pixels for high {zdev}\n", nd_zdev_P);  
    fprintf(stderr, "rejected %d pixels for low {zave}\n", nd_zave_P);  
    fprintf(stderr, "rejected %d pixels for {zave} too far from {zFoc}\n", nd_zrav_P);  
    fprintf(stderr, "wrote %d pixels for analysis\n", NL_dat); 

    demand(NL_dat > 0, "ABORTED");

    /* Prefix for outopt image files: */
    char *outImagePrefix = NULL;
    asprintf(&outImagePrefix, "%s-img-st%s-%04dx%04d-%s", outPrefix, imo->sceneType, NX, NY, imo->pattern);

    /* Write out the squared basis coefficient images: */
    for (int32_t kb = 0; kb < NB; kb++)
      { char *frameBelTag = NULL;
        asprintf(&frameBelTag, "%s-kb%03d", frameTag, kb);
        multifok_test_write_basis_coeff_squared_image(bqimg[kb], outImagePrefix, frameBelTag);
        free(frameBelTag);
      }

    /* Write out the quadratic term images: */
    for (int32_t kt = 0; kt < NT; kt++)
      { char *frameTermTag = NULL;
        asprintf(&frameTermTag, "%s-kt%03d", frameTag, kt);
        multifok_test_write_quadratic_term_image(tmimg[kt], outImagePrefix, frameTermTag);
        free(frameTermTag);
      }

    /* Write the average, deviation, and normalized images: */
    multifok_test_write_window_average_image(avimg, outImagePrefix, frameTag);
    multifok_test_write_window_gradient_image(gvimg, outImagePrefix, frameTag);
    multifok_test_write_window_deviation_image(dvimg, outImagePrefix, frameTag);
    multifok_test_write_normalized_image(nrimg, outImagePrefix, frameTag);

    /* Write the pixel selection mask image: */
    multifok_test_write_pixel_mask_image(mkimg, outImagePrefix, frameTag);

    free(inPrefix);
    free(frameTag);
    free(outImagePrefix);
  }
            
bool_t mfdo_pixel_is_useful
  ( double vave, 
    double vgrd, 
    double vdev, 
    double zave, 
    double zdev, 
    double zFoc, 
    double zDep,
    int32_t *nd_vdev_P,
    int32_t *nd_zdev_P,
    int32_t *nd_zave_P,
    int32_t *nd_zrav_P
  )
  { double vdev_min = 0.02; /* Reject pixels with {vdev} below this value. */
    if (vdev < vdev_min) { (*nd_vdev_P)++; return FALSE; } 
    double zdev_max = 0.75; /* Ignore pixels which more than this variance. */
    if (zdev > zdev_max) { (*nd_zdev_P)++; return FALSE; }
    double zave_min = 1.1;       /* Ignore pixels with actual scene {Z} below this level. */
    if (zave < zave_min) { (*nd_zave_P)++; return FALSE; }
    double zrav_max = 10.0*zDep; /* Ignore pixels with abs relative {Z} above this value. */
    double zrav = zave - zFoc;
    if (fabs(zrav) > zrav_max) { (*nd_zrav_P)++; return FALSE; }
    return TRUE;
  }

mfdo_options_t *mfdo_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    mfdo_options_t *o = notnull(malloc(sizeof(mfdo_options_t)), "no mem");

    argparser_get_keyword(pp, "-inDir");
    o->inDir = argparser_get_next(pp);
    
    o->image = mfdo_image_spec_vec_new(20);
    int32_t NI = 0; /* Number of input images. */
    while (argparser_keyword_present(pp, "-image"))
      { mfdo_image_spec_t imo;
        imo.sceneType = argparser_get_next_non_keyword(pp);
        imo.NX = (int32_t)argparser_get_next_int(pp,50,4096);
        imo.NY = (int32_t)argparser_get_next_int(pp,50,4096);
        imo.pattern = argparser_get_next_non_keyword(pp);
        imo.zDep = argparser_get_next_double(pp, 0.01,99.99);
        imo.zFoc = argparser_get_next_double(pp, 0.0,2*ZMAX);
        
        mfdo_image_spec_vec_expand(&(o->image), NI);
        o->image.e[NI] = imo;
        NI++;
      }
    mfdo_image_spec_vec_trim(&(o->image), NI);
    if (NI == 0) { argparser_error(pp, "must specify at least one \"-image\""); }

    argparser_get_keyword(pp, "-basisName");
    char *bName = argparser_get_next_non_keyword(pp);
    o->basisType = multifok_basis_type_from_text(bName, FALSE);
    if (bName < 0) { argparser_error(pp, "invalid basis type"); }
      
    argparser_get_keyword(pp, "-termSetFile");
    o->termSetFile = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-noise");
    o->noise = argparser_get_next_double(pp, 0.0, 100.0);  

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
