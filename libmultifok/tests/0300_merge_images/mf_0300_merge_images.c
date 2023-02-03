#define PROG_NAME "mf_0300_merge_images"
#define PROG_DESC "Merges several registered images with focus blur at different heights"
#define PROG_VERS "1.0"

/* Last edited on 2023-01-29 13:58:31 by stolfi */ 
/* Created on 2023-01-19 by J. Stolfi, UNICAMP */

#define mf_0300_merge_images_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "Duh?"

#define PROG_INFO \
  "SYNOPSIS" \
  "\n" \
  "  Merges a stack of registerd images taken with same lighting and" \
  " view geometry but with short depth of focus and different focus plane position\n" \
  "\n" \
  "INPUTS\n" \
  "\n" \
  "  BLURRED IMAGES\n" \
  "    The program reads a stack of {NI} images {csimg[ki]}, {shimg[ki]}, where\n" \
  "\n" \
  "      {csimg[ki]} is sythetic image a scene with simulated focus blurring.\n" \
  "      {shimg[ki]} specifies the actual sharpness of {csimg[k]} at each pixel.\n" \
  "\n" \
  "    These images are read from files \"{inPrefix}-z{ZZZ.ZZZZ}{tail}\" where {ZZZ.ZZZZ} is" \
  " the coordinate {zFoc[ki]} of the focus plane, formated as \"%8.4f\", and" \
  " {tail} is \"-cs.ppm\" or \"-sh.pgm\", respectively.  The number {NI} of" \
  " images and the values of {zFoc[0..NI-1]} are" \
  " specified by the \"-Zrange\" and \"-zStep\" options.\n" \
  "\n" \
  "  REFERENCE IMAGES\n" \
  "    Also, for comparison purposes, reads \"unblurred\" images {azimg}, {dzimg}, where\n" \
  "\n" \
  "      {timg} an unblurred color view of the scene.\n" \
  "      {azimg} specifies the average {Z}-coord of the scene within each pixel.\n" \
  "      {dzimg} specifies the deviation of that {Z} coordinate within each pixel.\n" \
  "\n" \
  "    These images are read from files \"{inPrefix}-sharp{tail}\" where" \
  " {tail} is \"-cs.ppm\", \"-az.pgm\" or \"-g.pgm\", respectively.\n" \
  "\n" \
  "  All images must be of same scene and must have the same size, lighting, and viewing" \
  " parameters. The {csimg} images are in color, the other images are greyscale with linear" \
  " encoding (gamma = 1).  The images {csimg[ki]} and {shimg[ki]} must have the" \
  " same depth of focus {zDep}, with the focus plane at" \
  " the {Z} coordinate {zFoc[ki]}.  The samples of the image {azimg} are assumed" \
  " to be {Z} coordinates relative to {sharpImage_zFoc}. " \
  " scaled as described in {mutifok_test_write_zave_image}.\n" \
  "\n" \
  "  TERM NAMES AND WEIGHTS\n" \
  "    The term weights are read from a file \"term-weights/ptBEST-bt{BBBB}-trm{T}-qsh0-ns{NNNN}\" where" \
  " {BBBB} and {T} are chosen by the program, and {NNNN} is the value of {noise} formatted with \"%05.2f\".\n" \
  "\n" \
  "OUTPUTS\n" \
  "\n" \
  "  RECONSTRUCTED COLOR IMAGE\n" \
  "    Writes to \"{outPrefix}-cr.ppm\" an image {crimg} reconstructed from the images {csimg[0..NI-1]}, using" \
  " at each pixel a set of weights determined from the sharpness scores obtained at that pixel on those" \
  " images.\n" \
  "\n" \
  "  ESTIMATED Z IMAGE\n" \
  "    Also writes to \"{outPrefix}-cz.pgm\" an image {czimg} with the {zEst} for each pixel, the" \
  " estimated {Z} depth of the scene, based" \
  " on the sharpness scores and the focus plane positions {zFoc[ki]}.\n" \
  "\n" \
  "  ESTIMATED Z ERROR IMAGE\n" \
  "    Also writes to \"{outPrefix}-ez.pgm\" an image {ezimg}  with the difference" \
  " between the estimated {Z} depth of the scene {zEst} and the \"actual\" depth {zAve} as read from the {azimg} image.\n" \
  "\n" \
  "  COLOR ERROR IMAGE\n" \
  "    Also writes to \"{outPrefix}-er.ppm\" a color image {erimg}  with the difference" \
  " between the reconstructed color estimated at each pixel, as in the {crimg} image,  and the \"actual\" color as in the {timg}.\n" \
  "\n" \
  "  TERMS AND WEIGHTS\n" \
  "    Also writes to \"{outPrefix}-terms.txt\" the term names and weights used, for documentation."

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
#include <lsq.h>

#include <multifok_focus_op.h>
#include <multifok_test.h>

typedef struct mfmi_options_t 
  { char *inPrefix;      /* Prefix for input filenames. */
    double zRange_lo;    /* {Z} value of lowest frame. */
    double zRange_hi;    /* {Z} value of highest frame. */
    double zStep;        /* {Z} increment between frames. */
    double focDepth;     /* Depth of focus of blurred images. */
    bool_t actualSharp;  /* Use actual sharpness instead of estimated one. */
    double noise;        /* Assumed noise level in image samples. */
    char *outPrefix;     /* Prefix for output filenames. */
  } mfmi_options_t;
  /* Command line parameters. */

#define ZMAX multifok_scene_ZMAZ

int32_t main(int32_t argc, char **argv);

mfmi_options_t *mfmi_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

void mfmi_choose_basis_and_terms
  ( int32_t NW,
    double ws[],
    int32_t *NB_P,
    double ***bas_P,
    char ***belName_P,
    double noise, 
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P,
    int32_t *NT_P, 
    double **wt_P, 
    char ***termName_P,
    char *outPrefix
  );
  /* Chooses and creates an operator basis for an {NW} by {NW} window
    with with window sample weights {ws[0..NS-1]}, where {NS=NW*NW}. See
    {multifok_basis_make}. The basis will be orrthonormal. The
    number {NB} of basis elements is returned in {*NB_P}, the basis
    elements {bas[0..NB-1][0..NS-1]} are returned in {*bas_P}, and their
    names {belName[0..NB-1]} are returned in {*belName_P}.
    
    Also prints the basis and writes its names out to file
    "{outPrefix}-belnames.txt" with
    {multifok_test_write_basis_elem_names(outPrefix,NP,belName)}.
    
    Then it defines the quadratic formula to be used to
    compute the sharpness score. That includes the number of independent quadratic terms {NT}, their formulas {termName[0..NT-1]}, and
    and their respective weights {wt[0..NT-1]}. It also defines the number of distinct canonical coeff products {NP}
    and the table {prix[0..NP-1]} that tells which products are part of which terms.  These data are returned in
    {*NP_P,*prix_P,*NT_P,*wt_P,*termName_P}.
    
    The term names and weights are read from a file
    "term-weights/ptBEST-bt{BASNAME}-trm{T}-qsh0-ns{NOISE}.txt" where
    {BASNAME} is the basis name ("LAPL", "DIFF", etc.), {T} is the term
    subset to use ("1" for all diff ops squared plus unity, "2" for the
    squared order 1 and 2 diff ops only, etc.) and {NOISE} is the
    {noise} bias formatted with "%05.2f".
    
    Also writes the terms, weights, and basis index pairs to {stderr} and writes the names and weights to a file "{outPrefix}-termName.txt" with
    {multifok_test_write_term_names(outPrefix,NT,termName)}. */

void mfmi_choose_terms
  ( int32_t NB, 
    char *belName[], 
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P,
    int32_t *NT_P, 
    double **wt_P, 
    char ***termName_P,
    char *outPrefix
  )
  ( char *outPrefix
  );
  /* Given the basis element names {belName[0..NB-1]}, ?? */

void mfmi_process_image_stack
  ( FILE *wr_cmp,
    int32_t NI, 
    float_image_t *csimg[], 
    float_image_t *grimg[], 
    float_image_t *shimg[],
    double zFoc[],
    float_image_t *timg, 
    float_image_t *azimg, 
    float_image_t *dzimg, 
    bool_t actualSharp,
    int32_t NW,
    double ws[],
    double noise,
    int32_t NB,
    double **bas,
    int32_t NT,
    double wt[],
    int32_t NP,
    multifok_term_prod_t prix[],
    char *outPrefix
  );
  /* Processes the image stack {csimg[0..NI-1]} 
    
    Writes images "{outPrefix}-{KKK}-bc.pgm" with the basis 
    coeff values at each pixel, where {KKK} is the basis element index 
    in {0..NB-1} zero-padded to 3 digits.
    
    Also writes images {scimg} = "{outPrefix}-sc.pgm" with the sharpness {score} computed
    from the basis coeffs, and {esimg} = "{outPrefix}-es.pgm" with the error {score-sharp}.
    See {multifok_score_from_basis_coeffs} for how the {score} is computed.
    
    Writes out to {wr_reg} and {wr_cmp} one line for each pixel considered
    useful for the regression, as explained in {PROG_INFO}. */

void mfmi_compute_pixel_scores
  ( int32_t NI, 
    float_image_t *grimg[],
    int32_t ix, int32_t iy,
    int32_t NW,
    double ws[],
    double noise,
    int32_t NB,
    double **bas,
    int32_t NT,
    double wt[],
    int32_t NP,
    multifok_term_prod_t prix[],
    double score[]
  );
  /* Estimates the sharpness scores {score[0..NI-1]} of grayscale images {grimg[0..NI-1]} at the pixel {ix,iy},
    as per {multifok_score_from_basis_coeffs}. */

void mfmi_get_pixel_actual_sharpness
  ( int32_t NI, 
    float_image_t *shimg[],
    int32_t ix, int32_t iy,
    double score[]
  );
  /* Obtains {score[0..NI-1]} from the actual sharpness images {shimg[0..NI-1]} at 
    pixel in column {ix} and row {iy}. */
    
void mfmi_estimate_Z_and_color_from_scores
  ( int32_t NI, 
    double score[],
    double zFoc[],
    double *zEst_P, 
    int32_t ix, 
    int32_t iy, 
    float_image_t *csimg[],
    frgb_t *clrEst_P
  );
  /* Computes an estimate the value {zEst} of the surface's {Z} at pixel
    {ix,iy} from the scores {score[0..NI-1]} of that pixel in a stack of
    {NI} images and the respective positions {zFoc[0..NI-1]} of the
    focus planes. Then sets {clr[0..NC-1]} to the color of pixel as
    reconstructed from colors of that image stack {csimg[0..NI-1]},
    assuming that each image {csimg[ki]} has the focus plane at
    {zFoc[ki]} and its sharpness scores is {score[ki]}. The results are
    returned in {*zEst_P} and {*clrEst_P}. */
  
void mfmi_fit_quadratics
  ( int32_t NP, 
    double zFoc[],  
    double score[], 
    double *Asc_P,  
    double *Bsc_P,  
    double *Csc_P,
    int32_t NC, 
    float csmp[], 
    double Aclr[],  
    double Bclr[],  
    double Cclr[],
    bool_t verbose
  );
  /* Computes the coefficient of a quadratic
    {score_fit(Z) = Asc*Z^2 + Bsc*Z + Csc} fitted by least squares to the data points 
    {(zFoc[k],score[k])} with {k} in {0..NP-1}. 
    
    Also, for each {ic} in {0..NC-1], computes the coefficient of a quadratic
    {clr_fit[ic] = Aclr[ic]*Z^2 + Bclr[ic]*Z + Cclr[ic]} fitted by least squares 
    to the data points {(zFoc[k],clr[ic][k])} where {clr[ic][k] = csmp[k*NC+ic]}. 
    The number of points must be at least 3. The coefficients {Asc,Bsc,Csc} are returned in 
    {*Asc_P}, {*Bsc_P}, and {*Csc_P}. */

FILE *mfmi_open_text_file(char *outPrefix, char *tag);
  /* Returns the open handle of a file called "{outPrefix}{tag}.txt" for writing. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfmi_options_t *o = mfmi_parse_options(argc, argv);
    
    /* Get the window dimnsions {NW} and sample weights {ws[0..NS-1]}: */
    int32_t NW = 3;
    double *ws = multifok_window_sample_weights(NW);

    /* Create the basis and terms info: */
    int32_t NB;
    double **bas = NULL;
    char **belName = NULL;
    int32_t NT;
    char **termName = NULL; /* Term names. */
    double *wt = NULL;   
    multifok_term_prod_t *prix = NULL;
    mfmi_choose_basis_and_terms(NW, ws, &NB, &bas, &belName, o->noise, &NT, &termName, &wt, &prix, o->outPrefix);

    /* Files with data for plotting, regression, etc: */
    FILE *wr_cmp = mfmi_open_text_file(o->outPrefix, "-cdata");

    /* The multi-focus image stack: */
    double zMin = o->zRange_lo;
    double zMax = o->zRange_hi;
    double zStep = o->zStep;
    int32_t NI = (int32_t)ceil((zMax - zMin + 0.0001*zStep)/zStep);
    float_image_t *csimg[NI]; /* Sythetic images of scene with simulated focus blurring. */
    float_image_t *grimg[NI]; /* Grayscale versions of {csimg[0..NI-1]}. */
    float_image_t *shimg[NI]; /* Specify the actual sharpness of {csimg[k]} at each pixel. */
    double zFoc[NI];         /* The {Z} coordinates of focus planes. */
    double zDep = o->focDepth; /* Depth of focus. */
    
    int32_t NC; /* Number of channels of {csimg[ki]}. */
    int32_t NX, NY; /* Image dimensions. */

    for (int32_t ki = 0; ki < NI; ki++)
      { zFoc[ki] = fmin(zMax, zMin + ki*zStep);
        assert((zFoc[ki] >= zMin) && (zFoc[ki] <= zMax));
        char *zTag = NULL;
        asprintf(&zTag, "-fd%08.4f-z%08.4f", zDep, zFoc[ki]);
        csimg[ki] = multifok_test_read_scene_color_image(o->inPrefix, zTag); 
        if (ki == 0)
          { float_image_get_size(csimg[ki], &NC, &NX, &NY); }
        else
          { float_image_check_size(csimg[ki], NC, NX, NY); }

        /* Convert to grayscale: */
        grimg[ki] = float_image_new(1, NX, NY);
        float_image_map_channels_RGB_to_YUV(csimg[ki], grimg[ki]);

        shimg[ki] = multifok_test_read_sharpness_image(o->inPrefix, zTag); 
        float_image_check_size(shimg[ki], 1, NX, NY);
      } 

    /* Read {timg}, the unblurred true scene colors, for comparison: */
    float_image_t *timg = multifok_test_read_scene_color_image(o->inPrefix, "-sharp");
    float_image_check_size(timg, NC, NX, NY);

    /* Read {azimg}, the unblurred true scene {Z} map, for comparison: */
    float_image_t *azimg = multifok_test_read_zave_image(o->inPrefix, "-sharp");
    float_image_check_size(azimg, 1, NX, NY);
   
    /* Read {dzimg}, the map with deviation of {Z} inside each pixel: */
    float_image_t *dzimg = multifok_test_read_zdev_image(o->inPrefix, "-sharp"); 
    float_image_check_size(dzimg, 1, NX, NY);

    /* Process images: */
    mfmi_process_image_stack
      ( wr_cmp, 
        NI, csimg, grimg, shimg, zFoc,
        timg, azimg, dzimg,
        o->actualSharp,
        NW, ws, o->noise, 
        NB, bas, 
        NT, wt, prix,
        o->outPrefix
      );
      
    fclose(wr_cmp);
      
    return 0;
  }
    
void mfmi_process_image_stack
  ( FILE *wr_cmp,
    int32_t NI, 
    float_image_t *csimg[], 
    float_image_t *grimg[], 
    float_image_t *shimg[],
    double zFoc[],
    float_image_t *timg, 
    float_image_t *azimg, 
    float_image_t *dzimg, 
    bool_t actualSharp,
    int32_t NW,
    double ws[],
    double noise,
    int32_t NB,
    double **bas,
    int32_t NT,
    double wt[],
    multifok_term_prod_t prix[],
    char *outPrefix
  )
  {
    int32_t NC, NX, NY;
    float_image_get_size(csimg[0], &NC, &NX, &NY);
    assert(NC == 3);

    /* Allocate the image with estimated {Z}: */
    fprintf(stderr, "allocating images...\n");
    float_image_t *czimg = float_image_new(1, NX, NY);

    /* Allocate the {Z} error image: */
    float_image_t *ezimg = float_image_new(1, NX, NY);

    /* Allocate the reconstructed (de-blurred) color image: */
    float_image_t *crimg = float_image_new(NC, NX, NY);

    /* Enumerate all windows in in domain and estimate the {Z} coordinate: */
    int32_t HW = NW/2;
    double score[NI];
    double sum_wp_dz = 0;   /* Weighted sum of error (comp {score} minus actual {sharp}). */
    double sum_wp_dz2 = 0;  /* Weighted sum of squared error. */
    double sum_wp = 0;      /* Sum of pixel weights. */
    int32_t NP_scan = 0;    /* Number of pixels processed. */
    int32_t NP_edge = 0;    /* Number of pixels discarded because of high {zdev}. */
    int32_t NP_reg = 0;     /* Number of pixels considered in regression. */
    for (int32_t ix = HW; ix < NX-HW; ix++)
      { for (int32_t iy = HW; iy < NY-HW; iy++) 
          { /* Compute sharpness scores {score[0..NI-1]} of pixel in all images: */
            if (actualSharp)
              { mfmi_get_pixel_actual_sharpness(NI, shimg, ix, iy, score); }
            else
              { mfmi_compute_pixel_scores(NI, grimg, ix, iy, NW, ws, noise, NB, bas, NT, wt, prix, score); }
            
            /* Estimate the surface {Z} and color {clr} at the pixel from {score[]} and {zFoc[]}: */
            double zEst;
            frgb_t clrEst;
            mfmi_estimate_Z_and_color_from_scores(NI, score, zFoc, &zEst, ix, iy, csimg, &clrEst);
            
            /* Save the estimated height {zEst} in the image {czimg}: */
            float_image_set_sample(czimg, 0, ix, iy, (float)zEst);
            
            /* Save the estimated sharp color {clrEst} in the image {crimg}: */
            float_image_set_pixel(crimg, ix, iy, clrEst.c);
            
            /* Compare estimated {Z} with actual {Z}: */
            double zave = float_image_get_sample(azimg, 0, ix, iy);
            double zErr = zEst - zave;
            float_image_set_sample(ezimg, 0, ix, iy, (float)zErr);

            /* Decide whether the pixel is useful for the regression: */
            double zdev = float_image_get_sample(dzimg, 0, ix, iy);
            double zdev_max = 0.05*multifok_scene_ZMAX; /* Ignore pixels which more than this variance. */
            if (zdev > zdev_max) 
              { NP_edge++; }
            else
              { double wp = 1/(1 + zdev*zdev);

                /* Accumulate RMS error: */
                sum_wp_dz += wp*zErr; 
                sum_wp_dz2 += wp*zErr*zErr;
                sum_wp += wp;

                /* Write the data for comparisons: */
                fprintf(wr_cmp, "P%d.%d ", ix, iy);
                fprintf(wr_cmp, " %12.6f ", wp); /* Pixel weight. */
                fprintf(wr_cmp, " %+12.6f ", zave); /* Avg scene height rel to focus plane in pixel. */
                fprintf(wr_cmp, " %12.6f ", zdev); /* Deviation of scene height in pixel. */
                fprintf(wr_cmp, " %+12.6f ", zEst); /* Estimated {Z} of surface in pixel. */
                fprintf(wr_cmp, " %+12.6f ", zErr); /* Error of estimated {Z} in pixel. */
                fprintf(wr_cmp, "\n");

                NP_reg++;
              }
          }
        NP_scan++;
      }
    fprintf(stderr, "ignored %d pixels out of %d for high {zdev}, used %d\n", NP_edge, NP_scan, NP_reg); 

    /* Compute the RMS {Z} estimate error: */
    double avg_err = sum_wp_dz/sum_wp;
    fprintf(stderr, "Average {Z} error (computed minus actual) = %12.4f\n", avg_err);
    double rms_err = sqrt(sum_wp_dz2/sum_wp);
    fprintf(stderr, "RMS {Z} error (computed minus actual) = %12.4f\n", rms_err);

    /* Compute the image {erimg = crimg-timg}, the error in the reconstructed image {crimg}: */
    float_image_t *erimg = float_image_new(NC, NX, NY);
    for (int32_t ix = HW; ix < NX-HW; ix++)
      { for (int32_t iy = HW; iy < NY-HW; iy++) 
          { for (int32_t ic = 0; ic < NC; ic++)
              { float rv = float_image_get_sample(crimg, ic, ix, iy);
                float tv = float_image_get_sample(timg, ic, ix, iy);
                float_image_set_sample(erimg, ic, ix, iy, (float)(rv - tv));
              }
          }
      }
    double sum_var = 0;
    for (int32_t ic = 0; ic < NC; ic++)
      { double avg, dev;
        float_image_compute_sample_avg_dev(erimg, ic, &avg, &dev);
        sum_var += dev*dev;
      }
    double qrms = sqrt(sum_var/NC); /* RMS sample error. */
    fprintf(stderr, "RMS color sample error = %16.8f\n", qrms);

    /* Write the estimated {Z} image {czimg}: */
    multifok_test_write_zavg_image(czimg, outPrefix, "");

    /* Write the estimated {Z} error image {ezimg}: */
    multifok_test_write_Z_error_image(ezimg, outPrefix, "");
      
    /* Write the reconstructed (deblurred) image {crimg}: */
    multifok_test_write_reconstructed_color_image(crimg, outPrefix, "");
      
    /* Write the color error image {erimg}: */
    multifok_test_write_color_error_image(erimg, outPrefix, "");

  }
  
void mfmi_estimate_Z_and_color_from_scores
  ( int32_t NI, 
    double score[],
    double zFoc[],
    double *zEst_P, 
    int32_t ix, 
    int32_t iy, 
    float_image_t *csimg[],
    frgb_t *clrEst_P
  )
  {
    int32_t NC = (int32_t)csimg[0]->sz[0];
    
    bool_t debug = ((ix == 100) && (iy == 100));
    bool_t verbose = debug;
    
    /* Find the largest score: */
    int32_t ki_max = 0;       /* Index of max in {score[0..NI-1]}. */
    for (int32_t ki = 1; ki < NI; ki++)
      { if (score[ki] > score[ki_max]) { ki_max = ki; } }
    
    double zEst;
    frgb_t clrEst;
    if ((ki_max == 0) || (ki_max == NI-1))
      { /* Maximum is at extremity of stack. Just return from that image: */
        zEst = NAN;
      }
    else 
      { /* Fit quadratic functions of {score} and color depending on {zFoc}: */
        int32_t NP; /* Number of data points in regression. */
        int32_t ki0; /* First data point for regression. */
        if ((ki_max == 1) || (ki_max == NI-2))
           { /* Maximum is 1 away form stack end. Fit quadratic to 3 points: */
             ki0 = ki_max - 1;
             NP = 3;
           }
         else
           { /* Maximum is at least 2 away from ends: */
             ki0 = ki_max - 2;
             NP = 5;
           }
         float csmp[NP*NC]; /* Channel {ic} of pixel from {cimh[ki0+kp]} is {csmp[NC*kp + ic]} */
         for (int32_t kp = 0; kp < NP; kp++)
           { float_image_get_pixel(csimg[ki0+kp], ix, iy, &(csmp[NC*kp])); }
         double Asc, Bsc, Csc; /* Coefficients of {score} as function of {zFoc}. */
         double Aclr[NC], Bclr[NC], Cclr[NC]; /* Coefficients of color as function of {zFoc}. */
         if (verbose)
           { fprintf(stderr, "fitting quadratics at pixel (%d,%d):\n", ix, iy);
             fprintf(stderr, "ki_max = %d ki0 = %d NP = %d\n", ki_max, ki0, NP);
           }
         mfmi_fit_quadratics
           ( NP, &(zFoc[ki0]), &(score[ki0]), &Asc, &Bsc, &Csc, 
             NC, csmp, Aclr, Bclr, Cclr,
             verbose
           );
         
         /* Find value {zEst} of {Z} for which fitted {score(Z)} is maximum: */
         if (Asc >= 0)
           { fprintf(stderr, "!! pixel (%d,%d): fitted score is concave\n", ix, iy);
             zEst = NAN;
           }
         else
           { zEst = -Bsc/Asc/2; 
             if ((zEst < zFoc[ki0]) || (zEst > zFoc[ki0 + NP - 1]))
               { fprintf(stderr, "!! pixel (%d,%d): fitted score has faraway max\n", ix, iy);
                 zEst = NAN;
               }
             else
               { /* Compute color from quadratic fit: */
                 for (int32_t ic = 0; ic < NC; ic++) 
                   { clrEst.c[ic] = (float)(Aclr[ic]*zEst*zEst + Bclr[ic]*zEst + Cclr[ic]); } 
               }
           }
       }
     if (isnan(zEst))
       { /* Quadratic fitting failed. Use the maximum: */
         zEst = zFoc[ki_max];
         float_image_get_pixel(csimg[ki_max], ix, iy, clrEst.c);
       }
          
    (*zEst_P) = zEst;
    (*clrEst_P) = clrEst;
  }
   
void mfmi_fit_quadratics
  ( int32_t NP, 
    double zFoc[],
    double score[], 
    double *Asc_P,  
    double *Bsc_P,  
    double *Csc_P,
    int32_t NC, 
    float csmp[], 
    double Aclr[],  
    double Bclr[],  
    double Cclr[],
    bool_t verbose
  )
  {
    if (verbose)
      { fprintf(stderr, "  --------------------------------------------------\n");
        fprintf(stderr, "  assembling and solving the least squares system...\n");
      }
      
    demand(NP >= 3, "can't fit less than 3 points");
    int32_t NX = 3; /* Number of terms in the quadratic. */
    int32_t NF = 1 + NC; /* First func is score, others are color components. */

    auto void gen_point(int32_t k, int32_t nx, double Tk[], int32_t nf, double Fk[], double *WkP);
  
    double U[NX*NF];
    int32_t rank = lsq_fit(NP, NX, NF, gen_point, U, verbose);
    if (rank != 3) { fprintf(stderr, "!! fitting rank = %d\n", rank); }
    
    double Asc = U[0*NF + 0];
    double Bsc = U[1*NF + 0];
    double Csc = U[2*NF + 0];
    
    (*Asc_P) = Asc;
    (*Bsc_P) = Bsc;
    (*Csc_P) = Csc;
    
    for (int32_t ic = 0; ic < NC; ic++)
      { int32_t jf = ic + 1;
        Aclr[ic] = U[0*NF + jf];
        Bclr[ic] = U[1*NF + jf];
        Cclr[ic] = U[2*NF + jf];
      }
        
    if (verbose)
      { fprintf(stderr, "  --------------------------------------------------\n");
        fprintf(stderr, "  fitted quadratic formulas:\n");
        fprintf(stderr, "    score(Z)      = %+.6f * Z^2 %+.6f * Z %+.6f\n", Asc, Bsc, Csc);
        for (int32_t ic = 0; ic < NC; ic++)
          { fprintf(stderr, "    color.c[%d](Z) = %+.6f * Z^2 %+.6f * Z %+.6f\n", ic, Aclr[ic], Bclr[ic], Cclr[ic]); }
        fprintf(stderr, "  --------------------------------------------------\n");
        fprintf(stderr, "  data and results:\n");
        fprintf(stderr, "  scores:\n");
        for (int32_t kp = 0; kp < NP; kp++)
          { double Z = zFoc[kp];
            double scD = score[kp];
            double scF = Asc*Z*Z + Bsc*Z + Csc;
            fprintf(stderr, "    Z = %+12.6f given = %+12.6f fitted = %+12.6f error =  %+12.6f\n", Z, scD, scF, scF-scD);
          }
        for (int32_t ic = 0; ic < NC; ic++)
          { fprintf(stderr, "\n");
            fprintf(stderr, "  color.c[%d]:\n", ic);
            for (int32_t kp = 0; kp < NP; kp++)
              { double Z = zFoc[kp];
                double clrD = csmp[kp*NC + ic];
                double clrF = Aclr[ic]*Z*Z + Bclr[ic]*Z + Cclr[ic];
                fprintf(stderr, "    Z = %+12.6f given = %7.4f fitted = %7.4f error =  %+7.4f\n", Z, clrD, clrF, clrF-clrD);
              }
           }
         fprintf(stderr, "--------------------------------------------------\n");
       }

    return;
    
    void gen_point(int32_t k, int32_t nx, double Tk[], int32_t nf, double Fk[], double *WkP)
      { assert(nx == NX);
        assert(nf == NF);

        /* Independent variables (powers of {zFoc}): */
        double xk = zFoc[k];
        Tk[0] = xk*xk;
        Tk[1] = xk;
        Tk[2] = 1.0;

        /* Dependent variables: */
        Fk[0] = score[k];
        for (int32_t ic = 0; ic < NC; ic ++)
          { int32_t jf = ic + 1;
            Fk[jf] = (double)csmp[k*NC + ic];
          }

        /* Weight: */
        double ak = (k + 0.5)/NP;
        (*WkP) = ak*(1.0 - ak);
      }
  }

void mfmi_get_pixel_actual_sharpness
  ( int32_t NI, 
    float_image_t *shimg[],
    int32_t ix, int32_t iy,
    double score[]
  )
  { for (int32_t ki = 0; ki < NI; ki++)
      { /* Get the "true" sharpness {sharp} at this pixel. */
        double sharp = float_image_get_sample(shimg[ki], 0, ix, iy);
        assert((sharp >= 0.0) && (sharp <= 1.0));
        score[ki] = sharp;
      }
  }

void mfmi_compute_pixel_scores
  ( int32_t NI, 
    float_image_t *grimg[],
    int32_t ix, int32_t iy,
    int32_t NW,
    double ws[],
    double noise,
    int32_t NB,
    double **bas,
    int32_t NT,
    double wt[],
    multifok_term_prod_t prix[],
    double score[]
  )
  {
    int32_t NS = NW*NW; /* Number of samples in window. */
    
    float fsmp[NS];
    double dsmp[NS];
    double coeff[NB]; /* Coefficients of normalized window in be basis. */
    double term[NT]; /* Score terms derived from {coeff[0..NB-1]}. */
    bool_t squared = FALSE;


    for (int32_t ki = 0; ki < NI; ki++)
      { /* Get the samples in the window and normalize them for brightness and contrast: */
        float_image_get_window_samples(grimg[ki], 0,ix,iy, NW, NW, FALSE, fsmp);
        for (int32_t ks = 0; ks < NS; ks++) { dsmp[ks] = fsmp[ks]; }
        double ave, dev;
        multifok_window_normalize_samples(NW, dsmp, ws, noise, &avg, &dev); 

        /* Compute coefficients of window in basis and estimated sharpness score: */

        /* Compute coefficients of window in basis and estimated sharpness score: */
        multifok_basis_compute_coeffs(NW, dsmp, ws, NB, bas, coeff);
        score[ki] = multifok_score_from_basis_coeffs(NW, NB, coeff, NT, prix, wt, term, squared);
      }
  }      

void mfmi_choose_basis_and_terms
  ( int32_t NW,
    double ws[],
    int32_t *NB_P,
    double ***bas_P,
    char ***belName_P,
    double noise, 
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P,
    int32_t *NT_P, 
    double **wt_P, 
    char ***termName_P,
    char *outPrefix
  )
  {
    fprintf(stderr, "choosing and computing the local operator basis...\n");

    char *basisName = "LAPL";
    int32_t termSet = 2;
    multifok_basis_type_t basisType = multifok_basis_type_LAPL;
    
    int32_t NB;
    double **bas = NULL;
    char **belName = NULL;
    bool_t ortho = TRUE;
    multifok_basis_make(NW, ws, basisType, ortho, &NB, &bas, &belName);
    fprintf(stderr, "obtained NB = %d independent basis elements\n", NB);
    multifok_basis_print(stderr, NW, NB, bas, belName);
    multifok_basis_ortho_check(stderr, NW, ws, NB, bas);
      
    multifok_test_write_basis_elem_names(outPrefix, NB, belName);
      
    (*NB_P) = NB;
    (*bas_P) = bas;
    (*belName_P) = belName;
    
    fprintf(stderr, "reading the score terms and weights file...\n");
    char *twfname = NULL;
    asprintf(&twfname, "term-weights/ptBEST-bt%s-trm%d-qsh0-ns%05.2f.txt", basisName, termSet, noise);
    multifok_test_read_term_names_and_weights(twfname, NB, belName, NT_P, termName_P, wt_P, prix_P, TRUE);
    multifok_test_write_term_names_and_weights(outPrefix, *NT_P, *termName_P, *wt_P);
  }

FILE *mfmi_open_text_file(char *outPrefix, char *tag)
  { char *fname = NULL;
    asprintf(&fname, "%s%s.txt", outPrefix, tag);
    FILE *wr = open_write(fname, TRUE);
    return wr;
    free(fname);
  }

mfmi_options_t *mfmi_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    mfmi_options_t *o = notnull(malloc(sizeof(mfmi_options_t)), "no mem");

    argparser_get_keyword(pp, "-inPrefix");
    o->inPrefix = argparser_get_next(pp);
    
    argparser_get_keyword(pp, "-focDepth");
    o->focDepth = argparser_get_next_double(pp, 0.1, 100.0);  
    
    argparser_get_keyword(pp, "-zStep");
    o->zStep = argparser_get_next_double(pp, 0.001, multifok_scene_ZMAX);  

    argparser_get_keyword(pp, "-zRange");
    o->zRange_lo = argparser_get_next_double(pp, 0.0, multifok_scene_ZMAX);  
    o->zRange_hi = argparser_get_next_double(pp, o->zRange_lo, multifok_scene_ZMAX);  

    argparser_get_keyword(pp, "-actualSharp");
    o->actualSharp = argparser_get_next_bool(pp);  

    argparser_get_keyword(pp, "-noise");
    o->noise = argparser_get_next_double(pp, 0.0, 100.0);  

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
