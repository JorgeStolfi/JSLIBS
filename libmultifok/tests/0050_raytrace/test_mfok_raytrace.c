#define PROG_NAME "test_mfok_raytrace"
#define PROG_DESC "test of {multifok_sampling.h} and {multifok_raytrace.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-02-05 15:14:24 by stolfi */
/* Created on 2023-01-05 by J. Stolfi, UNICAMP */

#define test_mfok_raytrace_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>

#include <bool.h>
#include <frgb.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <r2.h>
#include <r3.h>
#include <i2.h>
#include <affirm.h>
#include <argparser.h>
#include <float_image.h>
#include <float_image_test.h>
#include <float_image_interpolate.h>
#include <float_image_read_gen.h>
#include <float_image_write_gen.h>

#include <multifok_sampling.h>
#include <multifok_raytrace.h>
#include <multifok_frame.h>

#define zFoc_MAX 80.0
  /* Maximum {zFoc} to consider. */

#define zDep_DEFAULT 4.0
  /* The depth-of-focus to use. */
  

/* For each frame, the program writes a file "{frameFolder}/pixel-rays-{XXXX}-{YYYY}.txt",
where {XXXX} and {YYYY} are the indices of the pixel {o.debigPixel} 
formated as "%04d".   The file will will have one line for each ray that was 
cast for that pixel, as per {multifok_raytrace_write_ray_data_INFO}.

The program will also write a single file
"{stackFolder}/pixel-data-{XXXX}-{YYYY}.txt" with one line for each frame,
showing the computed data for that pixel, as per
multifok_raytrace_write_pixel_data_INFO}. */

typedef struct mfrs_options_t
  { int32_t imageSize_X;    /* Image width. */
    int32_t imageSize_Y;    /* Image height. */
    char *imageType;        /* Image type: "bullsex", "bullsqr", "noise01", etc. */
    uint32_t pixSampling;   /* Determines the sampling points per axis and pixel. */
    uint32_t dirSampling;   /* Min rays per sampling point, or 0 for sharp view. */
    i2_t debugPixel;        /* Pixel to debug. */
    char *stackFolder;         /* Directory where to put output frame folders. */
  } mfrs_options_t;
  /* Command line parameters. */

int32_t main(int32_t argn, char **argv);

mfrs_options_t *mfrs_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

void mfrs_make_and_write_image
  ( int32_t NX,
    int32_t NY,
    float_image_test_generator_t *fimg,
    uint32_t HS,
    uint32_t KR_min,
    double zDep,
    double zFoc,
    multifok_raytrace_debug_pred_t *debug_pixel,
    multifok_raytrace_report_ray_proc_t report_ray,
    multifok_raytrace_report_pixel_proc_t report_pixel,
    char *frameFolder,
    bool_t verbose
  );
  /* Creates a monochrme test image {oimg} with size {NX,NY} showing the
    color at each pixel for of a focus-blurred version of the image
    {fimg}.  Writes the image to "{frameFolder}/img-blur.png".

    using {(2*HS+1)^2} subsampling points at each pixel aof
    {oimg} and about {3*HR^2} rays through each sampling point.

    The output image {oimg} and the procedural image {fimg}
    are assumed to have the same {Z} axis and parallel {X} and {Y}
    axes, with {Z=0} of the former corresponding to {Z=zFoc} in the latter.
    The {X} and {Y} image coordinates for both images are assumed
    to range in {[0 _ NX] × [0 _ NY]}.

    If {HS} is zero, each pixel of {oimg} will have only one sampling
    point, at its center. Otherwise the sampling points will extend
    outside that pixel, with 2D Hann window weighting.

    If {HR} is zero, there will be only one ray per sampling point,
    directed straight down. Otherwise, besides that ray there will be
    one or more other oblique rays; their directions will be chosen so
    that at distance {zDep/2} below the plane of {oimg} plane their
    spread (defined as the RMS value of their distance from vertical)
    will be 1.

    As a special case, if {zFoc} is zero or {zDep} is {+INF}, then
    {HR} is ignored, and the scene view {sVal} will be just
    {fimg} with antialiasing by the sub-pixel sampling. */

void mfrs_write_image(float_image_t *oimg, char *frameFolder);
  /* Writes the image {oimg} to file
    "{frameFolder}/img-blur.png". */
    
    
float_image_t *mfrs_read_pattern_image(char *patName);
  /* Reads an image from file "in/{patName}.png". */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfrs_options_t *o = mfrs_parse_options(argc, argv);
    
    int32_t NX = o->imageSize_X;
    int32_t NY = o->imageSize_Y;
    
    float_image_t *img_pat = NULL;
    
    auto void pattern_from_image(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[]);
      /* A "procedural" image function that actually samples the image {img_pat} at 
      the point {p} with bilinear interpolation. */

    /* Choose the test image to be blurred: */
    float_image_test_generator_t *fimg;
    if (strcmp(o->imageType, "bullsex") == 0)
      { fimg = &float_image_test_gen_bullsex; }
    else if (strcmp(o->imageType, "bullsqr") == 0)
      { fimg = &float_image_test_gen_bullsqr; }
    else if (strcmp(o->imageType, "noise01") == 0)
      { img_pat = mfrs_read_pattern_image(o->imageType);
        fimg = &pattern_from_image;
      }
    else
      { demand(FALSE, "unrecognized image type"); }
 
    uint32_t HS = o->pixSampling;
    uint32_t KR_min = o->dirSampling;

    uint32_t kf = 0; /* Images generated so far. */

    /* Define the output directory {stackFolder}: */
    char *stackFolder = o->stackFolder;

    i2_t *iPixDeb = &(o->debugPixel);
    
    char *fname_pixels = jsprintf("%s/pixel-data-%04d-%04d.txt", stackFolder, iPixDeb->c[0], iPixDeb->c[1]);
    FILE *wr_pixels = open_write(fname_pixels, TRUE);
    free(fname_pixels);

    auto bool_t debug_pixel(i2_t *iPix);
      /* A procedure suitable for the {debug_pixel} argument of
        {multifok_raytrace_make_frame}. It returns {TRUE} iff
        the pixel {iPix} is {iPixDeb}. */

    auto void report_pixel
      ( i2_t *iPix, r3_t *pCtr, double zFoc, double zDep,
        double shrp, double hAvg, double hDev, 
        r3_t *sNrm, int32_t NC, float sVal[]
      );
      /* A procedure suitable for the {report_pixel} argument of
        {multifok_raytrace_make_frame}.  It writes the data to the file {wr_pixel}.
        It will be called only when {iPix} is {iPixDeb}. */

    double zDep = zDep_DEFAULT;
    double zFoc = 0.0;
    while (zFoc < zFoc_MAX + 0.0000001)
      { bool_t verbose_frame = (kf == 0) || (kf == 1);

        char *frameFolder = jsprintf("%s/zf%08.4f-df%08.4f", stackFolder, zFoc, zDep);
        mkdir(frameFolder, 0755);

        char *fname_rays = jsprintf("%s/pixel-rays-%04d-%04d.txt", frameFolder, iPixDeb->c[0], iPixDeb->c[1]);
        FILE *wr_rays = open_write(fname_rays, TRUE);
        free(fname_rays);

        auto void report_ray
          ( i2_t *iPix, i2_t *iSmp, double step, double wSmp,
            r3_t *pRay, r3_t *dRay, double wRay,
            r3_t *pHit, double hHit, double vBlr,
            r3_t *sNrm, int32_t NC, float sVal[]
          );
          /* A procedure suitable for the {report_ray} argument of
            {multifok_raytrace_make_frame}. It writes the data to the 
            file {wr_rays}. It will be called only when {iPix} is 
            equal to {iPixDeb}. */

        mfrs_make_and_write_image
          ( NX, NY, fimg, HS, KR_min, zFoc, zDep, 
            &debug_pixel, &report_ray, &report_pixel, frameFolder, verbose_frame
          );
        free(frameFolder);

        if (wr_rays != NULL) { fclose(wr_rays); }

        zFoc = (zFoc == 0 ? zDep/4.0 : sqrt(2)*zFoc);
        kf++;

        void report_ray
          ( i2_t *iPix, i2_t *iSmp, double step, double wSmp,
            r3_t *pRay, r3_t *dRay, double wRay,
            r3_t *pHit, double hHit, double vBlr,
            r3_t *sNrm, int32_t NC, float sVal[]
          )
          { demand((iPix->c[0] == iPixDeb->c[0]) && (iPix->c[1] == iPixDeb->c[1]), "wrong pixel");
            double pixSize = 1.0;
            multifok_raytrace_write_ray_data
              ( wr_rays, iPix, pixSize, 
                iSmp, step, wSmp, pRay, dRay, wRay, pHit, hHit, vBlr, sNrm, NC, sVal
              );
          }
      }
    
    fclose(wr_pixels);

    return 0;
    
    void pattern_from_image(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[])
      { assert(img_pat != NULL);
        double dfs[NC];
        float_image_interpolate_pixel(img_pat, p->c[0], p->c[1], 1, ix_reduce_mode_PXMIRR, dfs);
        for (int32_t k = 0; k < NC; k++) { fs[k] = (float)dfs[k]; }
      }

    bool_t debug_pixel(i2_t *iPix)
      { bool_t debug_col = iPix->c[0] == iPixDeb->c[0];
        bool_t debug_row = iPix->c[1] == iPixDeb->c[1];
        bool_t debug = (debug_col && debug_row);
        return debug;
      }
    
    void report_pixel
      ( i2_t *iPix, r3_t *pCtr, double zFoc, double zDep,
        double shrp, double hAvg, double hDev, 
        r3_t *sNrm, int32_t NC, float sVal[]
      )
      { double pixSize = 1.0;
        demand((iPix->c[0] == iPixDeb->c[0]) && (iPix->c[1] == iPixDeb->c[1]), "wrong pixel");
        multifok_raytrace_write_pixel_data
          ( wr_pixels, iPix, pixSize, pCtr, zFoc, zDep,
            shrp, hAvg, hDev, sNrm, NC, sVal
          );
      }
  }

void mfrs_make_and_write_image
  ( int32_t NX,
    int32_t NY,
    float_image_test_generator_t *fimg,
    uint32_t HS,
    uint32_t KR_min,
    double zFoc,
    double zDep,
    multifok_raytrace_debug_pred_t *debug_pixel,
    multifok_raytrace_report_ray_proc_t report_ray,
    multifok_raytrace_report_pixel_proc_t report_pixel,
    char *frameFolder,
    bool_t verbose
  )
  {
    fprintf(stderr, "generating blurred image zFoc = %12.6f pixels\n", zFoc);

    multifok_sampling_t *samp = multifok_sampling_choose(HS, KR_min, verbose);
    uint32_t NS = samp->NS;  /* Number of pixel sub-sampling points. */
    uint32_t KR = samp->KR;  /* Number of pixel sub-sampling points. */
    /* Paranoia: */
    assert(NS >= 1);
    assert(KR >= 1);

    auto void trace_ray(r3_t *pR, r3_t *dR, bool_t debug_loc, r3_t *pHit_P, r3_t *sNrm_P, int32_t NC_loc, float sVal[]);
      /* A ray-tracing func suitable for {multifok_raytrace_make_frame}.
        The procedure assumes that the scene is the input image {fimg},
        covering the rectangle {[0 _ NX] × [0 _ NY]} on the plane {Z=0}. */

    auto void map_point(r2_t *p2_img, r3_t *p3_scene);
      /* An image-to-scene coordinate conversion function, suitable for
        the {map} argument of {multifok_raytrace_make_frame}. The point
        {p2_img} is assumed to be in the 2D image coordinate system,
        inside the output image domain or just outside it.
        It assumes that the image domain is the rectalgle
        {[0 _ NX] × [0 _ NY]} located at {Z = Zfoc}. */

    int32_t NC = 1;
    r3_t dRef = (r3_t){{ 0, 0, -1 }};

    multifok_frame_t *fr = multifok_raytrace_make_frame
      ( NC, NX, NY, trace_ray, map_point,
        &dRef, zFoc, zDep, samp,
        verbose, debug_pixel, report_ray, report_pixel
      );

    mfrs_write_image(fr->sVal, frameFolder);

    multifok_sampling_free(samp);

    multifok_frame_free(fr);

    return;

    void trace_ray(r3_t *pR, r3_t *dR, bool_t debug_loc, r3_t *pHit_P, r3_t *sNrm_P, int32_t NC_loc, float sVal[])
      { assert(NC_loc == NC);
        assert(dR->c[2] < 0);

        double zscale = -pR->c[2]/dR->c[2];
        r2_t q = (r2_t){{ pR->c[0] + dR->c[0]*zscale, pR->c[1] + dR->c[1]*zscale }};
        r3_t sNrm = (r3_t){{0,0,0}};
        if ((q.c[0] < 0) || (q.c[0] > NX) || (q.c[1] < 0) || (q.c[1] > NY))
          { for (uint32_t ic = 0; ic < NC; ic++) { sVal[ic] = 0.0; } }
        else
          { fimg(&q, NC, NX, NY, sVal);
            for (int32_t j = 0; j < 3; j++) { sNrm.c[j] = (j < NC ? sVal[j] : 1.0e-6); }
            (void)r3_dir(&sNrm, &sNrm);
          }
        r3_t pHit = (r3_t){{ q.c[0], q.c[1], 0.0 }};
        (*pHit_P) = pHit;
        (*sNrm_P) = sNrm;
      }

    void map_point(r2_t *p2_img, r3_t *p3_scene)
      {
        (*p3_scene) = (r3_t){{ p2_img->c[0], p2_img->c[1], zFoc }};
      }
  }

void mfrs_write_image(float_image_t *oimg, char *frameFolder)
  { char *fname = jsprintf("%s/img-blur.png", frameFolder);
    image_file_format_t ffmt = image_file_format_PNG;
    bool_t yUp = TRUE;
    double gammaEnc = 1.0;
    double bias = 0.0;
    bool_t verbose = TRUE;
    float_image_write_gen_named(fname, oimg, ffmt, yUp, 0.0,1.0, gammaEnc,bias, verbose);
    free(fname);
  }
  
float_image_t *mfrs_read_pattern_image(char *patName)
  { 
    char *fileName = jsprintf("in/%s.png", patName);
    bool_t yUp = TRUE;
    double gammaDec, bias;
    bool_t verbose = TRUE;
    float_image_t *img = float_image_read_gen_named
      ( fileName, image_file_format_PNG, yUp, 
        0.0, 1.0, NULL, &gammaDec, &bias, 
        verbose
      );
    free(fileName);
    return img;
  }

mfrs_options_t *mfrs_parse_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, "HELP!");
    argparser_set_info(pp, "KNOW THYSELF");
    argparser_process_help_info_options(pp);

    mfrs_options_t *o = notnull(malloc(sizeof(mfrs_options_t)), "no mem");

    argparser_get_keyword(pp, "-imageSize");
    o->imageSize_X = (int32_t)argparser_get_next_int(pp, 30, 4096);
    o->imageSize_Y = (int32_t)argparser_get_next_int(pp, 30, 4096);
    argparser_get_keyword(pp, "-imageType");
    o->imageType = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-pixSampling");
    o->pixSampling = (uint32_t)argparser_get_next_int(pp, 0, 9999);

    argparser_get_keyword(pp, "-dirSampling");
    o->dirSampling = (uint32_t)argparser_get_next_int(pp, 1, 9999);

    argparser_get_keyword(pp, "-debugPixel");
    o->debugPixel.c[0] = (int32_t)argparser_get_next_int(pp, 0, o->imageSize_X-1);
    o->debugPixel.c[1] = (int32_t)argparser_get_next_int(pp, 0, o->imageSize_Y-1);

    argparser_get_keyword(pp, "-stackFolder");
    o->stackFolder = argparser_get_next_non_keyword(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
