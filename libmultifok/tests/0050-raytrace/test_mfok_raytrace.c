#define PROG_NAME "test_mfok_raytrace"
#define PROG_DESC "test of {multifok_sampling.h} and {multifok_raytrace.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-10-22 08:34:09 by stolfi */ 
/* Created on 2023-01-05 by J. Stolfi, UNICAMP */

#define test_mfok_raytrace_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
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
#include <r2.h>
#include <r3.h>
#include <i2.h>
#include <affirm.h>
#include <argparser.h>
#include <float_image.h>
#include <float_image_test.h>
#include <float_image_write_gen.h>

#include <multifok_sampling.h>
#include <multifok_raytrace.h>
#include <multifok_frame.h>

#define zFoc_MAX 80.0
  /* Maximum {zFoc} to consider. */
  
#define zDep_DEFAULT 4.0
  /* The depth-of-focus to use. */

typedef struct mfrs_options_t 
  { int32_t imageSize_X;   /* Image width. */
    int32_t imageSize_Y;   /* Image height. */
    char *imageType;       /* Image type: "bullsex", "bullsqr", etc. */
    int32_t pixSampling;   /* Determines the sampling points per axis and pixel. */
    int32_t dirSampling;   /* Min number of aperture rays per sampling point. */
  } mfrs_options_t;
  /* Command line parameters. */
  
int32_t main(int32_t argn, char **argv);

mfrs_options_t *mfrs_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */
  
void mfrs_make_and_write_image
  ( int32_t NX, 
    int32_t NY, 
    float_image_test_generator_t *fimg,
    int32_t HS,
    int32_t HR,
    double zDep,
    double zFoc,
    char *frameDir,
    bool_t verbose
  );
  /* Creates a monochrme test image {oimg} with size {NX,NY} showing the
    color at each pixel for of a focus-blurred version of the image
    {fimg}.  Writes the image to "{frameDir}/img-blur.png".
    
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
 
void mfrs_write_image(float_image_t *oimg, char *frameDir);
  /* Writes the image {oimg} to file
    "{frameDir}/img-blur.png". */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfrs_options_t *o = mfrs_parse_options(argc, argv);
    
    int32_t NX = o->imageSize_X;
    int32_t NY = o->imageSize_Y;
    
    /* Choose the test image to be blurred: */
    float_image_test_generator_t *fimg;
    if (strcmp(o->imageType, "bullsex") == 0)
      { fimg = &float_image_test_gen_bullsex; }
    else if (strcmp(o->imageType, "bullsqr") == 0)
      { fimg = &float_image_test_gen_bullsqr; }
    else
      { demand(FALSE, "unrecognized image type"); }

    mkdir("out", 0755); /* "-rwxr-xr-x" */

    int32_t HS = o->pixSampling;
    int32_t KR = o->dirSampling;
    double zDep = zDep_DEFAULT;

    /* Define the output directory {frameDir}: */
    char *stackDir = NULL; asprintf(&stackDir, "out/img-%s-%04dx%04d-hs%02d-kr%02d", o->imageType, NX, NY, HS, KR);
    mkdir(stackDir, 0755);

    int32_t kf = 0; /* Images generated so far. */
    double zFoc = 0.0;
    while (zFoc < zFoc_MAX + 0.0000001)
      { bool_t verbose = (kf == 0) || (kf == 1);

        char *frameDir = NULL; asprintf(&frameDir, "%s/zf%08.4f-df%08.4f", stackDir, zFoc, zDep); 
        mkdir(frameDir, 0755);

        mfrs_make_and_write_image(NX, NY, fimg, HS, KR, zFoc, zDep, frameDir, verbose);
        free(frameDir);
        
        zFoc = (zFoc == 0 ? zDep/4.0 : sqrt(2)*zFoc);
        kf++;
      }
    
    free(stackDir);
    return 0;
  }
  
void mfrs_make_and_write_image
  ( int32_t NX, 
    int32_t NY,
    float_image_test_generator_t *fimg,
    int32_t HS,
    int32_t KR,
    double zFoc,
    double zDep,
    char *frameDir,
    bool_t verbose
  )
  { 
    fprintf(stderr, "generating blurred image zFoc = %12.6f pixels\n", zFoc);
    
    /* Pixel sub-sampling points: */
    int32_t NS;  /* Number of pixel sub-sampling points. */
    r2_t *uSmp;
    double *wSmp;
    multifok_sampling_choose_pixel_sampling_points_and_weights(HS, &NS, &uSmp, &wSmp, verbose);

    /* Choose the relative ray tilts and weights: */
    int32_t NR; /* Actual number of aperture rays. */
    r2_t *tRay;
    double *wRay;
    multifok_sampling_choose_ray_tilts_and_weights(KR, NS, &NR, &tRay, &wRay, verbose);
    /* Make sure that {NR} is divisible by {NS}: */
    assert(NR >= NS);
    assert((NR % NS) == 0);
    
    i2_t iDeb = (i2_t){{ ix_DEBUG, iy_DEBUG }};

    auto void trace_ray(r3_t *pR, r3_t *dR, bool_t debug_loc, r3_t *pHit_P, int32_t NC_loc, float colr[]);
      /* A ray-tracing func suitable for {multifok_raytrace_make_frame}.
        The procedure assumes that the scene is the input image {fimg},
        covering the rectangle {[0 _ NX] × [0 _ NY]} on the plane {Z=0}. */
    
    auto void map_point(r2_t *p2_img, r3_t *p3_scene);
      /* An image-to-scene coordinate conversion fucntion, suitable for
        the {map} argument of {multifok_raytrace_make_frame}. The point
        {p2_img} is assumed to be in the 2D image coordinate system,
        inside the output image domain or just outside it.
        It assumes that the image domain is the rectalgle
        {[0 _ NX] × [0 _ NY]} located at {Z = Zfoc}. */

    auto bool_t debug_pix(i2_t *iPix);
      /* A procedure suitable for the {debug_pix} argument of
        {multifok_raytrace_make_frame}. It returns {TRUE} iff 
        the pixel {iPix} is {iDeb} or adjacent to it. */

    FILE *wr_ray = NULL;
    i2_t iPix_wr = (i2_t){{ -1, -1 }}; /* Pixel for which {wr_ray} was opened. */
        
    auto void report_ray
      ( i2_t *iPix, r2_t *pSmp, double wSmp,
        r3_t *pRay, r3_t *dRay, double wRay, 
        r3_t *pHit, double hHit, double vBlr
      );
      /* Type of a procedure suitable for the {report_ray} argument of
        {multifok_raytrace_make_frame}.
        
        The procedure prints the ray data on {stderr} and also writes it
        to the file {wr_pix}, saving {iPix} to {iPix_wr}. The file
        {wr_pix} is opened if needed, and closed and re-opened whenever
        {iPix} differs from {iPix_wr}. */
    
    int32_t NC = 1;
    r3_t dRef = (r3_t){{ 0, 0, -1 }};
    
    multifok_frame_t *fr = multifok_raytrace_make_frame
      ( NC, NX, NY, trace_ray, map_point, 
        &dRef, zFoc, zDep, NS, uSmp, wSmp, NR, tRay, wRay, 
        verbose, &debug_pix, &report_ray 
      );
     
    mfrs_write_image(fr->sVal, frameDir);
    
    if (wr_ray != NULL) { fclose(wr_ray); }

    free(tRay);
    free(wRay);

    multifok_frame_free(fr);
    
    return;
  
    void trace_ray(r3_t *pR, r3_t *dR, bool_t debug_loc, r3_t *pHit_P, int32_t NC_loc, float colr[])
      { assert(NC_loc == NC);
        assert(dR->c[2] < 0);
        
        double zscale = -pR->c[2]/dR->c[2];
        r2_t q = (r2_t){{ pR->c[0] + dR->c[0]*zscale, pR->c[1] + dR->c[1]*zscale }};
        if ((q.c[0] < 0) || (q.c[0] > NX) || (q.c[1] < 0) || (q.c[1] > NY))
          { for (int32_t ic = 0; ic < NC; ic++) { colr[ic] = 0.0; } }
        else
          { fimg(&q, NC, NX, NY, colr); }
        r3_t pHit = (r3_t){{ q.c[0], q.c[1], 0.0 }};
        (*pHit_P) = pHit;
      }
      
    void map_point(r2_t *p2_img, r3_t *p3_scene)
      {
        (*p3_scene) = (r3_t){{ p2_img->c[0], p2_img->c[1], zFoc }};
      }
      
    bool_t debug_pix(i2_t *iPix)
      { bool_t debug_col = (abs(iPix->c[0] - iDeb.c[0]) <= 1);
        bool_t debug_row = (abs(iPix->c[1] - iDeb.c[1]) <= 1);
        bool_t debug = (debug_col && debug_row);
        return debug;
      }
    
    void report_ray
      ( i2_t *iPix, r2_t *pSmp, double wSmp,
        r3_t *pRay, r3_t *dRay, double wRay, 
        r3_t *pHit, double hHit, double vBlr
      )
      { if ((wr_ray != NULL) && ((iPix->c[0] != iPix_wr.c[0]) || (iPix->c[1] != iPix_wr.c[1])))
          { fclose(wr_ray);
            iPix_wr = (i2_t){{ -1, -1 }};
            wr_ray = NULL;
          }
        if (wr_ray == NULL)
          { char *fname_ray = NULL;
            asprintf(&fname_ray, "%s/pixel-rays-%04d-%04d.txt", frameDir, iPix->c[0], iPix->c[1]);
            wr_ray = open_write(fname_ray, TRUE);
            free(fname_ray);
            iPix_wr = *iPix;
          }
        fprintf(stderr, "        ");
        double pixSize = 1.0;
        multifok_raytrace_show_ray_data(stderr, iPix, pixSize, pSmp, wSmp, pRay, dRay, wRay, pHit, hHit, vBlr);
        multifok_raytrace_write_ray_data(wr_ray, iPix, pixSize, pSmp, wSmp, pRay, dRay, wRay, pHit, hHit, vBlr);
      }
  }
 
void mfrs_write_image(float_image_t *oimg, char *frameDir)
  { char *fname = NULL;
    asprintf(&fname, "%s/img-blur.png", frameDir);
    image_file_format_t ffmt = image_file_format_PNG;
    double gammaEnc = 1.0;
    double bias = 0.0;
    bool_t verbose = TRUE;
    float_image_write_gen_named(fname, oimg, ffmt, 0.0,1.0, gammaEnc,bias, verbose);
    free(fname);
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
    o->pixSampling = (int32_t)argparser_get_next_int(pp, 0, 9999);  

    argparser_get_keyword(pp, "-dirSampling");
    o->dirSampling = (int32_t)argparser_get_next_int(pp, 1, 9999);  

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
