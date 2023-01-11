#define PROG_NAME "mf_0100_make_image"
#define PROG_DESC "test of {multifok_test_image_make.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-01-10 20:58:46 by stolfi */ 
/* Created on 2012-01-25 by J. Stolfi, UNICAMP */

#define mf_0100_make_image_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

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
#include <r3.h>
#include <frgb.h>
#include <float_image.h>
#include <float_image_waves.h>
#include <float_image_write_pnm.h>

#include <multifok_test.h>

#define SCENE_DEPTH 10.0
  /* Total depth of simulated scene. */

typedef struct mfmi_options_t 
  { int32_t imgSize_X;   /* Image width. */
    int32_t imgSize_Y;   /* Image height. */
    int32_t rays;        /* Number of rays per pixel. */
    bool_t dense;        /* If true, generates overlapping disks. */
    double focDepth;     /* Depth of focus. */
    double rPix;         /* Nominal radius of pixel. */
    double zRange_lo;    /* {Z} value of lowest frame. */
    double zRange_hi;    /* {Z} value of highest frame. */
    double zStep;        /* {Z} increment between frames. */
    char *outPrefix;     /* Prefix for output filenames. */
  } mfmi_options_t;
  /* Command line parameters. */

int32_t main(int32_t argn, char **argv);

mfmi_options_t *mfmi_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */
  
multifok_test_scene_t *mfmi_make_scene(int32_t NX, int32_t NY, double minSep, double rMin, double rMax);
  /* Generates a test scene with disks of random colors, random radii in
    {[rMIn _ rMax]}, and random {Z} coordinates between 0 and
    {SCENE_DEPTH}.
    
    The number of disks is chosen by the procedure.
    
    If {minSep} is negative, the disks may overlap and touch the edges
    of the image. If {minSep} is non-negative, the disks will be at
    least {minSep} away from each other in their {XY} projecttion and
    from the edges of the image. In this second case, the scene may have
    fewer than {ND} disks. */

void mfmi_make_single_Z_images
  ( int32_t NX, 
    int32_t NY, 
    multifok_test_scene_t *scene,
    multifok_test_pattern_t *pattern,
    int32_t NR, 
    double rPix,
    double zDep, 
    double zFoc,
    char *outPrefix
  );
  /* Creates test images with size {NX,NY} showing the color and nominal focus blur radius
    for the given {scene}, using {NR} rays per pixel.
    Assumes {zDep} as the nominal depth of focus and {zFoc} as the
    {Z} position of the focus plane.  The disks and backplane will be 
    colored with mixes of the respecive {bg} and {fg} colors, in proportions
    defined by the {pattern} function. */

void mfmi_write_color_image(float_image_t *cimg, char *outPrefix, char *tag);
  /* Writes the color image {cimg} to file "{outPrefix}{tag}-c.ppm". */

void mfmi_write_radius_image(float_image_t *rimg, double rPix, double zDep, char *outPrefix, char *tag);
  /* Writes the blur radius image {rimg} to file "{outPrefix}{tag}-r.pgm".
    The {zDep} is used to scale the sample values in a consistent way. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfmi_options_t *o = mfmi_parse_options(argc, argv);
    
    /* Chose the wave parameters for the {pattern} function: */
    fprintf(stderr, "choosing frequencies...\n");
    int32_t NF = 5;
    int32_t NS = 0;
    double amp[NF];
    double fx[NF];
    double fy[NF]; 
    double phase[NF];
    bool_t verbose = TRUE;
    float_image_waves_pick(NF, amp, fx, fy, phase, verbose);
    
    /* Choosing a squash amplitude: */
    double sum_a2 = 0;
    for (int32_t kf = 0; kf < NS + NF; kf++) { sum_a2 += amp[kf]*amp[kf]; }
    double amp_rms = sqrt(sum_a2); /* Root-mean-sum amplitude. */
    double squash = 0.5*amp_rms;
    
    auto double pattern(double x, double y, double z);
      /* Returns a grayscale pattern value as function of the point {(x,y,z)}.
        Currently a combination of the above waves on {(x,y)}, rotated 
        by an angle proportional to {z}  */

    int32_t NX = o->imgSize_X;
    int32_t NY = o->imgSize_Y;

    double rPix = o->rPix;
    double zDep = o->focDepth;
    double zMin = o->zRange_lo;
    double zMax = o->zRange_hi;
    double zStep = o->zStep;

    /* Choose suitable parameters: */
    double rMin = 0.5;  /* Min disk radius (pixels). */
    double rMax = 10.0; /* Max disk radius (pixels). */
    double minSep = (o->dense ? -1 : SCENE_DEPTH/(2*zDep));
    /* The max number of disks is adjusted to give a reasonable cover fraction in dense mode: */
    multifok_test_scene_t *scene = mfmi_make_scene(NX, NY, minSep, rMin, rMax);
    fprintf(stderr, "generated a scene with %d %sdisks\n", scene->ND, (o->dense ? "" : "non-overlapping "));
    
    double zFoc = zMin;
    while (zFoc <= zMax + 0.0001*zStep)
      { mfmi_make_single_Z_images(NX, NY, scene, pattern, o->rays, rPix, zDep, zFoc, o->outPrefix);
        zFoc += zStep;
      }
    int32_t NR_sharp = 30;
    double zDep_sharp = +INF;
    double zFoc_sharp = (zMin + zMax)/2;
    mfmi_make_single_Z_images(NX, NY, scene, pattern, NR_sharp, rPix, zDep_sharp, zFoc_sharp, o->outPrefix);
    return 0;

    double pattern(double x, double y, double z)
      {
        /* Rotate {x,y} proportional to {z}: */
        double rot = z*M_PI/5; 
        double cr = cos(rot); 
        double sr = sin(rot);
        double xr = + cr*x - sr*y;
        double yr = + sr*x + cr*y;
        
        /* Eval the waves: */
        double r = float_image_waves_eval(xr, yr, NF, amp, fx, fy, phase);
        
        /* Map {r} to {[-1 _ +1]}: */
        r = r/squash;
        r = r/hypot(1, r);
        
        /* Map {r} to {[0_1]}: */
        r = (r + 1.0)/2;
        return r;
      }

  }

multifok_test_scene_t *mfmi_make_scene(int32_t NX, int32_t NY, double minSep, double rMin, double rMax)
  {
    /* Compute the box that has to contain all disks: */
    double mrg = (minSep >= 0 ? minSep : 0); /* Margin to leave all around the image. */
    interval_t box[3];
    box[0] = (interval_t){{ mrg, NX - mrg }};
    box[1] = (interval_t){{ mrg, NY - mrg }};
    box[2] = (interval_t){{ 0, SCENE_DEPTH }};
    
    /* Compute the target number of disks: */
    double rr0 = rMin + (minSep >= 0 ? 0.5*minSep : 0); /* Min radius disk occup. */
    double rr1 = rMax + (minSep >= 0 ? 0.5*minSep : 0); /* Max radius disk occup. */
    double aDsk = M_PI*(rr1*rr1*rr1 - rr0*rr0*rr0)/(rr1 - rr0)/3; /* Average disk area. */
    double wx = box[0].end[1] - box[0].end[0];
    double wy = box[1].end[1] - box[1].end[0];
    double aBox = wx*wy; /* Total useful area of box. */
    int32_t ND_in = (int32_t)ceil(aBox/aDsk);
    
    /* Generate the disks: */
    bool_t verbose = TRUE;
    multifok_test_scene_t *scene = multifok_test_scene_throw(box, ND_in, rMin, rMax, minSep, verbose);
    return scene;
  }

void mfmi_make_single_Z_images
  ( int32_t NX, 
    int32_t NY, 
    multifok_test_scene_t *scene,
    multifok_test_pattern_t *pattern,
    int32_t NR, 
    double rPix,
    double zDep, 
    double zFoc,
    char *outPrefix
  )
  {
    fprintf(stderr, "generating images for zDep = %12.6f zFoc = %12.6f\n", zDep, zFoc);

    /* Filename tag: */
    char *tag = NULL;
    if (zDep == +INF)
      { asprintf(&tag, "-sharp-z%08.4f", zFoc); }
    else
      { asprintf(&tag, "-fd%08.4f-z%08.4f", zDep, zFoc); }
    
    /* Generate the images: */
    float_image_t *cimg = NULL;
    float_image_t *rimg = NULL;
    multifok_test_images_make(NX, NY, scene, pattern, rPix, zFoc, zDep, NR, &cimg, &rimg);

    mfmi_write_color_image(cimg, outPrefix, tag);
    mfmi_write_radius_image(rimg, rPix, zDep, outPrefix, tag);

    free(tag);
    float_image_free(cimg);
    float_image_free(rimg);
  }

#define mfmi_image_gamma 1.000
  /* Assumed encoding gamma of input and output images. */

#define mfmi_image_bias 0.0327
  /* Assumed encoding bisa of input and output images. */

void mfmi_write_color_image(float_image_t *cimg, char *outPrefix, char *tag)
  {
    assert(cimg->sz[0] == 3);
    char *fname = NULL;
    asprintf(&fname, "%s%s-c.ppm", outPrefix, tag);
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    double gamma = mfmi_image_gamma;
    double bias = mfmi_image_bias;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = TRUE;
    float_image_write_pnm_named(fname, cimg, isMask, gamma, bias, yup, warn, verbose);
    free(fname);
  }

void mfmi_write_radius_image(float_image_t *rimg, double rPix, double zDep, char *outPrefix, char *tag)
  {
    assert(rimg->sz[0] == 1);
    
    /* Determine the max blur radius assuming that {zFoc} is in {[0_SCENE_DEPTH]}: */
    double rbMax = hypot(rPix, SCENE_DEPTH/zDep);

    /* Remap image values to {[0_1]}: */
    float vMin = 0.0; /* Force image 0 to be radius 0. */
    float vMax = (float)rbMax;
    float_image_rescale_samples(rimg, 0, vMin, vMax, 0.0, 1.0);

    char *fname = NULL;
    asprintf(&fname, "%s%s-r.pgm", outPrefix, tag);
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    double gamma = mfmi_image_gamma;
    double bias = mfmi_image_bias;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = TRUE;
    float_image_write_pnm_named(fname, rimg, isMask, gamma, bias, yup, warn, verbose);
    free(fname);
  }

mfmi_options_t *mfmi_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, "HELP!");
    argparser_set_info(pp, "KNOW THYSELF");
    argparser_process_help_info_options(pp);
    
    mfmi_options_t *o = notnull(malloc(sizeof(mfmi_options_t)), "no mem");

    argparser_get_keyword(pp, "-imgSize");
    o->imgSize_X = (int32_t)argparser_get_next_int(pp, 30, 4096);
    o->imgSize_Y = (int32_t)argparser_get_next_int(pp, 30, 4096);

    argparser_get_keyword(pp, "-rays");
    o->rays = (int32_t)argparser_get_next_int(pp, 1, 9999);  

    argparser_get_keyword(pp, "-dense");
    o->dense = argparser_get_next_bool(pp);  

    argparser_get_keyword(pp, "-focDepth");
    o->focDepth = argparser_get_next_double(pp, 0.1, 1.0e200);  

    argparser_get_keyword(pp, "-rPix");
    o->rPix = argparser_get_next_double(pp, 0.1, 3.0);  

    argparser_get_keyword(pp, "-zStep");
    o->zStep = argparser_get_next_double(pp, 0.001, SCENE_DEPTH);  

    argparser_get_keyword(pp, "-zRange");
    o->zRange_lo = argparser_get_next_double(pp, 0.0, SCENE_DEPTH);  
    o->zRange_hi = argparser_get_next_double(pp, o->zRange_lo, SCENE_DEPTH);  

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
