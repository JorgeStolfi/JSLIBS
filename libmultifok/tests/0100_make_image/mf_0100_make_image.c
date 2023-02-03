#define PROG_NAME "mf_0100_make_image"
#define PROG_DESC "test of {multifok_test_image_make.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-02-01 06:23:13 by stolfi */ 
/* Created on 2023-01-05 by J. Stolfi, UNICAMP */

#define mf_0100_make_image_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

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
#include <ix.h>
#include <frgb.h>
#include <float_image.h>
#include <float_image_waves.h>
#include <float_image_interpolate.h>
#include <float_image_write_pnm.h>
#include <float_image_read_pnm.h>
#include <float_image_map_channels.h>

#include <multifok_test.h>
#include <multifok_scene.h>

typedef struct mfmi_options_t 
  { int32_t imgSize_X;   /* Image width. */
    int32_t imgSize_Y;   /* Image height. */
    char *sceneType;     /* Scene type: "R", "F", "T", etc. */
    int32_t pixSampling; /* Number of sampling points per axis and pixel. */
    int32_t dirSampling; /* Min number of aperture rays per point. */
    double focDepth;     /* Depth of focus. */
    double zRange_lo;    /* {Z} value of lowest frame. */
    double zRange_hi;    /* {Z} value of highest frame. */
    double zStep;        /* {Z} increment between frames. */
    char *patternFile;   /* File name of pattern to use to paint objects. */
    char *outPrefix;     /* Prefix for output filenames. */
  } mfmi_options_t;
  /* Command line parameters. */

int32_t main(int32_t argn, char **argv);

mfmi_options_t *mfmi_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */
  
multifok_scene_t *mfmi_make_scene(int32_t NX, int32_t NY, bool_t rampOnly, double minSep);
  /* Generates a test scene with disks, balls, and a back plane, roughly spanning the 
    image domain {[0_NX]×[0_NY]} in {X} and {Y}, and the interval {[1.0 _ZMAX-1.0]} in {Z},
    where  {ZMAX} is {multifok_scene_ZMAX}.
  
    If {rampOnly} is true, the scene will have just a tilted floor, 
    that rises from {Z=box[2].end[0]} at {X=0} to {Z=box[2].end[1]} at {X=NX}.  The {minSep} is
    ignored in this case.
    
    If {bak_tilt} is false, the floor will be horizontal a {Z=0}, 
    and there will be a number of textured horizontal disks or balls
    of random radii and centers in front of it.
    
    In this case, the number of objects is chosen by the procedure.
    
    If {minSep} is negative, the {XY} projections of the objects objects
    may overlap , and lie partly outside of the image's domain
    {[0_NX]×[0_NY]}. If {minSep} is non-negative, their {XY} projections
    will be at least {minSep} away from each other and from the edges of
    that rectangle. In any case, the objects {Z} projections will be
    entirely contained in the interval {box[2]}. */

void mfmi_make_images_from_pattern_image(mfmi_options_t *o, multifok_scene_t *scene);
  /* Makes blurred and  sharp images of the {scene}, colored with a given image. */

void mfmi_make_images_from_pattern_function
  ( mfmi_options_t *o, 
    multifok_scene_t *scene, 
    multifok_scene_pattern_t *pattern
  );
  /* Makes blurred and  sharp images of the {scene}, colored with the given pattern 
    function. */

void mfmi_make_single_Z_images
  ( int32_t NX, 
    int32_t NY, 
    multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    double zFoc,
    double zDep, 
    int32_t NP,
    int32_t NR_min,
    char *outPrefix
  );
  /* Creates test images with size {NX,NY} showing the color, nominal
    sharpness, scene {Z} average, and scene {Z} deviation at each pixel
    for the given {scene}, as described in {multifok_scene_images_make}.
    
    Writes the images out to files
    "(outPrefix}-fd{DD.DD}-z{ZZ.ZZ}{tail}" where {DD.DD} is
    {zDep}, {ZZ.ZZ} is {zFoc}, and {tail} is "-cs.ppm" for the color
    image, "-sh.pgm" for the sharpness image, "-az.pgm" for the
    {Z} average image, and "-dz.pgm" for the {Z} deviation image. If 
    {zDep} is infinite then the "-fd{DD.DD}-zf{ZZ.ZZ}" part is replaced by "-sharp". */

float_image_t *mfmi_read_pattern_image(char *fname);
  /* Reads an image from file {fname}, to be used as pattern to paint the scene.
    If the image is in color, converts it to grayscale. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfmi_options_t *o = mfmi_parse_options(argc, argv);

    int32_t NX = o->imgSize_X;
    int32_t NY = o->imgSize_Y;

    double zDep = o->focDepth;

    /* Create the test scene: */
    bool_t rampOnly; /* If true the scene is just a rampOnly, else it has balls and disks. */
    double minSep; /* Min {XY} sep of non-overlapping objs. */
    if (strcmp(o->sceneType, "R") == 0)
      { /* Tilted floor surface only: */
        rampOnly = TRUE;
        minSep = NAN;
      }
    else 
      { /* Disk, balls, etc: */
        rampOnly = FALSE;
        if (strcmp(o->sceneType, "F") == 0)
          { /* Non-overlapping objects: */
            minSep  = multifok_scene_ZMAX/(2*zDep);
          }
        else if (strcmp(o->sceneType, "T") == 0)
          { /* Overlapping objcts: */
            minSep = -1.0;
          }
        else
          { demand(FALSE, "invalid \"-sceneType\""); }
      }

    /* The max number of disks is adjusted to give a reasonable cover fraction in overlap mode: */
    multifok_scene_t *scene = mfmi_make_scene(NX, NY, rampOnly, minSep);
    
    mfmi_make_images_from_pattern_image(o, scene);
      
    return 0;
  }
  
void mfmi_make_images_from_pattern_image(mfmi_options_t *o, multifok_scene_t *scene)
  { 
    
    float_image_t *pimg = mfmi_read_pattern_image(o->patternFile);
    int32_t NC_pat, NX_pat, NY_pat;
    float_image_get_size(pimg, &NC_pat, &NX_pat, &NY_pat);
    demand(NC_pat == 1, "pattern image must be grayscale");

    int32_t NX_img = o->imgSize_X;
    int32_t NY_img = o->imgSize_Y;
     
    auto double photo_pattern(double x, double y, double z, int32_t iobj);
      /* Returns a grayscale pattern value as function of the point {(x,y,z)}.
        Currently the value of {pimg} interpolated at {x,y}, rotated
        by an angle that depends on {iobj}.  Used to interpolate between
        the object's {bg} and {fg} colors.  
        
        The coordinates are shifted so that the center of {pimg} 
        is centered in the new image. */

    mfmi_make_images_from_pattern_function(o, scene, photo_pattern);
    
    return;

    double photo_pattern(double x, double y, double z, int32_t iobj)
      {
        /* Rotate {x,y} proportional to {z}: */
        double rot = iobj*M_PI/5; 
        double cr = cos(rot); 
        double sr = sin(rot);
        double dx = x - 0.5*NX_img, dy = y - 0.5*NY_img;
        double xr = 0.5*NX_pat + cr*dx - sr*dy;
        double yr = 0.5*NY_pat + sr*dx + cr*dy;
        
        /* Eval the pattern image: */
        int32_t order = 1;  /* Bicubic C1 interpolation. */
        ix_reduction_t red = ix_reduction_PXMIRR;
        double r = float_image_interpolate_sample(pimg, 0, xr, yr, order, red);
        return r;
      }
  }

void mfmi_make_images_from_pattern_function
  ( mfmi_options_t *o, 
    multifok_scene_t *scene, 
    multifok_scene_pattern_t *pattern
  ) 
  {
    int32_t NX = o->imgSize_X;
    int32_t NY = o->imgSize_Y;

    double zDep = o->focDepth;
    
    int32_t NP = o->pixSampling;
    int32_t NR_min = o->dirSampling;
    double zMin = o->zRange_lo;
    double zMax = o->zRange_hi;
    double zStep = o->zStep;
    double zFoc = zMin;

    while (zFoc <= zMax + 0.0001*zStep)
      { mfmi_make_single_Z_images(NX, NY, scene, pattern, zFoc, zDep, NP, NR_min, o->outPrefix);
        zFoc += zStep;
      }
    int32_t NR_sharp = 1;
    double zDep_sharp = +INF;
    double zFoc_sharp = 0.5*multifok_scene_ZMAX;
    mfmi_make_single_Z_images(NX, NY, scene, pattern, zFoc_sharp, zDep_sharp, NP, NR_sharp, o->outPrefix);
  }

multifok_scene_t *mfmi_make_scene(int32_t NX, int32_t NY, bool_t rampOnly, double minSep)
  {
    /* Compute the box that has to contain all disks and balls: */
    double ZMAX = multifok_scene_ZMAX;
    interval_t box[3];
    box[0] = (interval_t){{ 0.0, (double)NX }};
    box[1] = (interval_t){{ 0.0, (double)NY }};
    box[2] = (interval_t){{ 1.0, ZMAX-1.0 }};
    
    bool_t verbose = TRUE;
    bool_t flatFloor = (! rampOnly);
    multifok_scene_t *scene = multifok_scene_new(box, flatFloor, verbose);
    if (rampOnly)
      { assert(scene->NO == 0);
        fprintf(stderr, "generated a ramp-only scene\n");
      }
    else
      { /* Disk size range: */
        double rMin = 2.5; /* Min object radius (pixels). */
        double rMax = 3.5; /* Max object radius (pixels). */
        multifok_scene_throw_objects(scene, rMin, rMax, minSep, verbose);
        char *olapx = (minSep < 0.0 ? "" : "non-overlapping ");
        fprintf(stderr, "generated a scene with %d %sdisks\n", scene->NO, olapx);
      }
    return scene;
  }

void mfmi_make_single_Z_images
  ( int32_t NX, 
    int32_t NY, 
    multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    double zFoc,
    double zDep, 
    int32_t NP,
    int32_t NR_min,
    char *outPrefix
  )
  {
    fprintf(stderr, "generating images for zDep = %12.6f zFoc = %12.6f\n", zDep, zFoc);

    /* Filename tag: */
    char *tag = NULL;
    if (zDep == +INF)
      { asprintf(&tag, "-sharp"); }
    else
      { asprintf(&tag, "-fd%05.2f-zf%05.2f", zDep, zFoc); }
    
    /* Generate the images: */
    float_image_t *csimg = NULL;
    float_image_t *shimg = NULL;
    float_image_t *azimg = NULL;
    float_image_t *dzimg = NULL;
    multifok_test_images_make(NX, NY, scene, pattern, zFoc, zDep, NP, NR_min, &csimg, &shimg, &azimg, &dzimg);

    multifok_test_write_scene_color_image(csimg, outPrefix, tag);
    multifok_test_write_zavg_image(azimg, outPrefix, tag);
    multifok_test_write_zdev_image(dzimg, outPrefix, tag);
    multifok_test_write_sharpness_image(shimg, outPrefix, tag);

    free(tag);
    float_image_free(csimg);
    float_image_free(azimg);
    float_image_free(dzimg);
    float_image_free(shimg);
  }

#define mfmi_image_gamma 1.000
  /* Assumed encoding gamma of input and output images. */

#define mfmi_image_bias 0.0327
  /* Assumed encoding bisa of input and output images. */

float_image_t *mfmi_read_pattern_image(char *fname)
  {
    bool_t isMask = FALSE;
    double gama = mfmi_image_gamma;
    double bias = mfmi_image_bias;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = FALSE;
    float_image_t *img = float_image_read_pnm_named(fname, isMask, gama, bias, yup, warn, verbose);
    int32_t NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    if (NC == 3)
      { /* Convert color image to grayscale: */
        float_image_t *grimg = float_image_new(1, NX, NY);
        float_image_map_channels_RGB_to_YUV(img, grimg);
        float_image_free(img);
        img = grimg;
        NC = 1;
      }
    demand(NC == 1, "pattern image must be RGB or grayscale");
    return img;
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

    argparser_get_keyword(pp, "-sceneType");
    o->sceneType = argparser_get_next_non_keyword(pp);  
    
    argparser_get_keyword(pp, "-pixSampling");
    o->pixSampling = (int32_t)argparser_get_next_int(pp, 1, 9999);  

    argparser_get_keyword(pp, "-dirSampling");
    o->dirSampling = (int32_t)argparser_get_next_int(pp, 1, 9999);  

    argparser_get_keyword(pp, "-focDepth");
    o->focDepth = argparser_get_next_double(pp, 0.1, 1.0e200);  

    argparser_get_keyword(pp, "-zStep");
    o->zStep = argparser_get_next_double(pp, 0.001, multifok_scene_ZMAX);  

    argparser_get_keyword(pp, "-zRange");
    o->zRange_lo = argparser_get_next_double(pp, 0.0, multifok_scene_ZMAX);  
    o->zRange_hi = argparser_get_next_double(pp, o->zRange_lo, multifok_scene_ZMAX);  

    argparser_get_keyword(pp, "-patternFile");
    o->patternFile = argparser_get_next(pp);

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
