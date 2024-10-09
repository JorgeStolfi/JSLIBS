#define PROG_NAME "test_mfok_make_image"
#define PROG_DESC "test of {multifok_test_image_make.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-10-07 21:18:24 by stolfi */ 
/* Created on 2023-01-05 by J. Stolfi, UNICAMP */

#define test_mfok_make_image_COPYRIGHT \
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
#include <jsmath.h>
#include <jsrandom.h>
#include <r3.h>
#include <ix.h>
#include <i2.h>
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
    char *patternFile;   /* File name of pattern to use to paint objects. */
    char *outPrefix;     /* Prefix for output filenames. */
    /* These parameters are specified in user units: */
    double sceneSize_X;  /* Nominal {X}-span of scene */
    double sceneSize_Y;  /* Nominal {Y}-span of scene */
    double focDepth;     /* Depth of focus. */
    double zRange_lo;    /* {Z} value of lowest frame. */
    double zRange_hi;    /* {Z} value of highest frame. */
    double zStep;        /* {Z} increment between frames. */
  } mfmi_options_t;
  /* Command line parameters. */
  
#define ZMAX multifok_scene_ZMAX
  /* Max scene {Z} coord. */

int32_t main(int32_t argn, char **argv);

mfmi_options_t *mfmi_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */
  
multifok_scene_t *mfmi_make_scene
  ( double WX,
    double WY,
    int32_t NX,
    int32_t NY,
    bool_t rampOnly, 
    double minSep
  );
  /* Generates a test scene with disks, balls, and a back plane, roughly
    spanning the rectanle {[0_WX]×[0_WY]} in {X} and {Y}, and the
    interval {[1.0 _ ZMAX-1.0]} in {Z}.
  
    If {rampOnly} is true, the scene will have just a tilted floor, that
    rises from {Z=box[2].end[0]} at {X=0} to {Z=box[2].end[1]} at
    {X=WX}. The {minSep} is ignored in this case.
    
    If {bak_tilt} is false, the floor will be horizontal a {Z=0}, 
    and there will be a number of textured horizontal disks or balls
    of random radii and centers in front of it.
    
    In this case, the number of objects is chosen by the procedure.
    
    The object radius will range from {rMin} to {rMax} in scene coordinates.
    
    If {minSep} is negative, the {XY} projections of the objects objects
    may overlap , and lie partly outside of the image's domain
    {[0_WX]×[0_WY]}. If {minSep} is non-negative, their {XY} projections
    will be at least {minSep} scene units away from each other and from the edges of
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
    char *outPrefix,
    i2_vec_t plotixy,
    FILE *plt_wr
  );
  /* Creates test images with size {NX,NY} showing the color, nominal
    sharpness, scene {Z} average, and scene {Z} deviation at each pixel
    for the given {scene}, as described in {multifok_scene_images_make}.
    
    Writes the images out to files
    "(outPrefix}-fd{DD.DD}-z{ZZ.ZZ}{tail}" where {DD.DD} is
    {zDep}, {ZZ.ZZ} is {zFoc}, and {tail} is "-cs.ppm" for the color
    image, "-sh.pgm" for the sharpness image, "-az.pgm" for the
    {Z} average image, and "-dz.pgm" for the {Z} deviation image. If 
    {zDep} is infinite then the "-fd{DD.DD}-zf{ZZ.ZZ}" part is replaced by "-sharp".
    
    Also, if {plot_wr} is not {NULL}, writes to it a line with format

      {zFoc} {zDep} {zave[0]} {zdev[0]} {sharp[0]} ... {zave[NQ-1]} {zdev[NQ-1]} {sharp[NQ-1]}

    where {zave[k]}, {zdev[k]}, and {sharp[k]} are the average and deviation 
    of the scene {Z} and the sharpness for the pixel whose indices are {plotixy.e[k]},
    and {NQ} is {plotixy.ne}.  */

float_image_t *mfmi_read_pattern_image(char *fname);
  /* Reads an image from file {fname}, to be used as pattern to paint the scene.
    If the image is in color, converts it to grayscale. */

i2_vec_t mfmi_select_plot_pixels(int32_t NX, int32_t NY, multifok_scene_t *scene);
  /* Seletcts a small number of pixels in the image, biased to where the scene {Z} is close 
    to {ZMAX/2}.  May return 0 pixels.  */
    
FILE *mfmi_open_pixel_plot_data_file(char *outPrefix);
  /* Opens the file {outPrefix-pixplot.txt} for writing. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfmi_options_t *o = mfmi_parse_options(argc, argv);

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
            minSep  = ZMAX/(2*zDep);
          }
        else if (strcmp(o->sceneType, "T") == 0)
          { /* Overlapping objcts: */
            minSep = -1.0;
          }
        else
          { demand(FALSE, "invalid \"-sceneType\""); }
      }

    /* The max number of disks is adjusted to give a reasonable cover fraction in overlap mode: */
    multifok_scene_t *scene = mfmi_make_scene(o->sceneSize_X, o->sceneSize_Y, o->imgSize_X, o->imgSize_Y, rampOnly, minSep);
    
    mfmi_make_images_from_pattern_image(o, scene);
      
    return 0;
  }
  
void mfmi_make_images_from_pattern_image(mfmi_options_t *o, multifok_scene_t *scene)
  { 
    float_image_t *pimg = mfmi_read_pattern_image(o->patternFile);
    int32_t NC_pat, NX_pat, NY_pat;
    float_image_get_size(pimg, &NC_pat, &NX_pat, &NY_pat);
    demand(NC_pat == 1, "pattern image must be grayscale");

    /* Scene size in user units:  */
    double WX_scene = 2*interval_rad(&(scene->dom[0]));
    double WY_scene = 2*interval_rad(&(scene->dom[1]));
    
    /* Scale from pattern pixel {XY} to scene {XY}: */
    double scale_pat = fmin(((double)NX_pat)/WX_scene, ((double)NY_pat)/WY_scene);
     
    auto double photo_pattern(double x, double y, double z, int32_t iobj);
      /* Returns a grayscale pattern value as function of the point with
        scene coordinates {(x,y,z)}. Currently the value of {pimg}
        interpolated at {x,y}, rotated by an angle that depends on
        the parameter {iobj}, typically the index of the scene's object. 
        Used to interpolate between the object's {bg} and {fg}
        colors.
        
        The coordinates are shifted so that the center of {pimg} 
        is centered in scene's domain. */

    mfmi_make_images_from_pattern_function(o, scene, photo_pattern);
    
    return;

    double photo_pattern(double x, double y, double z, int32_t iobj)
      {
        /* Rotate {x,y} proportional to {z}: */
        double rot = iobj*M_PI/5; 
        double cr = cos(rot); 
        double sr = sin(rot);
        double dx_scene = x - 0.5*WX_scene;
        double dy_scene = y - 0.5*WY_scene;
        double xr_pat = 0.5*NX_pat + scale_pat*(cr*dx_scene - sr*dy_scene);
        double yr_pat = 0.5*NY_pat + scale_pat*(sr*dx_scene + cr*dy_scene);
        
        /* Eval the pattern image: */
        int32_t order = 1;  /* Bicubic C1 interpolation. */
        ix_reduction_t red = ix_reduction_PXMIRR;
        double r = float_image_interpolate_sample(pimg, 0, xr_pat, yr_pat, order, red);
        return r;
      }
  }

void mfmi_make_images_from_pattern_function
  ( mfmi_options_t *o, 
    multifok_scene_t *scene, 
    multifok_scene_pattern_t *pattern
  ) 
  {
    /* Image size in pixels: */
    int32_t NX = o->imgSize_X;
    int32_t NY = o->imgSize_Y;

    /* Select the pixels for sharpness x {Z} plot: */
    i2_vec_t plotixy = mfmi_select_plot_pixels(NX, NY, scene);
    FILE *plot_wr = mfmi_open_pixel_plot_data_file(o->outPrefix);
    
    /* Generate the sharp image for reference: */
    int32_t NP_sharp = o->pixSampling;
    int32_t NR_sharp = 1;
    double zDep_sharp = +INF;
    double zFoc_sharp = 0.5*ZMAX;
    mfmi_make_single_Z_images
      ( NX, NY, 
        scene, pattern, 
        zFoc_sharp, zDep_sharp, 
        NP_sharp, NR_sharp, 
        o->outPrefix,
        plotixy, NULL
      );

    /* Generate the blurred image stack: */
    double zDep = o->focDepth;
    int32_t NP_frame = o->pixSampling;
    int32_t NR_min = o->dirSampling;
    double zMin = o->zRange_lo;
    double zMax = o->zRange_hi;
    double zStep = o->zStep;
    double zFoc_frame = zMin;
    while (zFoc_frame <= zMax + 0.0001*zStep)
      { mfmi_make_single_Z_images
        ( NX, NY, 
          scene, pattern, 
          zFoc_frame, zDep, 
          NP_frame, NR_min, 
          o->outPrefix,
          plotixy, plot_wr
        );
        zFoc_frame += zStep;
      }
    fclose(plot_wr);
  }

multifok_scene_t *mfmi_make_scene(double WX, double WY, int32_t NX, int32_t NY, bool_t rampOnly, double minSep)
  {
    /* Compute the box that has to contain all disks and balls: */
    interval_t box[3];
    box[0] = (interval_t){{ 0.0, WX }};
    box[1] = (interval_t){{ 0.0, WY }};
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
        double Wmin = fmin(WX, WY);
        double scale_img = fmax(NX/WX, NY/WY);
        /* The min obj radius must be at least 7 pixels and 1/100 of scene: */
        double rMin = fmax(7.0/scale_img, 0.01*Wmin); /* Min object radius. */
        
        /* The max obj radius must be 1/4 of scene X or Y size. */
        /* Don't worry about Z, will be fitted as needed if spherical: */
        double rMax = fmin(0.75*ZMAX, 0.25*Wmin);       /* Max object radius. */
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
    char *outPrefix,
    i2_vec_t plotixy,
    FILE *plot_wr
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
    
    if (plot_wr != NULL)
      { /* Write data of seected pixels for plot: */
        fprintf(plot_wr, "%10.6f %10.6f", zFoc, zDep);
        for (int32_t k = 0; k < plotixy.ne; k++)
          { i2_t ixy = plotixy.e[k];
            float zave = float_image_get_sample(azimg, 0, ixy.c[0], ixy.c[1]);
            float zdev = float_image_get_sample(dzimg, 0, ixy.c[0], ixy.c[1]);
            float sharp = float_image_get_sample(shimg, 0, ixy.c[0], ixy.c[1]);
            fprintf(plot_wr, "  %10.6f %10.6f %12.8f", zave, zdev, sharp);
          }
        fprintf(plot_wr, "\n");
      }

    multifok_test_write_scene_color_image(csimg, outPrefix, tag);
    multifok_test_write_zave_image(azimg, outPrefix, tag);
    multifok_test_write_zdev_image(dzimg, outPrefix, tag);
    multifok_test_write_sharpness_image(shimg, outPrefix, tag);
    multifok_test_write_pix_sel_image(csimg, plotixy, outPrefix, tag);

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
    float_image_t *pat = float_image_read_pnm_named(fname, isMask, gama, bias, yup, warn, verbose);
    int32_t NC_pat, NX_pat, NY_pat;
    float_image_get_size(pat, &NC_pat, &NX_pat, &NY_pat);
    if (NC_pat == 3)
      { /* Convert color image to grayscale: */
        float_image_t *grimg = float_image_new(1, NX_pat, NY_pat);
        float_image_map_channels_RGB_to_YUV(pat, grimg);
        float_image_free(pat);
        pat = grimg;
        NC_pat = 1;
      }
    demand(NC_pat == 1, "pattern image must be RGB or grayscale");
    return pat;
  }

i2_vec_t mfmi_select_plot_pixels(int32_t NX, int32_t NY, multifok_scene_t *scene)
  {
    bool_t debug = TRUE;
    
    fprintf(stderr, "selecting plot pixels...\n");
    int32_t NO = scene->NO; /* Number of objects inscene. */
    int32_t NQ_exp = 10; /* Number of sample pixels to generate.*/
    
    /* Select the pixels: */
    i2_vec_t ixy = i2_vec_new(NQ_exp);
    int32_t NQ = 0; /* Pixels actually selected. */
    
    auto void selpix(double x_scene, double y_scene);
      /* Converts the point {(x_scene,y_scene)} from scene coordinates to
        a pair of pixel indices, and appends that pair to the list of sample pixels. */
    
    if (NO > 0)
      { /* Pick {NQ_exp} object centers, preferably flat: */
        /* Count the number of disks {NO_disk}: */
        int32_t NO_disk = 0;
        for (int32_t ko = 0; ko < NO; ko++)
          { multifok_scene_object_t *objk = &(scene->objs[ko]);
            if (objk->flat) { NO_disk++; }
          }
        bool_t disk_only =  (NO_disk >= NQ_exp);

        int32_t NO_rem = (disk_only ? NO_disk : NO); /* Number of candidate objects remaining. */
        for (int32_t ko = 0; ko < NO; ko++)
          { if ((NQ >= NQ_exp) || (NO_rem == 0)) { break; }
            multifok_scene_object_t *objk = &(scene->objs[ko]);
            if ((! disk_only) || objk->flat)
              { int32_t NQ_rem = NQ_exp - NQ;  /* Number of pixels still to be selected. */
                double prpick = (NO_rem <= NQ_rem ? 1.0 : ((double)NQ_rem)/NO_rem);
                if (drandom() < prpick) 
                  { double x_scene = objk->ctr.c[0];
                    double y_scene = objk->ctr.c[1];
                    selpix(x_scene, y_scene);
                  }
                NO_rem--;
              }
          }
      }
    else if (! scene->flatFloor)
      { /* Scene has only a slanted floor: */
        double WX_scene = 2*interval_rad(&(scene->dom[0]));
        double WY_scene = 2*interval_rad(&(scene->dom[1]));
        for (int32_t kq = 0; kq < NQ_exp; kq++)
          { double fr = (kq + 0.5)/NQ_exp;
            double x_scene = fr*WX_scene;
            double y_scene = fr*WY_scene;
            selpix(x_scene, y_scene);
          }
      }
    
    i2_vec_trim(&ixy, NQ);
    return ixy;
      
    void selpix(double x_scene, double y_scene) 
      { 
        r2_t p_img = multifok_test_scene_to_pixel(x_scene, y_scene, scene->dom, NX, NY);
        int32_t ix = (int32_t)floor(p_img.c[0] + 0.5);
        int32_t iy = (int32_t)floor(p_img.c[1] + 0.5);
        if ((ix >= 0) && (ix < NX) && (iy >= 0) && (iy < NY))
          { if (debug) { fprintf(stderr, "  selected %d,%d\n", ix, iy); }
            i2_vec_expand(&ixy, NQ);
            ixy.e[NQ] = (i2_t){{ ix, iy }};
            NQ++;
          }
      }
  }    
  
FILE *mfmi_open_pixel_plot_data_file(char *outPrefix)  
  {
    char *fname = NULL;
    asprintf(&fname, "%s-pixplot.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    free(fname);
    return wr;
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

    argparser_get_keyword(pp, "-sceneSize");
    o->sceneSize_X = (int32_t)argparser_get_next_int(pp, 30, 4096);
    o->sceneSize_Y = (int32_t)argparser_get_next_int(pp, 30, 4096);

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
    o->zRange_hi = argparser_get_next_double(pp, o->zRange_lo, multifok_scene_ZMAX + o->zStep);  

    argparser_get_keyword(pp, "-patternFile");
    o->patternFile = argparser_get_next(pp);

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
