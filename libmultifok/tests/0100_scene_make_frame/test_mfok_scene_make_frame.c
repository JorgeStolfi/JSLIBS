#define PROG_NAME "test_mfok_scene_make_frame"
#define PROG_DESC "test of {multifok_test_image_make.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-16 12:49:29 by stolfi */ 
/* Created on 2023-01-05 by J. Stolfi, UNICAMP */

#define test_mfok_scene_make_frame_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>

#include <affirm.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <argparser.h>
#include <argparser_geo.h>
#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <ix_reduce.h>
#include <r3.h>
#include <i2.h>
#include <frgb.h>
#include <float_image.h>
#include <float_image_waves.h>
#include <float_image_interpolate.h>
#include <float_image_write_gen.h>
#include <float_image_read_gen.h>
#include <float_image_map_channels.h>

#include <multifok_test.h>
#include <multifok_scene.h>
#include <multifok_stack.h>
#include <multifok_frame.h>
#include <multifok_image.h>
#include <multifok_raytrace.h>
#include <multifok_scene_tree.h>
#include <multifok_scene_raytrace.h>
#include <multifok_scene_make_frame.h>

#define PROG_HELP \
  "  Please help yourself!"

#define PROG_INFO \
  "OUTPUTS\n" \
  "\n" \
  PROG_OUTPUT_INFO

typedef struct mfmi_options_t 
  { int32_t imageSize_X;       /* Image width. */
    int32_t imageSize_Y;       /* Image height. */
    char *sceneType;           /* Scene type: "R", "F", "T", etc. */
    uint32_t pixSampling;      /* Number of sampling points per axis and pixel. */
    uint32_t dirSampling;      /* Min number of aperture rays per point. */
    char *patternFile;         /* File name of pattern to use to paint objects. */
    r3_t *lightDir;            /* Light source direction, or {NULL}. */
    double ambient;            /* Ambient light fraction. */
    char *stackDir;            /* Parent folder for frames of the stack. */
    /* These parameters are specified in user units: */
    interval_t sceneSize[3];   /* Nominal coordinate ranges of scene. */
    double dephOfFocus;        /* Depth of focus. */
    double focusHeight_lo;     /* Lowest focus plane {Z}. */
    double focusHeight_hi;     /* Highest focus plane {Z}. */
    double focusHeight_step;   /* {Z} increment between focus planes. */
  } mfmi_options_t;
  /* Command line parameters. */
  
#define PROG_OUTPUT_INFO \
  "  RUN DIRECTORY\n" \
  "    The program writes all its output files to a specified folder \"{stackDir}\".\n" \
  "\n" \
  "  FRAME DIRECTORIES\n" \
  "    " multifok_FRAME_DIR_INFO "\n" \
  "\n" \
  "  FRAME IMAGE FILES\n" \
  "\n" \
  "    The run directory will contain one frame sub-folder {frameDir} for each frame of the stack.  " multifok_FRAME_FILES_INFO "\n" \
  "\n" \
  "  PIXEL V FOCUS PLOT FILES\n" \
  "\n" \
  "    The program also writes to the run directory \"{stackDir}\" zero or more" \
  " files \"pixel-data-{XXXX}-{YYYY}.txt\" where {XXXX} and {YYYY} are the column and" \
  " row indices {ix,iy} of a pixel that was selected for detailed debugging.  Each of" \
  " these files contain one line for each frame, excluding the sharp one.  Each line has" \
  " fields\n" \
  "\n" \
  "     \"{kf} {zFoc[kf]} {zDep[kf]} {hAvg[kf]} {hDev[kf]} {shrp[kf]} {sVal[kf][0]} .. {sVal[kf][NC-1]}\n" \
  "\n" \
  " where {kf} is a frame index (from 0), {zFoc[kf]} and {zDep[kf]} are the frame's" \
  " nominal {zFoc} and {zDep} parameters, and the other fields are the sample values of" \
  " pixel {ix,it} of the corresponding image from that frame.  These files can be used" \
  " to plot the variation of those pixel properties as a function of the" \
  " parameters {zFoc} and/or {zDep}.\n" \
  "\n" \
  "  RAY DATA FILES\n" \
  "\n" \
  "    The program will also write in each frame directory \"{stackDir}/{framedir}\" zero" \
  " or more files \"pixel-rays-{XXXX}-{YYYY}.txt\" where {XXXX} and {YYYY} are the column" \
  " and row indices {ix,iy} of a selected debigging pixel, as above.  This file contains on" \
  "e line for each ray that was cast in order to compute the values of the frame images" \
  " at that pixel. \n" \
  "    " multifok_raytrace_write_ray_data_INFO "\n" \
  "\n" \
  "    The fields {xPix}, {yPix}, and {pixSize} are the same on all lines of each" \
  " file. Since the in-focus plane is horizontal and {zFoc} is its {Z} scene" \
  " coordinate, the field {hHit} is the same as {pHit.z}, and {vBlr} is" \
  " {(pHit.x-pRay.x)^2 +(pHit.y-pRay.y)^2}."
  
#define ot_FLAT multifok_scene_object_type_FLAT
#define ot_RAMP multifok_scene_object_type_RAMP
#define ot_DISK multifok_scene_object_type_DISK
#define ot_BALL multifok_scene_object_type_BALL
#define ot_CONE multifok_scene_object_type_CONE

int32_t main(int32_t argn, char **argv);

mfmi_options_t *mfmi_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */
  
multifok_scene_t *mfmi_make_scene
  ( char *sceneType,
    interval_t dom[],
    int32_t NX,
    int32_t NY,
    double zDep,
    bool_t verbose
  );
  /* Generates a test scene with disks, balls, and a back plane, roughly
    spanning the box {dom[0]×dom[1]×dom[2]} in scene coordinates.
    
    If {sceneType} is "R", the scene will have no disks or balls, just a
    titled floor --- a ramp that rises from {Z=0} at {X=0} to
    {Z=WZ} at {X=WX}.
    
    If {sceneType} is "F" or "T", the floor will be flat at {Z=0.0},
    and there will be some number of disks and balls at random
    {XY} coordinates and contained between {Z=0} and {Z=WZ}.
    
    Specifically, if {sceneType} is "F", the {XY} projections of the
    disks and balls will be disjoint, even when out of focus (assuming
    that the depth of focus is {zDep}). If {sceneType} is "T", the {XY}
    projections of the disks and balls will often overlap. */

multifok_stack_t *mfmi_make_and_write_stack
  ( char *stackDir,
    int32_t NX,
    int32_t NY, 
    multifok_scene_t *scene,
    char *patternFile,
    r3_t *light_dir,
    double ambient,
    double zFoc_min,
    double zFoc_max,
    double zFoc_step,
    double zDep,
    uint32_t HS,
    uint32_t KR,
    uint32_t NQ,
    i2_t iDeb[]
  );
  /* Makes a stack of {NI} frames {fr[0..NI-1]} with simulated views of the {scene},
    where {NI} is {o->NI + 1}.  
    
    The scene is texturized with a monochrome image to be read
    from the given  {patternFile}.
    
    The first {NI-1} frames {fr[0..NI-2]} will be blurred, with the {Z}
    coordinate {fr[i].zFoc} of the in-focus plane varying from
    {zFoc_min} by {zFoc_step}, up to but not exceeding {zFoc_max}. The
    depth of focus {fr[i].zDep} wil be {zDep}. For the meaning of {HS}
    and {KR}, see {multifok_scene_make_frame}.

    The last frame {fr[NI-1]} will be sharp, with {fr[NI-1].zDep = +INF}
    and {fr[NI-1].zFoc} arbitrary A single vertical ray will be used for
    each image sampling point.
    
    Also writes debuging information about the rays used to compute the
    pixels {iDeb[0..NQ-1]}, to files called
    "{stackDir}/{frameDir}/pixel-rays-{XXXX}-{YYYY}.txt", where
    {XXXX,YYYY} are the pixel indices formatted as "%04d". Here
    {stackDir} is as explained in {PROG_INFO}. */

multifok_stack_t *mfmi_make_and_write_stack_from_pattern_function
  ( char *stackDir,
    int32_t NX,
    int32_t NY,
    multifok_scene_t *scene, 
    multifok_scene_raytrace_pattern_t *pattern,
    r3_t *light_dir,
    double ambient,
    double zFoc_min,
    double zFoc_max,
    double zFoc_step,
    double zDep,
    uint32_t HS,
    uint32_t KR,
    uint32_t NQ,
    i2_t iDeb[] 
  );
  /* Same as {mfmi_make_stack}, execpt that the scene is colorized
    with the given pattern function. */

multifok_frame_t *mfmi_make_and_write_frame
  ( char *frameDir,
    int32_t NX, 
    int32_t NY, 
    multifok_scene_t *scene,
    multifok_scene_tree_t *tree,
    multifok_scene_raytrace_pattern_t *pattern,
    r3_t *light_dir,
    double ambient,
    double zFoc,
    double zDep, 
    uint32_t HS,
    uint32_t KR, 
    bool_t verbose,
    uint32_t NQ,
    i2_t iDeb[]
  );
  /* Creates test images with size {NX,NY} showing the color, nominal
    blurring indicator, scene {Z} average, and scene {Z} deviation at each pixel
    for the given {scene}, as described in {multifok_scene_images_make}.
    
    As a special case, if {zDep} is infinite then the scene view {sVal}
    will be sharp everywhere, and {shrp} will be 1.0.
    
    Writes the images to "{frameDir}/{tag}.png" were {tag} is "sVal",
    "hAvg", "hDev", "shrp". Also writes pixel ray data files
    "{frameDir}/pixel-rays-{XXXX}-{YYYY}.txt" */

void mfmi_write_pixel_profiles
  ( multifok_stack_t *stack,
    uint32_t NQ,
    i2_t iDeb[],
    multifok_scene_t *scene,
    char *stackDir
  );
  /* Given a list of {NQ} of pixels {iDeb[0..NQ-1]} in the 
    image domain, writes a file "{stackDir}/pixel-data-{XXXX}-{YYYY}.txt"
    for each pixel, where {XXXX} an {YYYY} are the pixel column and row indices.
    
    Each file will have up to {NI = stack.NI} lines with format

      "{ki} {zFoc} {zDep}  {hAvg} {hDev} {shrp} {sVal[0]} ...  {sVal[NC-1]}"

    where {ki} is a frame number in {0..NI-1} and {hAvg}, and 
    {hDev}, {shrp}, and {sVal[0..NC-1]} are the values of
    that pixel in the images {hAvg} {hDev}, {shrp}, and {sVal} of
    {stack.frame[ki]}. Note that {zDep} may be {+INF}.  */
    
void mfmi_select_debug_pixels
  ( int32_t NX,
    int32_t NY,
    multifok_scene_t *scene,
    uint32_t NQ_max,
    uint32_t *NQ_P,
    i2_t **iDeb_P
  );
  /* Selects a certain number {NQ <= NQ_max} of pixel index pairs {iDeb[0..NQ-1]}.
    Returns {NQ} and {iDeb} in {*NQ_P} and {ideb_P}. 
    
    The pixels will be preferrably selected over the scene's objects, if there 
    are any. */

float_image_t *mfmi_read_pattern_image(char *fname);
  /* Reads an image from file {fname}, to be used as pattern to paint the scene.
    If the image is in color, converts it to grayscale. */
    
FILE *mfmi_open_pixel_plot_data_file(char *stackDir);
  /* Opens the file "{stackDir}/pixplot.txt" for writing. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfmi_options_t *o = mfmi_parse_options(argc, argv);
    
    int32_t NX = o->imageSize_X;
    int32_t NY = o->imageSize_Y;

    /* Create the test scene: */
    bool_t verbose_scene = FALSE;
    multifok_scene_t *scene = mfmi_make_scene
      ( o->sceneType, o->sceneSize,
        NX, NY, 
        o->dephOfFocus, verbose_scene
      );
      
    /* Select the debugging pixels: */
    uint32_t NQ_max = 10;
    uint32_t NQ;
    i2_t *iDeb;
    mfmi_select_debug_pixels(NX, NY, scene, NQ_max, &NQ, &iDeb);

    multifok_stack_t *stack = mfmi_make_and_write_stack
      ( o->stackDir,
        NX, NY, scene, o->patternFile, o->lightDir, o->ambient, 
        o->focusHeight_lo, o->focusHeight_hi, o-> focusHeight_step,
        o->dephOfFocus, o->pixSampling, o->dirSampling,
        NQ, iDeb
      );
    uint32_t NI = stack->NI;
 
    /* Write image showing selected pixels: */
    multifok_frame_t *fr_sharp = stack->frame[NI-1];
    assert(fr_sharp->zDep == +INF);
    float_image_t *bgrd = fr_sharp->sVal;
    multifok_image_selected_pixels_write(NQ, iDeb, bgrd, o->stackDir);
    
    mfmi_write_pixel_profiles(stack, NQ, iDeb, scene, o->stackDir);
    
    return 0;
  }

multifok_scene_t *mfmi_make_scene
  ( char *sceneType,
    interval_t dom[],
    int32_t NX,
    int32_t NY,
    double zDep,
    bool_t verbose
  )
  {
    fprintf(stderr, "creating scene of type %s", sceneType); 
    /* Define the scene's domain: */
    for (uint32_t j = 0;  j < 3; j++)
      { fprintf(stderr, "%s[%.3f _ %.3f]", (j == 0 ? " size = " : " × "), dom[j].end[0], dom[j].end[1]); }
    fprintf(stderr, "\n");

    /* Parse the {sceneType} defining {floorOnly,flatFloor,minSep}: */
    bool_t floorOnly; /* If true the scene is just the floor, else it has foreground objs. */
    bool_t flatFloor; /* If true the floor is a flat plane at {zMin}, else a ramp. */
    double minSep; /* Min {XY} sep of non-overlapping objs. */
    
    /* The {minSep} is ignored if {floorOnly} is true. */

    if (strcmp(sceneType, "R") == 0)
      { /* Tilted floor surface only: */
        floorOnly = TRUE;
        flatFloor = FALSE;
        minSep = NAN;   /* Irrelevant. */
      }
    else 
      { /* Disk, balls, etc: */
        floorOnly = FALSE;
        flatFloor = TRUE;
        if (strcmp(sceneType, "F") == 0)
          { /* Non-overlapping objects: */
            minSep = 0.5*zDep; 
          }
        else if (strcmp(sceneType, "T") == 0)
          { /* Overlapping objects: */
            minSep = -1.0;
          }
        else
          { demand(FALSE, "invalid scene type"); }
      }

    assert(floorOnly | flatFloor);
    multifok_scene_t *scene = multifok_scene_new(dom, verbose);
    /* Define the object radius range {rMin,rMax} (in scene units): */
    double WX = 2*interval_rad(&(dom[0]));
    double WY = 2*interval_rad(&(dom[1]));
    double WZ = 2*interval_rad(&(dom[2]));
    double Wmin = fmin(WX, WY);
    double scale_out = fmax(NX/WX, NY/WY);
    /* The min obj radius must be at least 7 pixels and 1/100 of scene: */
    double rMin = fmax(7.0/scale_out, 0.01*Wmin); /* Min object radius. */
    /* The max obj radius must be 1/4 of scene X or Y size. */
    /* Don't worry about Z, will be fitted as needed if spherical: */
    double rMax = fmin(0.30*WZ, 0.25*Wmin);       /* Max object radius. */
    multifok_scene_throw_objects
      ( scene, floorOnly, flatFloor, rMin, rMax, minSep, verbose );
    
    assert((scene->objs[0].type == ot_FLAT) || (scene->objs[0].type == ot_RAMP));
    char *floorX = multifok_scene_object_type_to_string(scene->objs[0].type);
    fprintf(stderr, "generated a scene with a %s floor", floorX);
    if (scene->NO == 1)
      { fprintf(stderr, " only.\n"); }
    else
      { char *olapX = (minSep < 0.0 ? "overlapping" : "non-overlapping");
        fprintf(stderr, " and %d %s foreground objects\n", scene->NO-1, olapX);
      }
    return scene;
  }
  
multifok_stack_t *mfmi_make_and_write_stack
  ( char *stackDir,
    int32_t NX,
    int32_t NY, 
    multifok_scene_t *scene,
    char *patternFile,
    r3_t *light_dir,
    double ambient,
    double zFoc_min,
    double zFoc_max,
    double zFoc_step, 
    double zDep,
    uint32_t HS,
    uint32_t KR,
    uint32_t NQ,
    i2_t iDeb[]
  )
  { 
    fprintf(stderr, "creating stack with zDep = %8.4f HS = %d KR = %d\n", zDep, HS, KR); 
    fprintf(stderr, "zFoc varying from %8.4f to %8.4f step %8.4f\n", zFoc_min, zFoc_max, zFoc_step); 
    r3_gen_print(stderr, light_dir, "%+8.5f", "light direction = ( ", " ", " )");
    fprintf(stderr, " ambient = %6.4f\n", ambient); 

    /* Scene size and center in user units:  */
    r3_t size_scene, ctr_scene;
    for (uint32_t j = 0;  j < 3; j++) 
      { interval_t *domj = &(scene->dom[j]);
        size_scene.c[j] = 2*interval_rad(domj);
        ctr_scene.c[j] = interval_mid(domj);
      }
 
    /* Compute the output image pixels per scene unit: */
    double scale_out = fmax(((double)NX)/size_scene.c[0], ((double)NY)/size_scene.c[1]);
    fprintf(stderr, "frame image size = %d × %d", NX, NY);
    fprintf(stderr, " (%.8f pixels per scene unit)\n", scale_out);
  
    /* Read the pattern image: */
    float_image_t *pimg = mfmi_read_pattern_image(patternFile);

    /* Pattern size and center in pixels: */
    int32_t NC_pat, NX_pat, NY_pat;
    float_image_get_size(pimg, &NC_pat, &NX_pat, &NY_pat);
    demand(NC_pat == 1, "pattern image must be grayscale");
    /* We want the scene surface texture on the output image to be close to the pattern texture: */
    double scale_pat = scale_out; /* Pattern image pixels per scene unit. */
    fprintf(stderr, "pattern file = %s size = %d × %d", patternFile, NX_pat, NY_pat); 
    fprintf(stderr, " (%.8f pixels per scene unit)\n", scale_pat); 

    r2_t ctr_pat = (r2_t){{ 0.5*NX_pat, 0.5*NY_pat }};
     
    auto double eval_pattern(double x, double y, double z);
      /* Returns a grayscale pattern value as function of the point with
        scene coordinates {(x,y,z)}. Currently returns the value of {pimg}
        interpolated at {x,y} relative to the center of the scene,
        scaled so that 1 pattern pixel is about 1 image pixel.
        The image {pimg} is replicated to cover the whole plane. */
    
    multifok_stack_t *stack = mfmi_make_and_write_stack_from_pattern_function
      ( stackDir,
        NX, NY, scene, eval_pattern, light_dir, ambient,
        zFoc_min, zFoc_max, zFoc_step,
        zDep, HS, KR, NQ, iDeb
      );
    
    return stack;

    double eval_pattern(double x, double y, double z)
      {
        /* Express {x,y} relative to center of scene: */
        double dx_scene = x - ctr_scene.c[0];
        double dy_scene = y - ctr_scene.c[1];
        
        /* Convert to image units relative to center of pattern: */
        double xr_pat = ctr_pat.c[0] + scale_pat*dx_scene;
        double yr_pat = ctr_pat.c[1] + scale_pat*dy_scene;
        
        /* Eval the pattern image, mirroring as needed: */
        int32_t order = 1;  /* Bicubic C1 interpolation. */
        ix_reduce_mode_t red = ix_reduce_mode_PXMIRR;
        double r = float_image_interpolate_sample(pimg, 0, xr_pat, yr_pat, order, red);
        return r;
      }
  }

multifok_stack_t *mfmi_make_and_write_stack_from_pattern_function
  ( char *stackDir,
    int32_t NX,
    int32_t NY, 
    multifok_scene_t *scene, 
    multifok_scene_raytrace_pattern_t *pattern,
    r3_t *light_dir,
    double ambient,
    double zFoc_min,
    double zFoc_max,
    double zFoc_step,
    double zDep,
    uint32_t HS,
    uint32_t KR,
    uint32_t NQ,
    i2_t iDeb[] 
  ) 
  {
    int32_t debug_level = -1;
    multifok_scene_tree_t *tree = multifok_scene_tree_build(scene->NO, scene->objs, debug_level);
    /* Check that the objects are still OK: */
    multifok_scene_check_object_IDs(scene->NO, scene->objs);

    /* Define the number of frames {NI} (including the sharp one): */
    uint32_t NI = (uint32_t)floor((zFoc_max - zFoc_min)/zFoc_step + 0.0001) + 2;
    
   /* Generate the blurred image stack: */
    multifok_stack_t *stack = multifok_stack_new(NI, 3,NX,NY);
    char *frameDir = NULL;
    for (uint32_t ki = 0;  ki < NI; ki++)
      { double zFoc_fri, zDep_fri;
        uint32_t KR_fri; /* Rays to trace per image sampling point. */
        if (ki == NI-1)
          { /* Sharp frame: */
            zFoc_fri = zFoc_min + 0.5*NI*zFoc_step;
            zDep_fri = +INF;
            KR_fri = 1;
            frameDir = jsprintf("%s/sharp", stackDir); 
          }
        else
          { /* Blurred frame: */
            zFoc_fri = zFoc_min + ki*zFoc_step;
            zDep_fri = zDep;
            KR_fri = KR;
            frameDir = jsprintf("%s/zf%08.4f-df%08.4f", stackDir, zFoc_fri, zDep_fri); 
          }

        /* Create the frame's directory: */
        char *mkdirCmd = jsprintf("mkdir -pv %s", frameDir);
        system(mkdirCmd);
        free(mkdirCmd);

        bool_t verbose = FALSE;
        multifok_frame_t *fri = mfmi_make_and_write_frame
          ( frameDir, NX, NY, scene, tree, pattern, light_dir, ambient,
            zFoc_fri, zDep_fri, HS, KR_fri, verbose,
            NQ, iDeb
          );
        free(frameDir);
          
        stack->frame[ki] = fri;
      }
    return stack;
  }

multifok_frame_t *mfmi_make_and_write_frame
  ( char *frameDir,
    int32_t NX, 
    int32_t NY, 
    multifok_scene_t *scene,
    multifok_scene_tree_t *tree,
    multifok_scene_raytrace_pattern_t *pattern,
    r3_t *light_dir,
    double ambient,
    double zFoc,
    double zDep, 
    uint32_t HS,
    uint32_t KR, 
    bool_t verbose,
    uint32_t NQ,
    i2_t iDeb[]
  )
  {
    double pixSize = multifok_scene_pixel_size(NX, NY, scene->dom); 
    fprintf(stderr, "generating frame images for zFoc = %12.6f  zDep = %12.6f\n", zFoc, zDep);
    fprintf(stderr, "HS = %d KR = %d pixSize = %.6f\n", HS, KR, pixSize);
    r3_gen_print(stderr, light_dir, "%+8.5f", "light direction = ( ", " ", " )");
    fprintf(stderr, " ambient = %6.4f\n", ambient); 

    auto bool_t debug_pixel(i2_t *iPix);
      /* A pixel debugging predicate suitable for 
       the {debug_pixel} argument of {multifok_raytrace_make_frame}.
       It returns {TRUE} if {iPix} is equal to
       one of the pixels selected for debugging. */

    FILE *wr_ray = NULL;
    i2_t iPix_wr = (i2_t){{ -1, -1 }}; /* Pixel for which {wr_ray} was opened. */
        
    auto void report_ray
      ( i2_t *iPix, i2_t *iSmp, double step, double wSmp,
        r3_t *pRay, r3_t *dRay, double wRay, 
        r3_t *pHit, double hHit, double vBlr
      );
      /* Type of a procedure suitable for the {report_ray} argument of
        {multifok_raytrace_make_frame}.
        
        The procedure prints the ray data on {stderr} and also writes it
        to the file {wr_pix}, saving {iPix} to {iPix_wr}. The file
        {wr_pix} is opened if needed, and closed and re-opened whenever
        {iPix} differs from {iPix_wr}. */

    multifok_frame_t *frame = multifok_scene_make_frame
      ( NX, NY, scene, tree, pattern, light_dir, ambient,
        zFoc, zDep, HS, KR, verbose, 
        &debug_pixel, &report_ray, NULL
      );
      
    /* Write the frame out: */
    double hMin = scene->dom[2].end[0];
    double hMax = scene->dom[2].end[1];
    multifok_frame_write(frame, frameDir, hMin, hMax);
    
    return frame;
       
    bool_t debug_pixel(i2_t *iPix)
      { if (! verbose) { return FALSE; }
        bool_t deb = FALSE;
        for (uint32_t kq = 0;  (kq < NQ) && (! deb); kq++)
          { bool_t debug_col = (iPix->c[0] == iDeb[kq].c[0]);
            bool_t debug_row = (iPix->c[1] == iDeb[kq].c[1]);
            deb = deb || (debug_col && debug_row);
          }
        return deb;
      }
    
    void report_ray
      ( i2_t *iPix, i2_t *iSmp, double step, double wSmp,
        r3_t *pRay, r3_t *dRay, double wRay, 
        r3_t *pHit, double hHit, double vBlr
      )
      { if ((wr_ray != NULL) && ((iPix->c[0] != iPix_wr.c[0]) || (iPix->c[1] != iPix_wr.c[1])))
          { fclose(wr_ray);
            iPix_wr = (i2_t){{ -1, -1 }};
            wr_ray = NULL;
          }
        if (wr_ray == NULL)
          { char *fname_ray = jsprintf("%s/pixel-rays-%04d-%04d.txt", frameDir, iPix->c[0], iPix->c[1]);
            wr_ray = open_write(fname_ray, TRUE);
            free(fname_ray);
            iPix_wr = *iPix;
          }
        multifok_raytrace_show_ray_data(stderr, 8, iPix, pixSize, iSmp, step, wSmp, pRay, dRay, wRay, pHit, hHit, vBlr);
        multifok_raytrace_write_ray_data(wr_ray, iPix, pixSize, iSmp, step, wSmp, pRay, dRay, wRay, pHit, hHit, vBlr);
      }
  }

#define mfmi_image_gamma 1.000
  /* Assumed encoding gamma of input and output images. */

#define mfmi_image_bias 0.0327
  /* Assumed encoding bisa of input and output images. */

float_image_t *mfmi_read_pattern_image(char *fname)
  {
    image_file_format_t ffmt = image_file_format_PNG;
    double gammaDec[1];
    double bias[1];
    bool_t verbose = TRUE;
    float_image_t *pat = float_image_read_gen_named(fname, ffmt, 0.0, 1.0, NULL, gammaDec, bias, verbose);
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

void mfmi_select_debug_pixels
  ( int32_t NX,
    int32_t NY,
    multifok_scene_t *scene,
    uint32_t NQ_max,
    uint32_t *NQ_P,
    i2_t **iDeb_P
  )
  {
    bool_t verbose = TRUE;
    
    fprintf(stderr, "selecting plot pixels...\n");
    uint32_t NO = scene->NO; /* Number of objects inscene. */
    
    /* Select the pixels: */
    i2_vec_t iDeb = i2_vec_new(NQ_max);
    uint32_t NQ = 0; /* Pixels actually selected. */
    
    auto void selpix(double x_scene, double y_scene);
      /* Converts the point {(x_scene,y_scene)} from scene coordinates to
        a pair of pixel indices, and appends that pair to the list of sample pixels. */
    
    if (NO >= 2)
      { uint32_t NQ_exp = (NQ_max < NO-1 ? NQ_max : NO-1);
        /* There is at least one foreground object. */
        /* Pick {NQ_exp} object centers, preferably disks. */
        /* Count the number of disks {NO_disk}: */
        uint32_t NO_disk = 0;
        for (uint32_t ko = 0;  ko < NO; ko++)
          { multifok_scene_object_t *objk = &(scene->objs[ko]);
            if (objk->type == ot_DISK) { NO_disk++; }
          }
        bool_t disk_only = (NO_disk >= NQ_exp); /* TRUE if we can choose only disks. */

        /* The number of candidate objects remaining is {NO_rem}: */
        uint32_t NO_rem = (disk_only ? NO_disk : NO-1); 
        for (uint32_t ko = 0;  ko < NO; ko++)
          { if ((NQ >= NQ_exp) || (NO_rem == 0)) { break; }
            multifok_scene_object_t *objk = &(scene->objs[ko]);
            if ((objk->type == ot_FLAT) || (objk->type == ot_RAMP)) { continue; }
            if ((! disk_only) || (objk->type == ot_DISK))
              { uint32_t NQ_rem = NQ_exp - NQ;  /* Number of pixels still to be selected. */
                double prpick = (NO_rem <= NQ_rem ? 1.0 : ((double)NQ_rem)/NO_rem);
                if (drandom() < prpick) 
                  { r3_t ctr = multifok_scene_box_center(objk->bbox);
                    selpix(ctr.c[0], ctr.c[1]);
                  }
                NO_rem--;
              }
          }
      }
    else
      { /* Scene has only the floor: */
        assert(NO == 1);
        multifok_scene_object_t *obj0 = &(scene->objs[0]);
        assert((obj0->type == ot_FLAT) || (obj0->type == ot_RAMP));
        uint32_t NQ_exp = (obj0->type == ot_FLAT ? 1 : 5); /* Redefine the goal... */
        interval_t *xr = &(scene->dom[0]);
        interval_t *yr = &(scene->dom[1]);
        for (uint32_t kq = 0;  kq < NQ_exp; kq++)
          { double fr = (NQ_exp == 0 ? 0.0 : 0.95*((2.0*kq)/(NQ_exp-1) - 1.0));
            double x_scene = interval_mid(xr) + fr*interval_rad(xr);
            double y_scene = interval_mid(yr) + fr*interval_rad(yr);
            selpix(x_scene, y_scene);
          }
      }
    
    i2_vec_trim(&iDeb, NQ);
    
    (*NQ_P) = NQ;
    (*iDeb_P) = iDeb.e;
      
    void selpix(double x_scene, double y_scene) 
      { 
        r3_t p_scene = (r3_t){{ x_scene, y_scene, 0.0 }};  /* The {Z} coordinate is irrelevant. */
        r3_t p_img = multifok_scene_coords_to_image_coords(&p_scene, 0.0, scene->dom, NX, NY);
        int32_t ix = (int32_t)floor(p_img.c[0] + 0.5);
        int32_t iy = (int32_t)floor(p_img.c[1] + 0.5);
        if ((ix >= 0) && (ix < NX) && (iy >= 0) && (iy < NY))
          { if (verbose) { fprintf(stderr, "  selected %d,%d\n", ix, iy); }
            i2_vec_expand(&iDeb, NQ);
            iDeb.e[NQ] = (i2_t){{ ix, iy }};
            NQ++;
          }
      }
  }    

void mfmi_write_pixel_profiles
  ( multifok_stack_t *stack,
    uint32_t NQ,
    i2_t iDeb[],
    multifok_scene_t *scene,
    char *stackDir
  )
  { 
    int32_t NC = stack->NC;
    int32_t NX = stack->NX;
    int32_t NY = stack->NY;
    
    double pixSize = multifok_scene_pixel_size(NX, NY, scene->dom); 

    /* Get the sharp frame's scene view for background: */
    uint32_t NI = stack->NI;
    
    for (uint32_t kq = 0;  kq < NQ; kq++)
      { i2_t *iPix = &(iDeb[kq]);
        int32_t ix = iPix->c[0];
        int32_t iy = iPix->c[1];
        
        if ((ix >= 0) && (ix < NX) && (iy >= 0) && (iy < NY))
          { 
            /* Open the file for this pixel: */
            char *fname = jsprintf("%s/pixel-data-%04d-%04d.txt", stackDir, ix, iy);
            FILE *wr = open_write(fname, TRUE);

            for (uint32_t ki = 0;  ki < NI; ki++)
              { multifok_frame_t *fri = stack->frame[ki];
                if (isfinite(fri->zDep))
                  { float hAvg = float_image_get_sample(fri->hAvg, 0, ix, iy);
                    float hDev = float_image_get_sample(fri->hDev, 0, ix, iy);
                    float shrp = float_image_get_sample(fri->shrp, 0, ix, iy);
                    float colr[3];
                    for (int32_t ic = 0;  ic < 3; ic++)
                      { colr[ic] = (float)(ic < NC ? float_image_get_sample(fri->sVal, ic, ix, iy) : 0.0); }
                    r3_t ctr = (r3_t){{ ix+0.5, iy+0.5, 0.0 }};
                    r3_t pCtr = multifok_scene_coords_from_image_coords(&ctr, NX, NY, scene->dom, fri->zFoc);
                    multifok_raytrace_write_pixel_data
                      ( wr, iPix, pixSize, &pCtr, fri->zFoc, fri->zDep,
                        shrp, hAvg, hDev, NC, colr
                      );
                    fprintf(wr, "\n");
                  }
              }
            fclose(wr);
            free(fname);
          }
      }
  }

mfmi_options_t *mfmi_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    mfmi_options_t *o = notnull(malloc(sizeof(mfmi_options_t)), "no mem");

    argparser_get_keyword(pp, "-imageSize");
    o->imageSize_X = (int32_t)argparser_get_next_int(pp, 30, 4096);
    o->imageSize_Y = (int32_t)argparser_get_next_int(pp, 30, 4096);

    argparser_get_keyword(pp, "-sceneSize");
    for (uint32_t j = 0;  j < 3; j++)
      { double lo = argparser_get_next_double(pp, -4096.0, +4096.0);
        double hi = argparser_get_next_double(pp, lo + 1.0, +4096.0);
        o->sceneSize[j] = (interval_t){{ lo, hi }};
      }

    argparser_get_keyword(pp, "-sceneType");
    o->sceneType = argparser_get_next_non_keyword(pp);  
    
    argparser_get_keyword(pp, "-pixSampling");
    o->pixSampling = (uint32_t)argparser_get_next_int(pp, 0, 9999);  

    argparser_get_keyword(pp, "-dirSampling");
    o->dirSampling = (uint32_t)argparser_get_next_int(pp, 1, 9999);  

    argparser_get_keyword(pp, "-dephOfFocus");
    o->dephOfFocus = argparser_get_next_double(pp, 0.1, 1.0e200);  

    argparser_get_keyword(pp, "-focusHeight");
    o->focusHeight_lo = argparser_get_next_double(pp, -4096.0, +4096.0);  
    argparser_get_keyword_next(pp, "to");
    o->focusHeight_hi = argparser_get_next_double(pp, o->focusHeight_lo, +4096.0);  
    argparser_get_keyword_next(pp, "step");
    o->focusHeight_step = argparser_get_next_double(pp, 0.001, 4096.0);  

    argparser_get_keyword(pp, "-patternFile");
    o->patternFile = argparser_get_next(pp);
    
    argparser_get_keyword(pp, "-lightDir");
    r3_t ldir = argparser_get_next_r3_dir(pp);
    if (isnan(ldir.c[0]) || isnan(ldir.c[1]) || isnan(ldir.c[2]))
      { o->lightDir = NULL; }
    else
      { o->lightDir = talloc(1, r3_t);
        (*(o->lightDir)) = ldir;
      }

    argparser_get_keyword(pp, "-ambient");
    o->ambient = argparser_get_next_double(pp, 0.0, 1.0);  

    argparser_get_keyword(pp, "-stackDir");
    o->stackDir = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
