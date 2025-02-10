#define PROG_NAME "test_mfok_scene_make_frame"
#define PROG_DESC "test of {multifok_test_image_make.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-02-10 09:06:28 by stolfi */ 
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
#include <interval.h>
#include <box.h>
#include <ix_reduce.h>
#include <r3.h>
#include <i2.h>
#include <frgb.h>
#include <frgb_ops.h>
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
    int32_t imageSize_Y;       /* Image height (not including the grayscale chart). */
    char *sceneType;           /* Type of scene. See below. */
    uint32_t pixSampling;      /* Number of sampling points per axis and pixel. */
    uint32_t dirSampling;      /* Min number of aperture rays per point. */
    char *patternFile;         /* File name of pattern to use to paint objects. */
    r3_t *lightDir;            /* Light source direction, or {NULL}. */
    double ambient;            /* Intensity of "ambient" (isotropic) light. */
    double gloss;              /* Intensity of "glossy" surface finish. */
    char *stackFolder;         /* Parent folder for frames of the stack. */
    /* These parameters are specified in user units: */
    interval_t sceneSize[3];   /* Nominal coordinate ranges of scene. */
    double dephOfFocus;        /* Depth of focus. */
    double focusHeight_lo;     /* Lowest focus plane {Z}. */
    double focusHeight_hi;     /* Highest focus plane {Z}. */
    double focusHeight_step;   /* {Z} increment between focus planes. */
  } mfmi_options_t;
  /* Command line parameters. 
    The {sceneType} may be:
     "R" Ramp only.
     "D" Single disk.
     "B" Single ball.
     "C" Single cone.
     "P" Single square pyramid.
     "F" Non-overlapping mixed objects.
     "T" Overlapping mixed objects.
  */
  
#define PROG_OUTPUT_INFO \
  "  RUN DIRECTORY\n" \
  "    The program writes all its output files to a specified" \
  " folder \"{stackFolder}\".  This folder has one sub-folder for each frame of" \
  " the multifocus stack," \
  " called \"{stackFolder}/{frameFolder}\"." \
  "  " multifok_FRAME_DIR_INFO "\n" \
  "\n" \
  "  \n" \
  "  FRAME IMAGE FILES\n" \
  "\n" \
  "    " multifok_FRAME_FILES_INFO "\n" \
  "\n" \
  "  PIXEL V FOCUS PLOT FILES\n" \
  "\n" \
  "    The program also writes to the \"{stackFolder}\" zero or more" \
  " files \"pixel-data-{XXXX}-{YYYY}.txt\" where {XXXX} and {YYYY} are" \
  " the column and" \
  " row indices {ix,iy} of a pixel that was selected for detailed" \
  " debugging.  Each of" \
  " these files contain one line for each frame, excluding the" \
  " sharp one.  " multifok_raytrace_write_pixel_data_INFO "\n" \
  "    These files can be used" \
  " to plot the variation of pixel properties as a function of the" \
  " parameters {zFoc} and/or {zDep}.\n" \
  "\n" \
  "  RAY DATA FILES\n" \
  "\n" \
  "    The program will also write in each frame directory \"{stackFolder}/{frameFolder}\" zero" \
  " or more files \"pixel-rays-{XXXX}-{YYYY}.txt\" where {XXXX} and {YYYY} are the column" \
  " and row indices {ix,iy} of a selected debugging pixel, as above.  This file contains one" \
  " line for each ray that was cast in order to compute the values of the frame images" \
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
#define ot_PYRA multifok_scene_object_type_PYRA

typedef enum { 
    mfmi_scene_kind_EMPTY,
    mfmi_scene_kind_SINGU,
    mfmi_scene_kind_MULTI,
    mfmi_scene_kind_QUART
  } mfmi_scene_kind_t;
  
#define sk_EMPTY mfmi_scene_kind_EMPTY
#define sk_SINGU mfmi_scene_kind_SINGU
#define sk_MULTI mfmi_scene_kind_MULTI
#define sk_QUART mfmi_scene_kind_QUART

int32_t main(int32_t argn, char **argv);

mfmi_options_t *mfmi_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */
  
multifok_scene_t *mfmi_make_scene
  ( interval_t dom[],
    char *sceneType,
    double gloss,
    double pixelSize,
    double zDep,
    bool_t verbose
  );
  /* Generates a test scene with some number of objects, roughly
    spanning the box {dom[0]×dom[1]×dom[2]} in scene coordinates.
    
    If {sceneType} is "R", the scene have only a floor object
    of type {ot_RAMP}. See {mfmi_fill_scene_empty}.
    
    If {sceneType} is "D", "B", "C", or "P", the scene will have a flat floor (an
    object of type {ot_FLAT}) at {Z=0.0} plus a single foreground object
    of type {ot_DISK}, {ot_BALL}, {ot_CONE}, or {ot_PYRA}, respectively. 
    See {mfmi_fill_scene_singu}.
    
    If {sceneType} is "F" or "T", the scene will have a flat floor at
    {Z=0.0}, and there will be some number of foreground objects at
    random {XY} coordinates. The parameters {pixelSize} and {zDep} are
    relevant only for this kind of scene. See {mfmi_fill_scene_multi}.
    The {overlap} flag will be true if {sceneType} is "T", and false if
    {sceneType} is "F". 
    
    If {sceneType} is "Q", the scene will be like "F", but there will be
    only four foreground objects, one of each type {DISK,BALL,CONE,PYRA}.
    See {mfmi_fill_scene_quartet}*/

void mfmi_fill_scene_empty
  ( multifok_scene_t *scene,
    multifok_scene_object_type_t type,
    bool_t verbose
  );
  /* The {scene} must have no objects.  Adds a floor object of the
    specified type ({ot_FLAT} or {ot_RAMP}). */

void mfmi_fill_scene_singu
  ( multifok_scene_t *scene,
    multifok_scene_object_type_t type,
    frgb_t *fgGlo,
    bool_t verbose
  );
  /* The {scene} must have no objects. Adds a floor object of the type
    {ot_FLAT} at {Z=0} and a single foreground object of the specified type
    (which must not be {ot_FLAT} or {ot_RAMP}). The object will be
    centered in the scene's domain {scene.dom}. */
  
void mfmi_fill_scene_multi
  ( multifok_scene_t *scene,
    frgb_t *fgGlo,
    double pixelSize,
    bool_t overlap,
    double zDep,
    bool_t verbose
  );
  /* The {scene} must have no objects.  Adds a floor object of the
    type {ot_FLAT} and zero or more foreground objects of random
    types, sizes, and positions.
    
    The {Z} spans of all these foreground objects will be strictly
    contained in the {Z} range {scene.dom[2]}. 
    
    If the {overlap} flag is true, the {XY} projections of the
    foreground objects will be disjoint and contained in the {XY}
    projection of {scene.dom}, even when out of focus (assuming that the
    depth of focus is {zDep}). Otherwise the {XY} projections of the
    foreground objects will probably overlap and extend outside the {XY}
    projection of {scene.dom}; only their centers will be in that
    rectangle.
    
    The {pixelSize} parameter is the size of one pixel in SCUs. 
    It affects the min size of objects. */
  
void mfmi_fill_scene_quartet
  ( multifok_scene_t *scene,
    frgb_t *fgGlo,
    double pixelSize,
    double zDep,
    bool_t verbose
  );
  /* The {scene} must have no objects.  Adds a floor object of the
    type {ot_FLAT} and four foreground objects of types
    {DISK,BALL,CONE,PYRA}, as large as possible.
    
    The {Z} range of each foreground object will be limited to 
    the {Z} range of the scene ({scene.dom[2]}).  Their {XY} projections
    will be disjoint and contained in the {XY}
    projection of {scene.dom}, even when out of focus (assuming that the
    depth of focus is {zDep}).
    
    The {pixelSize} parameter is the size of one pixel in SCUs. 
    It affects the min size of objects. */

multifok_stack_t *mfmi_make_and_write_stack
  ( char *stackFolder,
    int32_t NX,
    int32_t NY, 
    multifok_scene_t *scene,
    char *patternFile,
    r3_t *uLit,
    frgb_t *cLit,
    frgb_t *cIso,
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
    "{stackFolder}/{frameFolder}/pixel-rays-{XXXX}-{YYYY}.txt", where
    {XXXX,YYYY} are the pixel indices formatted as "%04d". Here
    {stackFolder} is as explained in {PROG_INFO}. */

multifok_stack_t *mfmi_make_and_write_stack_from_pattern_function
  ( char *stackFolder,
    int32_t NX,
    int32_t NY,
    multifok_scene_t *scene, 
    multifok_pattern_double_proc_t *patternGlo,
    multifok_pattern_double_proc_t *patternLam,
    r3_t *uLit,
    frgb_t *cLit,
    frgb_t *cIso,
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
  ( char *frameFolder,
    int32_t NX, 
    int32_t NY, 
    multifok_scene_t *scene,
    multifok_scene_tree_t *tree,
    multifok_pattern_double_proc_t *patternGlo,
    multifok_pattern_double_proc_t *patternLam,
    r3_t *uLit,
    frgb_t *cLit,
    frgb_t *cIso,
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
    
    Writes the images to "{frameFolder}/{tag}.png" were {tag} is "sVal",
    "hAvg", "hDev", "sNrm", and "shrp". Also writes pixel ray data files
    "{frameFolder}/pixel-rays-{XXXX}-{YYYY}.txt" */

void mfmi_write_pixel_profiles
  ( multifok_stack_t *stack,
    uint32_t NQ,
    i2_t iDeb[],
    multifok_scene_t *scene,
    char *stackFolder
  );
  /* Given a list of {NQ} of pixels {iDeb[0..NQ-1]} in the 
    image domain, writes a file "{stackFolder}/pixel-data-{XXXX}-{YYYY}.txt"
    for each pixel, where {XXXX} an {YYYY} are the pixel column and row indices.
    
    Each file will have up to {NI = stack.NI} lines with format

      "{ki} {zFoc} {zDep}  {hAvg} {hDev} {shrp} {sNrm.x} {sNrm.y} {sNrm.z} {sVal[0]} ...  {sVal[NC-1]}"

    where {ki} is a frame number in {0..NI-1} and {hAvg}, and 
    {hDev}, {shrp}, {sNrm}, and {sVal[0..NC-1]} are the values of
    that pixel in the images {hAvg} {hDev}, {shrp}, and {sNrm}, and 
    {sVal} of {stack.frame[ki]}. Note that {zDep} may be {+INF}.  */
    
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
    
FILE *mfmi_open_pixel_plot_data_file(char *stackFolder);
  /* Opens the file "{stackFolder}/pixplot.txt" for writing. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfmi_options_t *o = mfmi_parse_options(argc, argv);
    
    int32_t NX = o->imageSize_X;
    int32_t NY = o->imageSize_Y;
    
    interval_t *dom = o->sceneSize;
    double wX = 2*interval_rad(&(dom[0]));
    double wY = 2*interval_rad(&(dom[1]));
    double pixelSize = fmin(wX/NX, wY/NY);

    /* Create the test scene: */
    bool_t verbose_scene = TRUE;
    multifok_scene_t *scene = mfmi_make_scene
      ( o->sceneSize, o->sceneType, o->gloss,
        pixelSize, o->dephOfFocus, 
        verbose_scene
      );
      
    /* Select the debugging pixels: */
    uint32_t NQ_max = 10;
    uint32_t NQ;
    i2_t *iDeb;
    /* !!! {mfmi_select_debug_pixels} calls {drandom} !!! */
    /* !!! So call it as before but then reset {NQ} !!! */
    mfmi_select_debug_pixels(NX, NY, scene, NQ_max, &NQ, &iDeb);
    NQ = 0;
    
    frgb_t lightColor = (frgb_t){{ 0.900f, 0.850f, 0.800f }}; /* Color of point source. */
    frgb_t ambColor = (frgb_t){{ 0.800f, 0.900f, 1.000f }}; 
    ambColor = frgb_scale(o->ambient, &ambColor);

    multifok_stack_t *stack = mfmi_make_and_write_stack
      ( o->stackFolder,
        NX, NY, scene, o->patternFile, o->lightDir, &(lightColor), &(ambColor),
        o->focusHeight_lo, o->focusHeight_hi, o-> focusHeight_step,
        o->dephOfFocus, o->pixSampling, o->dirSampling,
        NQ, iDeb
      );
    uint32_t NI = stack->NI;
 
    /* Write image showing selected pixels: */
    multifok_frame_t *fr_sharp = stack->frame[NI-1];
    assert(fr_sharp->zDep == +INF);
    float_image_t *bgrd = fr_sharp->sVal;
    multifok_image_selected_pixels_write(NQ, iDeb, bgrd, o->stackFolder);
    
    mfmi_write_pixel_profiles(stack, NQ, iDeb, scene, o->stackFolder);
    
    return 0;
  }

multifok_scene_t *mfmi_make_scene
  ( interval_t dom[],
    char *sceneType,
    double gloss,
    double pixelSize,
    double zDep,
    bool_t verbose
  )
  {
    if (verbose)
      { fprintf(stderr, "creating scene of kind '%s'", sceneType); 
        for (uint32_t j = 0;  j < 3; j++)
          { char *ps = (j == 0 ? " size = " : " × ");
            fprintf(stderr, "%s[%.3f _ %.3f]", ps, dom[j].end[0], dom[j].end[1]);
          }
        fprintf(stderr, "\n");
      }
    demand(strlen(sceneType) == 1, "the scene type should be a single char");
 
    multifok_scene_t *scene = multifok_scene_new(dom, verbose);
    
    frgb_t fgGlo = (frgb_t){{ 0.300f, 1.000f, 0.500f }}; 
    fgGlo = frgb_scale(gloss, &(fgGlo));

    multifok_scene_object_type_t type = ot_FLAT;
    mfmi_scene_kind_t kind = sk_EMPTY;
    bool_t overlap = FALSE;
    switch (sceneType[0])
      { case 'R': kind = sk_EMPTY; type = ot_RAMP;  break;
        case 'F': kind = sk_MULTI; overlap = FALSE; break;
        case 'T': kind = sk_MULTI; overlap = TRUE;  break;
        
        case 'Q': kind = sk_QUART; break;

        case 'D': kind = sk_SINGU; type = ot_DISK; break;
        case 'B': kind = sk_SINGU; type = ot_BALL; break;
        case 'C': kind = sk_SINGU; type = ot_CONE; break;
        case 'P': kind = sk_SINGU; type = ot_PYRA; break;
        default: demand(FALSE, "invalid scene type");
      }
        
    switch (kind)
      { case sk_EMPTY:
          mfmi_fill_scene_empty(scene, type, verbose); break;
         
        case sk_SINGU:
          mfmi_fill_scene_singu(scene, type, &(fgGlo), verbose); break;
          
        case sk_MULTI:
          mfmi_fill_scene_multi(scene, &(fgGlo), pixelSize, overlap, zDep, verbose); break;
          
        case sk_QUART:
          mfmi_fill_scene_quartet(scene, &(fgGlo), pixelSize, zDep, verbose); break;
        
        default: assert(FALSE);
      }

    assert((scene->objs[0].type == ot_FLAT) || (scene->objs[0].type == ot_RAMP));
    if (verbose) { multifok_scene_print(stderr, scene); }
    return scene;
  }
      
void mfmi_fill_scene_empty
  ( multifok_scene_t *scene,
    multifok_scene_object_type_t type,
    bool_t verbose
  )
  {
    if (verbose) 
      { char *typeX = multifok_scene_object_type_to_string(type);
        fprintf(stderr, "generating a scene with only floor of type %s\n", typeX);
      }
    multifok_scene_add_floor(scene, type, verbose);
  }
  
void mfmi_fill_scene_singu
  ( multifok_scene_t *scene,
    multifok_scene_object_type_t type,
    frgb_t *fgGlo,
    bool_t verbose
  )  
  { 
    if (verbose) 
      { char *typeX = multifok_scene_object_type_to_string(type);
        fprintf(stderr, "filling scene with floor and one foreground object of type %s\n", typeX); 
      }
    
    /* In principle the object can use the whole {scene.dom} minus a margin: */
    interval_t *dom = scene->dom;
    double wX = interval_width(&(dom[0]));
    double wY = interval_width(&(dom[1]));
    
    double wSceneXY = fmin(wX, wY);
    double wSceneZ = interval_width(&(dom[2]));
    double marginXY = 0.05*wSceneXY; 
    double marginZ = 0.0005*wSceneZ;
    double slackX = wX - wSceneXY;
    double slackY = wY - wSceneXY;

    double wObjXY = wSceneXY - marginXY;
    if (verbose) { fprintf(stderr, "object {XY} size = %12.6f\n", wObjXY); }
    
    interval_t bbox[3]; 

    LO(bbox[0]) = 0.5*slackX + marginXY;
    HI(bbox[0]) = LO(bbox[0]) + wObjXY;

    LO(bbox[1]) = 0.5*slackY + marginXY;
    HI(bbox[1]) = LO(bbox[1]) + wObjXY;

    bbox[2] = dom[2];
    interval_widen(&(bbox[2]), -marginZ);
    
    /* Pick colors: */
    frgb_t fgLam = (frgb_t){{ 1.000f, 0.700f, 0.300f }};
    frgb_t bgLam = (frgb_t){{ 0.600f, 0.500f, 0.250f }};
      
    frgb_t bgGlo = (frgb_t){{ 0.000f, 0.000f, 0.000f }};
      
    if (verbose) { fprintf(stderr, "adding FLAT floor ...\n"); }
    multifok_scene_add_floor(scene, ot_FLAT, verbose);

    if (verbose) { fprintf(stderr, "adding the foreground object ...\n"); }
    multifok_scene_add_foreground_object(scene, type, bbox, fgGlo, &bgGlo, &fgLam, &bgLam, verbose);
  }
  
void mfmi_fill_scene_multi
  ( multifok_scene_t *scene,
    frgb_t *fgGlo,
    double pixelSize,
    bool_t overlap,
    double zDep,
    bool_t verbose
  ) 
  { 
    if (verbose) 
      { char *olapX = (overlap ? "overlapping" : "non-overlapping");
        fprintf(stderr, "filling scene with floor and multiple random %s foreground objects\n", olapX);
      }

    double minSep = (overlap ? -1.0 : 0.5*zDep);
    
    /* Define the object {XY} width range {wMinXY,wMaxXY} (in scene units): */
    interval_t *dom = scene->dom;
    double wX = interval_width(&(dom[0]));
    double wY = interval_width(&(dom[1]));
    double wXYMin = fmin(wX, wY);

    /* Define the range {wMinXY,wMaxXY} of object XY widths: */
    /* Don't worry about Z, will be fitted as needed: */
    double wMinXY = fmax(30.0*pixelSize, 0.05*wXYMin);
    double wMaxXY = 0.33*wXYMin;
    if (verbose) { fprintf(stderr, "object {XY} size range = [ %12.6f _ %12.6f ]\n", wMinXY, wMaxXY); }

    if (verbose) { fprintf(stderr, "adding FLAT floor ...\n"); }
    multifok_scene_add_floor(scene, ot_FLAT, verbose);

    if (verbose) { fprintf(stderr, "throwing foreground objects ...\n"); }
    multifok_scene_throw_foreground_objects(scene, wMinXY, wMaxXY, fgGlo, minSep, verbose);
  }
  
void mfmi_fill_scene_quartet
  ( multifok_scene_t *scene,
    frgb_t *fgGlo,
    double pixelSize,
    double zDep,
    bool_t verbose
  ) 
  { 
    if (verbose) 
      { fprintf(stderr, "filling scene with floor and four foreground objects to scene\n"); }

    double minSep = 0.5*zDep;

    if (verbose) { fprintf(stderr, "adding FLAT floor ...\n"); }
    multifok_scene_add_floor(scene, ot_FLAT, verbose);
    
    interval_t *dom = scene->dom;
    
    /* Define the foreground object's {XY} width {wXY} (in scene units): */
    double wX = interval_width(&(dom[0]));
    double wY = interval_width(&(dom[1]));
    double wZ = interval_width(&(dom[2]));
    double wSceneXY = fmin(wX, wY);
    double slackX = wX - wSceneXY;
    double slackY = wY - wSceneXY;

    double marginXY = 0.05*wSceneXY + minSep; 
    double marginZ = 0.0005*wZ; 

    /* Don't worry about Z, will be fitted as needed: */
    double wObjXY = (wSceneXY - 3*marginXY)/2;
    if (verbose) { fprintf(stderr, "object {XY} size = %12.6f\n", wObjXY); }
    
    auto void place_object(int32_t kx, int32_t ky, multifok_scene_object_type_t type);
    
    place_object(0,0, ot_DISK);
    place_object(1,0, ot_BALL);
    place_object(0,1, ot_PYRA);
    place_object(1,1, ot_CONE);
    
    return;
    
    void place_object(int32_t kx, int32_t ky, multifok_scene_object_type_t type)
      {
        interval_t bbox[3];
        
        LO(bbox[0]) = 0.5*slackX + marginXY + kx*(wObjXY + marginXY);
        HI(bbox[0]) = LO(bbox[0]) + wObjXY;
        
        LO(bbox[1]) = 0.5*slackY + marginXY + ky*(wObjXY + marginXY);
        HI(bbox[1]) = LO(bbox[1]) + wObjXY;
        
        bbox[2] = dom[2];
        interval_widen(&(bbox[2]), -marginZ);

        /* Pick colors: */
        frgb_t fgLam = (frgb_t){{ (float)drandom(), 1.000, 0.600f }}; frgb_from_HTY(&(fgLam));
        frgb_t bgLam = (frgb_t){{ (float)drandom(), 1.000, 0.200f }}; frgb_from_HTY(&(bgLam));

        frgb_t bgGlo = (frgb_t){{ 0.000f, 0.000f, 0.000f }};

        if (verbose) { fprintf(stderr, "adding the foreground object ...\n"); }
        multifok_scene_add_foreground_object(scene, type, bbox, fgGlo, &bgGlo, &fgLam, &bgLam, verbose);
      }
  }
  
multifok_stack_t *mfmi_make_and_write_stack
  ( char *stackFolder,
    int32_t NX,
    int32_t NY, 
    multifok_scene_t *scene,
    char *patternFile,
    r3_t *uLit,
    frgb_t *cLit,
    frgb_t *cIso,
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
    r3_gen_print(stderr, uLit, "%+8.5f", "light direction = ( ", " ", " )");
    frgb_print(stderr, " color = ( ", cLit, 3, "%7.5f", " )\n"); 
    frgb_print(stderr, " ambient = ( ", cIso, 3, "%7.5f", " )\n"); 

    /* Scene size and center in user units:  */
    r3_t size_scene;
    for (uint32_t j = 0;  j < 3; j++) 
      { interval_t *domj = &(scene->dom[j]);
        size_scene.c[j] = 2*interval_rad(domj);
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
     
    auto double eval_pattern(r3_t *q, float_image_t *img);
      /* Return the value of {img} interpolated at {q.x,q.y}, scaled so
        that 1 pattern pixel is about 1 image pixel, and relative to the
        center of the pattern. The image {img} is replicated to cover
        the whole plane. */

    auto double patternGlo(r3_t *q);
    auto double patternLam(r3_t *q);
      /* Each of these procedures returns a grayscale pattern value as
        function of the point with scene coordinates {q}.  Currently both 
        call {eval_pattern(q,pimg)} with the same image {pimg}. */
    
    multifok_stack_t *stack = mfmi_make_and_write_stack_from_pattern_function
      ( stackFolder,
        NX, NY, scene, patternGlo, patternLam, uLit, cLit, cIso,
        zFoc_min, zFoc_max, zFoc_step,
        zDep, HS, KR, NQ, iDeb
      );
    
    return stack;

    double patternGlo(r3_t *q)
      { return eval_pattern(q, pimg); }
      
    double patternLam(r3_t *q)
      { return eval_pattern(q, pimg); }
      
    double eval_pattern(r3_t *q, float_image_t *img)
      {
        double dx = q->c[0];
        double dy = q->c[1];
        
        /* Convert to image units relative to center of pattern: */
        double xr_pat = ctr_pat.c[0] + scale_pat*dx;
        double yr_pat = ctr_pat.c[1] + scale_pat*dy;
        
        /* Eval the pattern image, mirroring as needed: */
        int32_t order = 1;  /* Bicubic C1 interpolation. */
        ix_reduce_mode_t red = ix_reduce_mode_PXMIRR;
        double r = float_image_interpolate_sample(img, 0, xr_pat, yr_pat, order, red);
        return r;
      }
  }

multifok_stack_t *mfmi_make_and_write_stack_from_pattern_function
  ( char *stackFolder,
    int32_t NX,
    int32_t NY, 
    multifok_scene_t *scene, 
    multifok_pattern_double_proc_t *patternGlo,
    multifok_pattern_double_proc_t *patternLam,
    r3_t *uLit,
    frgb_t *cLit,
    frgb_t *cIso,
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
    char *frameFolder = NULL;
    for (uint32_t ki = 0;  ki < NI; ki++)
      { double zFoc_fri, zDep_fri;
        uint32_t KR_fri; /* Rays to trace per image sampling point. */
        if (ki == NI-1)
          { /* Sharp frame: */
            zFoc_fri = zFoc_min + 0.5*NI*zFoc_step;
            zDep_fri = +INF;
            KR_fri = 1;
            frameFolder = jsprintf("%s/sharp", stackFolder); 
          }
        else
          { /* Blurred frame: */
            zFoc_fri = zFoc_min + ki*zFoc_step;
            zDep_fri = zDep;
            KR_fri = KR;
            frameFolder = jsprintf("%s/zf%08.4f-df%08.4f", stackFolder, zFoc_fri, zDep_fri); 
          }

        /* Create the frame's directory: */
        char *mkdirCmd = jsprintf("mkdir -pv %s", frameFolder);
        system(mkdirCmd);
        free(mkdirCmd);

        bool_t verbose = FALSE;
        multifok_frame_t *fri = mfmi_make_and_write_frame
          ( frameFolder, NX, NY, scene, tree, patternGlo, patternLam, uLit, cLit, cIso,
            zFoc_fri, zDep_fri, HS, KR_fri, verbose,
            NQ, iDeb
          );
        free(frameFolder);
          
        stack->frame[ki] = fri;
      }
    return stack;
  }

multifok_frame_t *mfmi_make_and_write_frame
  ( char *frameFolder,
    int32_t NX, 
    int32_t NY, 
    multifok_scene_t *scene,
    multifok_scene_tree_t *tree,
    multifok_pattern_double_proc_t *patternGlo,
    multifok_pattern_double_proc_t *patternLam,
    r3_t *uLit,
    frgb_t *cLit,
    frgb_t *cIso,
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
    r3_gen_print(stderr, uLit, "%+8.5f", "light direction = ( ", " ", " )");
    frgb_print(stderr, " color = ( ", cLit, 3, "%7.5f", " )\n"); 
    frgb_print(stderr, " ambient = ( ", cIso, 3, "%7.5f", " )\n"); 

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
        r3_t *pHit, double hHit, double vBlr, 
        r3_t *sNrm, int32_t NC, float sVal[]
      );
      /* Type of a procedure suitable for the {report_ray} argument of
        {multifok_raytrace_make_frame}.
        
        The procedure prints the ray data on {stderr} and also writes it
        to the file {wr_pix}, saving {iPix} to {iPix_wr}. The file
        {wr_pix} is opened if needed, and closed and re-opened whenever
        {iPix} differs from {iPix_wr}. */

    multifok_frame_t *frame = multifok_scene_make_frame
      ( NX, NY, scene, tree, patternGlo, patternLam, uLit, cLit, cIso,
        zFoc, zDep, HS, KR, verbose, 
        &debug_pixel, &report_ray, NULL
      );
      
    /* Write the frame out: */
    double hMin = scene->dom[2].end[0];
    double hMax = scene->dom[2].end[1];
    multifok_frame_write(frame, frameFolder, hMin, hMax);
    
    return frame;
       
    bool_t debug_pixel(i2_t *iPix)
      { bool_t deb = FALSE;
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
        r3_t *pHit, double hHit, double vBlr, 
        r3_t *sNrm, int32_t NC, float sVal[]
      )
      { if ((wr_ray != NULL) && ((iPix->c[0] != iPix_wr.c[0]) || (iPix->c[1] != iPix_wr.c[1])))
          { fclose(wr_ray);
            iPix_wr = (i2_t){{ -1, -1 }};
            wr_ray = NULL;
          }
        if (wr_ray == NULL)
          { char *fname_ray = jsprintf("%s/pixel-rays-%04d-%04d.txt", frameFolder, iPix->c[0], iPix->c[1]);
            wr_ray = open_write(fname_ray, TRUE);
            free(fname_ray);
            iPix_wr = *iPix;
          }
        multifok_raytrace_show_ray_data(stderr, 8, iPix, pixSize, iSmp, step, wSmp, pRay, dRay, wRay, pHit, hHit, vBlr, sNrm, NC, sVal);
        multifok_raytrace_write_ray_data(wr_ray, iPix, pixSize, iSmp, step, wSmp, pRay, dRay, wRay, pHit, hHit, vBlr, sNrm, NC, sVal);
      }
  }

#define mfmi_image_gamma 1.000
  /* Assumed encoding gamma of input and output images. */

#define mfmi_image_bias 0.0327
  /* Assumed encoding bisa of input and output images. */

float_image_t *mfmi_read_pattern_image(char *fname)
  {
    image_file_format_t ffmt = image_file_format_PNG;
    bool_t yUp = FALSE;
    double gammaDec[1];
    double bias[1];
    bool_t verbose = TRUE;
    float_image_t *pat = float_image_read_gen_named(fname, ffmt, yUp, 0.0, 1.0, NULL, gammaDec, bias, verbose);
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
                  { r3_t ctr; box_center(3, objk->bbox, ctr.c);
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
            i2_vec_expand(&iDeb, (vec_index_t)NQ);
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
    char *stackFolder
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
            char *fname = jsprintf("%s/pixel-data-%04d-%04d.txt", stackFolder, ix, iy);
            FILE *wr = open_write(fname, TRUE);

            for (uint32_t ki = 0;  ki < NI; ki++)
              { multifok_frame_t *fri = stack->frame[ki];
                if (isfinite(fri->zDep))
                  { float hAvg = float_image_get_sample(fri->hAvg, 0, ix, iy);
                    float hDev = float_image_get_sample(fri->hDev, 0, ix, iy);
                    float shrp = float_image_get_sample(fri->shrp, 0, ix, iy);
                    r3_t sNrm;
                    for (int32_t j = 0;  j < 3; j++)
                      { sNrm.c[j] = float_image_get_sample(fri->sNrm, j, ix, iy); }
                    float sVal[3];
                    for (int32_t ic = 0;  ic < 3; ic++)
                      { sVal[ic] = (float)(ic < NC ? float_image_get_sample(fri->sVal, ic, ix, iy) : 0.0); }
                    r3_t ctr = (r3_t){{ ix+0.5, iy+0.5, 0.0 }};
                    r3_t pCtr = multifok_scene_coords_from_image_coords(&ctr, NX, NY, scene->dom, fri->zFoc);
                    multifok_raytrace_write_pixel_data
                      ( wr, iPix, pixSize, &pCtr, fri->zFoc, fri->zDep,
                        shrp, hAvg, hDev, &sNrm, NC, sVal
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
    if (strlen(o->sceneType) != 1)
      { argparser_error(pp, "invalid \"-sceneType\""); }
    
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

    argparser_get_keyword(pp, "-gloss");
    o->gloss = argparser_get_next_double(pp, 0.0, 1.0);

    argparser_get_keyword(pp, "-stackFolder");
    o->stackFolder = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
