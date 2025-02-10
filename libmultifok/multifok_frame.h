/* Stack of images with limited depth of focus. */
/* Last edited on 2025-02-01 16:55:48 by stolfi */

#ifndef multifok_frame_H
#define multifok_frame_H

#include <stdint.h>

#include <float_image.h>
  
typedef struct multifok_frame_t 
  { int32_t NC;             /* Number of channels in the scene image {sVal}. */
    int32_t NX, NY;         /* Image dimensions of all images. */
    float_image_t *sVal;    /* Simulated camera image. */
    float_image_t *shrp;    /* Actual sharpness map. */
    float_image_t *hAvg;    /* Actual scene {Z} average map. */
    float_image_t *hDev;    /* Actual scene {Z} deviation map. */
    float_image_t *sNrm;    /* Surface normal direction. */
    double zFoc;            /* Nominal {Z} of focus plane in each image. */
    double zDep;            /* Nominal depth of focus. */
  } multifok_frame_t;
  /* A a simulated limited-focus camera view image {sVal} and associated
    data images {shrp,hAvg,hDev,sNrm}. All five images must have the
    same number of columns {NX} and of rows {NY}. The simulated image
    {sVal} must have {NC} channels, the normal image {sNrm} must have 3
    channels, while the others must have a single channel. */
     
multifok_frame_t *multifok_frame_from_images
  ( float_image_t *sVal,    /* Simulated camera image. */
    float_image_t *shrp,    /* Actual sharpness map. */
    float_image_t *hAvg,    /* Actual scene {Z} average map. */
    float_image_t *hDev,    /* Actual scene {Z} deviation map. */
    float_image_t *sNrm,    /* Surface normal direction. */
    double zFoc,            /* Nominal {Z} of focus plane in each image. */
    double zDep             /* Nominal depth of focus. */
  );
  /* Creates a {multifok_frame_t} record from the given data. Checks
    that all images have the same {NX,NY}, and all have one channel
    except {sNrm} that must have three channels and {sVal} which may
    have one or more. */

void multifok_frame_free(multifok_frame_t *frame);
  /* Releases all storage used by {frame} including the images and
    the {*frame} record itself. */
    
multifok_frame_t *multifok_frame_read
  ( char *frameFolder,
    bool_t gray,
    double zFoc,
    double zDep,
    double hMin,
    double hMax
  );
  /* Reads a set of images {frame.{sVal,shrp,hAvg,hDev,sNrm}} from files
    "{frameFolder}/{sub}.png" where {sub} is "sVal", "shrp", 
    "hAvg", "hDev", and "sNrm", respectively.
    
    The samples in the {hAvg} and {hDev} files are implicit converted
    from their {0..maxval} range to {[hMin _ hMax]}. Typically
    {hMin} and {hMax} are the scene's min and max {Z} coordinates.
    
    The {sNrm} image will have the file samples in {0..maxval} mapped to the
    range {[-1 _ +1]}.
    
    The {shrp} image will have samples in the range {[0 _ 1]} as implied
    by the file.
    
    If {gray} is true, the camera view image {sVal} is converted to grayscale
    (with a single channel) as it is read.
    
    Returns those images in a {multifok_frame_t} object. */

void multifok_frame_write
  ( multifok_frame_t *frame,
    char *frameFolder,
    double hMin,
    double hMax
  ); 
  /* Writes the images of the {frame} to the directory "{frameFolder}"
    as described in {multifok_FRAME_DIR_INFO} and {multifok_FRAME_FILES_INFO}. */
  
#define multifok_FRAME_DIR_INFO \
  "The frame sub-folder name is" \
  " generally \"zf{FFF.FFFF}-df{DDD.DDDD}\", where" \
  " {FFF.FFFF} and {DDD.DDDD} are the values of the frame's in-focus plane" \
  " position {zFoc} and nominal depth of focus {zDep}" \
  " fields, both formatted as \"%08.4d\".  As a special case, the" \
  " sub-folder \"sharp\" will have a frame without any focus" \
  " blurring ({zDep} infinite)."
  
#define multifok_FRAME_FILES_INFO \
  "Each frame sub-folder will contain four images:\n" \
  "\n" \
  "      \"sVal.png\" The simulated snapshot of the scene with" \
  " depth-of-focus blur.  The samples are encoded linearly with " \
  " file samples {0..maxval} corresponding to {[0 _ 1]}.\n" \
  "\n" \
  "      \"sNrm.png\" The average normal direction of the" \
  " portion of the scene's surface visible within each pixel, in scene" \
  " coordinates.  All three channels will be encoded linearly" \
  " with file samples {0..maxval} corresponding to {[-1 _ +1]}.\n" \
  "\n\n" \
  "\n" \
  "      \"zAgv.png\" The average {hAvg} of the height of the" \
  " scene's surface visible within each pixel, accounting for" \
  " focus blur, measured perpendicularly to the in-focus plane" \
  " with position {zFoc = 0}, in the direction towards" \
  " the camera.  In particular, if the frame's in-focus plane" \
  " is horizontal at {Z = zFoc} and the camera is looking down, the" \
  " height of the scene's sutface is is just its {Z} coordinate.\n" \
  "\n" \
  "      \"hDev.png\" The standard deviation of the surface height" \
  " in each pixel (the RMS value of the height minus the average {hAvg}).\n" \
  "\n" \
  "      \"shrp.png\" An indicator of how blurry is the image within" \
  " each pixel.  It is defined as {1/vBlr}, where {vBlr} is the RMS value" \
  " of the distance between the" \
  " points of the scene surface that are visible in" \
  " the pixel and the center of the pixel, measured in" \
  " a direction parallel to the in-focus plane.   The" \
  " samples are encoded linearly with " \
  " file samples {0..maxval} corresponding to {[0 _ 1]}.\n" \
  "\n" \
  "    The pixel values of the {hAvg} and {hDev} images are" \
  " implicitly scaled so that some range {[hMin _ hMax]} is mapped" \
  " to {[0 _ 1]}.  Typically {hMin} and {hMax} are the nominal" \
  " min and max {Z} coordinates of the scene's surface."

#endif
