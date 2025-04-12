/* Stack of images with limited depth of focus. */
/* Last edited on 2025-04-10 18:40:45 by stolfi */

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
    as described in {multifok_FRAME_DIR_INFO} and
    {multifok_FRAME_FILES_INFO}.
    
    The height values in the {frame.hAvg} map are converted to sample values
    in the "hAvg.png" file by affinely scaling the range {[hMin _ hMax]}
    to {[0 _ 1]}.  
    
    The height deviation values in the {frame.hDev} map are converted to sample values
    in the "hDev.png" file by affinely scaling the range {[0_(hMax-hMin)/2]} to {[0_1]}. */
  
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
  "      \"sVal.png\" The simulated snapshot of the scene (in RGB color) with" \
  " depth-of-focus blur.  The R, G, and B samples are encoded linearly with " \
  " file samples {0..maxval} corresponding to {[0 _ 1]}.\n" \
  "\n" \
  "      \"hAvg.png\" A grayscale image with the average {hAvg} of the height of the" \
  " scene's surface visible within each pixel, accounting for" \
  " focus blur, measured perpendicularly from the {zFoc = 0} focusing plane" \
  " in the direction towards the camera.  In particular, if the frame's in-focus plane" \
  " is horizontal at {Z = zFoc} and the camera is looking down, the" \
  " height of the scene's sutface is is just its {Z} coordinate.\n" \
  "\n" \
  "      \"hAvg.fni\" Same as \"sNrm.png\" but in FNI format, with" \
  " two channels: heights in channel 0, with no conversion, and an" \
  " estimated reliability weight in {[0_1]} added as channel 1.\n" \
  "\n" \
  "      \"hDev.png\" A greyscale image with the standard deviation of the surface height" \
  " in each pixel (the RMS value of the height minus the average {hAvg}).\n" \
  "\n" \
  "      \"sNrm.png\" An RGB image with the average normal direction of the" \
  " portion of the scene's surface visible within each pixel.  The R, G, and B" \
  " channels are the {X}, {Y}, and {Z} coordinates of the normal vector, in the" \
  " scene coordinate system.  All three channels will be encoded linearly" \
  " with file samples {0..maxval} corresponding to coordinates {[-1 _ +1]}.\n" \
  "\n" \
  "      \"sNrm.fni\" Same as \"sNrm.png\" but in FNI format, with" \
  " four channels: normal vector components {X}, {Y}, {Z} in" \
  " channels 0..2, with no conversion, and an estimated reliability" \
  " weight in {[0_1]} added as channel 3.\n" \
  "\n" \
  "      \"shrp.png\" An indicator of how blurry is the image within" \
  " each pixel.  It is defined as {1/vBlr}, where {vBlr} is the RMS value" \
  " of the distance between the" \
  " points of the scene surface that are visible in" \
  " the pixel and the center of the pixel, measured in" \
  " a direction parallel to the in-focus plane.   The" \
  " samples are encoded linearly with " \
  " file samples {0..maxval} corresponding to {[0 _ 1]}."

#endif
