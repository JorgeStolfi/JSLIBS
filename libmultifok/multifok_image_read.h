/* Test tools for reading images related to multifocus stereo. */
/* Last edited on 2025-04-11 08:52:38 by stolfi */

#ifndef multifok_image_read_H
#define multifok_image_read_H

#include <stdint.h>

#include <bool.h>
#include <float_image.h>

#include <multifok_image.h>
 
/* FRAME-SPECIFC PNG FILES */

float_image_t *multifok_image_read_scene_view(char *frameFolder);
  /* Reads the simulated or real scene view image {sVal} from file "{frameFolder}/sVal.png".
    It may be color (3 channels) or grayscale (1 channel).
    The sample values will be to be in {[0 _ 1]}.   */ 

float_image_t *multifok_image_read_height_average(char *frameFolder, double hMin, double hMax);
  /* Reads the scene average {Z} image {hAvg} from file "{frameFolder}/hAvg.png"
    and scales the samples from {[0_1]} to {[hMin _ hMax]} */ 

float_image_t *multifok_image_read_height_deviation(char *frameFolder, double dMax);
  /* Reads the scene Z deviation image {hDev} from file "{frameFolder}/hDev.png"
    and scales the samples from {[0_1]} to {[0 _ dMax]}. */ 

float_image_t *multifok_image_read_normal_average(char *frameFolder);
  /* Reads the scene average normal vector image {sNrm} from the color
    file "{frameFolder}/sNrm.png". The red, green, and blue values are
    converted from {[0_1]} to {[-1 _ +1]} and used as {X}, {Y}, and {Z}
    coordinates, respectively. */ 

float_image_t *multifok_image_read_sharpness(char *frameFolder);
  /* Reads the actual sharpness image {shrp} from file "{frameFolder}/shrp.png".
    Sample values will be in {[0 _ 1]}. */ 
 
/* FRAME-SPECIFC FNI FILES */

float_image_t *multifok_image_read_fni_height_average(char *frameFolder);
  /* Reads from file "{frameFolder}/hAvg.fni" the scene average {Z} image {hAvgWht}. 
    It must have two channels, channel 1 being the reliability weight. */ 

float_image_t *multifok_image_read_fni_normal_average(char *frameFolder);
  /* Reads from file "{frameFolder}/sNrm.fni" the scene average {Z} image {sNrmWht}.  
  It must have four channels, channel 3 being the reliability weight. */ 
   
/* GENERIC IMAGE FILE READING */

float_image_t *multifok_image_read_png(char *dir, char *name, float vMin, float vMax);
  /* Reads an image from file "{dir}/{name}.png". It may have any number of  channels. 
    After reading, all samples are remapped from {[0 _ 1]} to {[vMin _ vMax]}. */

float_image_t *multifok_image_read_fni(char *dir, char *name, int32_t chns);
  /* Reads a float image from file "{dir}/{name}.fni" in FNI format, 
    without any conversion.  If {chns} is non-negative, the
    image must have either {chns} or {chns+1} channels. */

#endif
