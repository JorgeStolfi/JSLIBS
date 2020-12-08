/* Functions to summarize the camera calibration errors. */
/* Last edited on 2011-05-15 00:42:41 by stolfi */

#ifndef tf_errors_H
#define tf_errors_H

#include <stdio.h>
#include <ctype.h>
#include <r3.h> 
#include <r2.h> 
#include <r4x4.h> 
#include <tf_camera.h> 
#include <tf_calib.h> 

void tf_calib_summarize_errors (tf_calib_data_t * cdat, tf_camera_params_t *cpar, FILE *ferr);
  /* This routine prints summaries of the image-space and world-space errors of the 
    camera model {cpar}, for the calibration data {cdat}. */

void tf_calib_summarize_image_errors (tf_calib_data_t * cdat, tf_camera_params_t *cpar, FILE *ferr);
/* This routine prints to {ferr} the mean, standard deviation, max,
  and sum-of-squared errors in distorted image coordinates, as
  computed by {tf_camera_compute_image_error}, for the camera model
  {cpar} and calibration data {cdat}. */

void tf_calib_summarize_world_errors (tf_calib_data_t * cdat, tf_camera_params_t *cpar, FILE *ferr);
/* This routine prints to file {ferr} the mean, standard deviation, max, and sum-of-squared
* error of the distance of closest approach (i.e. 3D error) between the point   *
* in object space and the line of sight formed by back projecting the measured  *
* 2D coordinates out through the camera model.                                  *
* The calculation is for all of the points in the calibration data set.         *
 */
void tf_object_space_error_stats (  tf_calib_data_t * cdat, 
                                    tf_camera_params_t *cpar, 
                                    FILE *ferr  );
#endif
