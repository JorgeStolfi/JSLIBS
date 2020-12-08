/* Calibration error functions. */
/* Last edited on 2011-05-15 00:46:53 by stolfi */

#ifndef tf_calib_quality_H
#define tf_calib_quality_H

#define _GNU_SOURCE
#include <stdio.h>
#include <tf_camera.h>
#include <tf_calib.h>

void tf_distorted_image_plane_error
  ( tf_calib_data_t * cdat, 
    tf_camera_params_t *cpar, 
    FILE *ferr );
/*
 This routine computes the mean, standard deviation, max, and sum-of-squared 
 error of the magnitude of the error, in distorted image coordinates, between  
 the measured location of a feature point in the image plane and the image of  
 the 3D feature point as projected through the calibrated model.               
 The calculation is for all of the points in the calibration data set.         
 */

void tf_undistorted_image_plane_error 
  ( tf_calib_data_t * cdat, 
    tf_camera_params_t *cpar, 
    FILE *ferr );
/*
 This routine computes the mean, standard deviation, max, and sum-of-squared 
 error of the magnitude of the error, in undistorted image coordinates, between
 the measured location of a feature point in the image plane and the image of  
 the 3D feature point as projected through the calibrated model.               
 The calculation is for all of the points in the calibration data set.         
*/

void tf_object_space_error_stats 
  ( tf_calib_data_t * cdat, 
    tf_camera_params_t *cpar, 
    FILE *ferr );
/*
 This routine computes the mean, standard deviation, max, and sum-of-squared 
 error of the distance of closest approach (i.e. 3D error) between the point   
 in object space and the line of sight formed by back projecting the measured  
 2D coordinates out through the camera model.                                  
 The calculation is for all of the points in the calibration data set.         
 */

void tf_normalized_calibration_error
  ( tf_calib_data_t * cdat, 
    tf_camera_params_t *cpar, 
    FILE *ferr  );
/*
 This routine computes an error measure proposed by Weng in IEEE PAMI,
 October 1992.  */

void tf_norm_calibration_error_measures 
  ( tf_calib_data_t * cdat, tf_camera_params_t *cpar, double *mean, double *stddev );
/*
  This routine computes the mean and stadard deviation of 
  the ??? error metric over all marks in {cdat}.  */

double tf_norm_calibration_point_error
  ( int point,
    tf_calib_data_t * cdat,
    tf_camera_params_t *cpar );
/*
  This routine computes the ??? error metric for the mark
  with number {point} in the data set {cdat}.  */

void tf_show_optimization_errors
( double err[],
  double weight[],
  int m,
  FILE *ferr );
/* Shows the (weighted) undistorted sensor coordinate 
   errors {err[0..m-1]} and the weights {weight[0..m/2-1]}. 
   Assumes that {err[2*i]} and {err[2*i+1]} are the signed X and Y errors
   of mark number {i}*/

#endif
