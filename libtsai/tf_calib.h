/* Last edited on 2022-10-20 05:55:12 by stolfi */

#ifndef tf_calib_H
#define tf_calib_H

/* General-purpose data and functions for calibration routines. */
#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include <values.h>
#include <r3.h> 
#include <r2.h> 
#include <r4x4.h> 
#include <tf_camera.h>
#include <tf_camera_specs.h>
#include <tf_calib_data.h>

/* ---------------------------------------------------------------------- */
/* AUXILIARY DATA TYPES */

typedef struct tf_optimization_choice_t
{
  bool_t R;
  bool_t v_w[3];
  bool_t f;
  bool_t kappa;
  bool_t sx;
  bool_t Cx;
  bool_t Cy;
} tf_optimization_choice_t;
/* A record that specifies which camera parameters in {cpar} are
  to be optimized. */

typedef void tf_gather_params_proc_t
  ( int32_t nparams,
    double params[],
    tf_camera_specs_t *cspec,
    tf_calib_data_t * cdat,
    tf_camera_params_t *cpar,
    tf_optimization_choice_t *which );
/* 
  Type of a procedure that extracts from {cpar} the fields selected by
  {which} and packs them as {params[0..nparams-1]}.  The option
  {which->R} should extract the rotation part of the matrix
  {cpar->S} and pack it as the three equivalent euler angles. */

typedef void tf_scatter_params_proc_t
  ( int32_t nparams,
    double params[],
    tf_camera_specs_t *cspec,
    tf_calib_data_t * cdat,
    tf_camera_params_t *cpar,
    tf_optimization_choice_t *which );
/* 
  Type of a procedure that unpacks the parameters
  {params[0..nparams-1]} and stores them into the fields of {cpar}
  selected by {which}.  The option {which->R} should unpack three
  euler angles, convert them to an orthonormal matrix, and store 
  the latter into the rotation part of the matrix {cpar->S}. */
  
typedef void tf_guess_params_proc_t
  ( tf_calib_data_t * cdat, 
    tf_camera_specs_t *cspec, 
    tf_optimization_choice_t *which,   
    tf_camera_params_t *cpar,
    bool_t affine_only,
    bool_t debug);
/* 
  Type of a procedure that generates an initial guess
  for the camera parameter record {*cpar}, given the camera
  specification {cspec} and the fixed/variable flags {which}.
  If {affine_only} is true produces an affine transformation instead
  of a Euclidean transformation.*/

/* ---------------------------------------------------------------------- */
/* MAIN CALIBRATION PROCEDURE */

void tf_calib_generic
  ( tf_calib_data_t * cdat, 
    tf_camera_specs_t *cspec, 
    r2_t p_i_dev,
    bool_t use_cpar_dev_terms,
    tf_camera_params_t *cpar,
    tf_gather_params_proc_t gather_params,
    tf_scatter_params_proc_t scatter_params,
    bool_t debug
  );
/*
  This procedure calibrates all camera parameters in {cpar} that are
  variable according to {cspec}, using the data contained in {cdat}.
  The procedure automatically adapts the fitting according to the
  set of variable parameters (those that have non-trivial range in {cspec}).
  
  The data consists of the measured world coordinates
  {cdat->world[k]} (mm) and the measured image coordinates
  {cdat->image[k]} (pixels) of a certain set of reference
  points (/marks/) fixed on the scene, where {k} ranges in
  {0..cdat->np}. Note that each combination of variable parameters
  require a certain minimum number of data pairs, possibly
  non-coplanar.  
  
  The procedure uses a confidence weight {cdat->weight[k]} asociated
  to each data pair, and an in thes.

  An initial guess must be given in the {cpar} parameter.  This guess
  is refined by non-linear minimization of the
  {tf_calib_generic_optimization_error} goal function, parameterized
  by the {cspec} intervals, the given data weights
  {cdat->weight[0..cdat->np-1]}, and image coordinate
  uncertainties {p_i_dev}, using {gather_params/scatter_params}
  procedures to map the variable fields of {cpar} to/from an unbounded
  param vector.  The refined parameters are returned in {cpar} itself.
 */

/* ---------------------------------------------------------------------- */
/* GENERIC CAMERA PARAMETER OPTIMIZATION */

void tf_show_values_of_selected_parameters
  ( tf_camera_params_t *cpar,
    tf_optimization_choice_t *which, 
    FILE *ferr );
/* Prints to {ferr} all parameters of {cpar} that
   are variable according to {which}.  If {cpar} is null,
   prints only their names. */

void tf_show_which_parameters_are_selected(tf_optimization_choice_t *which, FILE *ferr);
/* Prints to {ferr} the names of all parameters that are selected by {which}. */

void tf_select_no_parameters(tf_optimization_choice_t *which);
/* Sets all fields {which.XXX} to FALSE. */

void tf_select_all_parameters(tf_optimization_choice_t *which);
/* Sets all fields {which.XXX} to TRUE. */

void tf_select_all_variable_parameters(tf_camera_specs_t *cspec, tf_optimization_choice_t *which);
/* Sets {which.XXX} to true only if the camera parameter {cspec.XXX} is variable.
  Leaves the other fields unchanged. */

void tf_unselect_all_fixed_parameters(tf_camera_specs_t *cspec, tf_optimization_choice_t *which);
/* Sets {which.XXX} to FALSE if the camera parameter {cspec.XXX} is fixed.
  Leaves the other fields unchanged. */

int32_t tf_generic_optimization_count_params 
  ( tf_optimization_choice_t *which );
/* Computes the number of {lmdif} parameters that are needed to 
  calibrate the parameters selected by {which}. */

int32_t tf_generic_optimization_count_error_terms 
  ( tf_calib_data_t * cdat,
    bool_t use_cpar_dev_terms,
    tf_optimization_choice_t *which );
/* Computes the number of {lmdif} error terms that arise when
  optimizing the parameters selected by {which}, given the calibration
  data {cdat}.  Assumes that the selected parameters are all variable.
  If {use_cpar_dev_terms} is TRUE, assumes also that those parameters have Kalman error
  terms associated with them. */

void tf_generic_optimization_compute_error_terms
  ( double *err,
    int32_t nerr,
    tf_calib_data_t * cdat,
    tf_camera_specs_t *cspec,
    r2_t p_i_dev,
    bool_t use_cpar_dev_terms,
    tf_camera_params_t *cpar,
    tf_optimization_choice_t *which );
/* Sets {err[0..n-1]} to the error terms that arise when optimizing
  the parameters selected by {which}, given the calibration data {cdat}.
  Assumes that the selected parameters are all variable and have 
  Kalman error terms associated with them.

  Each error term is computed by {(val - avg)/dev}, where {val}
  is the value of some quantity determined by {cpar}, {avg}
  is its expected value, and {dev} is its deviation.   The quantities
  are: 

     (a) mark position error terms: {avg} is the X or Y image
     coordinate of a mark {i}, as predicted from the world coordinates
     {cdat->world[i]} with parameters {cpar} (including radial
     distortion); {val} is the corresponding coordinate
     {cdat->image[i]}; and {dev} is {p_i_dev.c[axis]/sqrt(cdat->weight[i])}.

     (b) parameter strangeness terms: {val} is the 
     value of that parameter in {cpar}; {avg} is the midpoint of 
     that parameter's range in {cspec}; {dev} is 1/3 of the 
     half-width of that interval.  These terms are included
     only if {use_cpar_dev_terms} is true.
     
  For the {f} parameter, the computation is done with {log(f)} instead of {f}
  itself.  For the {Vx_w,Vy_w,Vz_w} parameters, the value {val} is taken from 
  elements {[1][0],[2][0]}, and {[3][0]} of the inverse of the matrix {S}.
  For the {R} parameter, rotation part of {S} is converted to Euler angles,
  which are compared to the corresponding ranges in {cspec}. */

void tf_generic_optimization
  ( tf_calib_data_t * cdat,
    tf_camera_specs_t *cspec,
    r2_t p_i_dev,
    bool_t use_cpar_dev_terms,
    tf_camera_params_t *cpar,
    tf_optimization_choice_t *which,
    tf_gather_params_proc_t gather_params,
    tf_scatter_params_proc_t scatter_params,
    bool_t debug
  );
/* Optimizes the fields of {cpar} which are selected by {which} and are
  variable in {cspec}, leaving the remaining fields unchanged.  
  The error metric is the weighted sum of squares of undistorted
  sensor coordinate mismatches, plus the {cpar}
  Kalman terms, computed by {tf_generic_optimization_error}

  At present, performs an unconstrained minimization on the 
  variable parameters, mapped to/from an unbounded real vector by
  {gather_params/scatter_params}.  */

void tf_clip_cpar_to_cspec
  ( tf_camera_specs_t *cspec,
    tf_camera_params_t *cpar,
    tf_optimization_choice_t *which );  
/* 
  Ensures that the parameters in {cpar} lie inside the ranges allowed
  by {cspec}.  If a parameter is variable according to {which}, the
  procedure merely checks that condition, and fails if it is false.
  If the parameter is variable, clips it to the specified range.
  
  This is a desperate measure intended to avoid NaNs and other nasty things.
  If this procedure needs to change any parameter, the situation is probably 
  hopeless already.  */

void tf_clip_to_range (char *name, double *v, bool_t variable, interval_t *range);
/* If {variable} is true, clips the value {*v} to the given interval {range}.
   Warns if the value changes. */
/* ---------------------------------------------------------------------- */
/* UTILITY FUNCTIONS */

void tf_compute_undistorted_obs_coordinates(int32_t n, r2_t p_i[], tf_camera_params_t *cpar, r2_t p_u[]);
/* Computes the observed undistorted sensor
  coordinates of all features.  Stores into {p_u[i]} the
  observed image coordinates {p_i[i]}, mapped to
  distorted sensor corodinates using
  {tf_image_coords_to_sensor_coords} with parameter
  {cpar->Cx,cpar->Cy,cpar->sx,cpar->dpx,cpar->dpy}, and then to
  undistorted sensor coordinates using
  {tf_dis_sensor_coords_to_und_sensor_coords} with
  {cpar->kappa}. */

void tf_compute_world_barycenter(int32_t n, r3_t p_w[], double weight[], r3_t *b_w);
/* Computes the weighted barycenter of the world coordinates of all
  {n} featuers.  Assumes that the coordinates of feature {i} are
  {p_w[i]}, and its weight is {weight[i]}. */
  
void tf_compute_undistorted_coords_barycenter(int32_t n, r2_t p_u[], double weight[], r2_t *b_u);
/* Computes the weighted barycenter of the undistorted sensor
  coordinates of all {n} features.  Assumes that the coordinates of
  feature {i} are {p_u[i]}, and its weight is {cdat->weight[i]}. */

void tf_show_optimization_errors
( double err[],
  double weight[],
  int32_t m,
  FILE *ferr );
/*  */

void tf_write_cpar_and_errors
  ( char *name, 
    tf_calib_data_t * cdat,
    tf_camera_specs_t *cspec,
    r2_t p_i_dev, 
    bool_t use_cpar_dev_terms,
    tf_camera_params_t *cpar,
    char *out_dir );
  /* If {cpar} is NULL, does nothing. */

void tf_generic_optimization_print_error_terms
  ( FILE *ferr,
    double *err,
    int32_t nerr,
    bool_t use_cpar_dev_terms,
    tf_calib_data_t * cdat,
    tf_optimization_choice_t *which );
  /* */

void tf_recompute_target_weights
  ( tf_camera_params_t *cpar,
    int32_t ntargets,
    r3_t p_w[],
    r2_t p_i[],
    double p_wgt[],
    r2_t dev_gud,
    r2_t dev_bad,
    bool_t debug );
/* Recomputes the weights {p_wgt[k]} of every mark {k} in {0..ntargets-1},
  based on the discrepancy between the given image coordinates {p_i[k]}
  and the image coordinates predicted from the world coordinates {p_w[k]}
  by the camera model {cpar}.  

  On input, assumes that the mark is a correct match (`inlier') with
  /a priori/ probability {p_wgt[k]}, and a bogus match (`outlier')
  otherwise. Assumes also that the found position {p_i[k]} for both
  inliers and outliers has a Gaussian distribution around the position
  predicted by {p_w} and {cpar} with standard deviations {dev_gud} and
  {dev_bad}, respectively.  On output, {p_wgt[k]} will be the /a
  posteriori/ probability of mark {k} being an inlier, taking into
  account its error. */

#endif
