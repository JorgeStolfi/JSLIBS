/* Tsai's algorithm for guessing the initial camera calibration. */
/* Last edited on 2011-05-17 02:01:17 by stolfi */

#ifndef tf_calib_guess1_H
#define tf_calib_guess1_H

#define _GNU_SOURCE
#include <stdio.h>
#include <ctype.h>
#include <r3.h> 
#include <r2.h> 
#include <r4x4.h> 
#include <rmxn.h> 

#include <tf_camera.h>
#include <tf_camera_specs.h>
#include <tf_matrix.h>
#include <tf_calib.h>

/* MAIN CALIBRATION GUESSING PROCEDURE */

void tf_calib_guess1_initial_camera_parameters 
  ( tf_calib_data_t * cdat, 
    tf_camera_specs_t *cspec, 
    r2_t p_i_dev,
    tf_camera_params_t *cpar,
    bool_t debug);
/*
  This routine stores into {cpar} approximate values for the camera parameters,
  intended as an initial guess for the general optimization.  It uses a variant
  of Tsai's guessing algorithm:

   @article{136938,                                                               
     author = {Roger Y. Tsai},                                                      
     title = {A versatile camera calibration technique for high-accuracy 3D machine 
            vision metrology using off-the-shelf TV cameras and lenses},          
     journal = {Journal of robotics and automation IEEE},                           
     year = {1987},                                                                 
     isbn = {0-86720-294-7},                                                        
     pages = {323--344},                                                            
     publisher = {IEEE},                                                            
     address = {, USA},                                                             
   }                                                                              

  The modifications made to that version are (1) automatic adaptation of the optimization steps
  according to the number of unknown parameters, and (2) the use of 
  a confidence weight {cdat->weight[k]} asociated to each data pair. */

/* AUXILIARY CALIBRATION DATA STRUCTURE */

typedef struct _calibration_aux_data_t {
    int    nmarks;           /* number of marks */
    r2_t   *distorted_coords;  /* computed image coordinates of each mark [mm] */
    double *r_squared;         /* squared distance from image center to each mark [mm] */
} *calibration_aux_data_t;

calibration_aux_data_t tf_calib_guess1_create_auxiliary_calibration_aux_data_structure (int nmarks);

void tf_calib_guess1_free_auxiliary_calibration_aux_data_structure (calibration_aux_data_t caux);

void tf_calib_guess1_print_calibration_aux_data (calibration_aux_data_t caux, FILE *wr);

/* AUXILIARY CALIBRATION FUNCTIONS */

void tf_calib_guess1_compute_xd_yd_and_r_squared
  ( tf_calib_data_t * cdat, 
    tf_camera_params_t *cpar, 
    calibration_aux_data_t caux );
/*
  This routine computes the distorted sensor-plane coordinates
  {caux->distorted_coords[i] = (Xd, Yd)} for each mark {i}, from its
  image coords {cdat->image[i] = (Xf, Yf)}.  Uses the
  {dpx,dpy,sx,Cx,Cy} camera parameters from {cdat}.  Does not use the
  {S,f,kappa} nor the world positions.  (see equation (6a) and (6b)
  (pp 328) in the article [1]). */

rm_t tf_calib_guess1_compute_U
  ( tf_calib_data_t * cdat, 
    calibration_aux_data_t caux);
/*
 This routine takes each 3D point and the corresponding image coordinate, and
 computes the U vector.  Uses initial estimates for the parameters {???}.
(see equation (16) (pp 333) in the article [1])
*/

void tf_calib_guess1_compute_Tx_and_Ty
  ( tf_calib_data_t * cdat, 
    tf_camera_specs_t *cspec, 
    tf_camera_params_t *cpar, 
    calibration_aux_data_t caux, 
    rm_t U );
/*
 This routine computes the approximate translation in x and y direction from the U vector.
  Uses initial estimates for the parameters {???}
 (see equations (17) and (iii) step (pp 331) in the article [1]) 
 */

void tf_calib_guess1_compute_sx (tf_camera_params_t *cpar, rm_t U);
/*
 This routine recomputes the X scale correction factor sx, given the
  previously computed tx,ty estimates, and the vector U.
 (see equations (17) and (iii) step (pp 331) in the article [1]) 
 */

void tf_calib_guess1_compute_R_and_store_in_S
  ( tf_camera_params_t *cpar, 
    rm_t U,  
    r3_t *R );
/*
 This routine computes the (3x3) rotation matrix
  from the U vector and previously computes parameters {???}.
 (see pp (333) in the article [1]) */

void tf_calib_guess1_compute_approximate_f_and_Tz 
  ( tf_camera_params_t *cpar, 
    tf_calib_data_t * cdat, 
    tf_camera_specs_t *cspec, 
    calibration_aux_data_t caux );

void tf_calib_guess1_compute_approximate_f_and_Tz_old
  ( tf_camera_params_t *cpar,
      tf_calib_data_t * cdat,
          tf_camera_specs_t *cspec,
	      calibration_aux_data_t caux );

/* This procedure computes the focal length f and the Z translation tz
  given previous estimates of {tx,ty,R,kappa,sx,Cx,Cy}. */

void tf_calib_guess1_compute_exact_f_Tz_kappa_error
  ( double *params, 
    double *err, 
    double *uwerr,
    tf_calib_data_t * cdat, 
    tf_camera_params_t *cpar, 
    calibration_aux_data_t caux );
/*
 The goal vector function for {tf_compute_exact_f_and_Tz}.
 It depends on three parameters {params[0..2]}, which 
 are assumed to be {f,Tz,kappa}, in that order.   
 */

 void tf_calib_guess1_compute_exact_f_Tz_kappa
 ( tf_calib_data_t * cdat,
   tf_camera_specs_t *cspec,
   r2_t p_i_dev,
   tf_camera_params_t *cpar ); 

/* This routine uses the lmdif package to compute better focal
  distance, Z translation and radial distortion {f,Tz,kappa} in
  {cpar}, assuming that the other parameters are fixed. */

#endif
