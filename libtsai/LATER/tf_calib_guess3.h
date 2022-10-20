/* Last edited on 2022-10-20 05:57:44 by stolfi */
/* Redesigned Tsai calibration routines (version 2). */

#ifndef tf_calib_guess3_H
#define tf_calib_guess3_H

#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include <r3.h> 
#include <r2.h> 
#include <r4x4.h> 
#include <rmxn.h> 
#include <tf_camera.h>
#include <tf_calib.h>
#include <tf_calib_guess1.h>

/* MAIN CALIBRATION GUESSING PROCEDURE */

void tf_calib_guess3_initial_camera_parameters 
  ( tf_calib_data_t * cdat, 
    tf_camera_specs_t *cspec,
    bool_t distance_weights,
    tf_camera_params_t *cpar,
    bool_t affine_only, 
    bool_t debug );
/* 
  This routine uses yet another method to compute the initial guess
  for the non-linear optimizer.  

  This routine ignores the intervals in {cspec} and tries
  to fit a general projective model to the calibration data points.

  The procedure currently works only if the world positions of the
  marks are non-coplanar and there are enough data points.

  The procedure uses the fields {kappa,Npx,Npy,dpx,dpy,sx,Cx,Cy} 
  of {cpar}. 

  If {distance_weights} is FALSE, each data point is weighted with
  {cdat->weight[mark]}.  If {distance_weights} is TRUE, the procedure
  uses {cpar->S} to compute the rough distance of each mark from the
  camera, and divides {cdat->weight[mark]} by that distance.  This
  option should have a (slight) effect only when the camera is very
  close to some marks. */

/* AUXILIARY ROUTINES */

void tf_calib_guess3_compute_f_of_fixed_camera
  ( int32_t n,
    r3_t p_w[],
    r3_t b_w,
    r2_t p_u[],
    r2_t b_u,
    double weight[], 
    double *f,
    bool_t debug);
/* 
  Special case of {tf_calib_guess3_initial_camera_parameters}
  when the camera has fixed position and orientation.
  Not yet implemented. */


void tf_calib_guess3_compute_S_f_of_fixed_position_camera
  ( int32_t n,
    r3_t p_w[],
    r3_t b_w,
    r2_t p_u[],
    r2_t b_u,
    double weight[],
    r3_t *v_w,  /* Position of camera. */
    bool_t affine_only, 
    bool_t debug,
    r4x4_t *S,  
    double *f);
/* 
  Special case of {tf_calib_guess3_initial_camera_parameters}
  when the camera has fixed position but unknown orientation.
  Not yet implemented. */

void tf_calib_guess3_compute_S_f_of_unrestricted_camera
  ( int32_t nmarks,
    r3_t p_w[],
    r3_t b_w,
    r2_t p_u[],
    r2_t b_u,
    double weight[],
    bool_t variable_R,
    bool_t affine_only, 
    bool_t debug, 
    r4x4_t *S,
    double *f );
  /* ??? */

void tf_calib_guess3_compute_S_f_of_far_away_camera
  ( int32_t n,
    r3_t p_w[],
    r3_t b_w,
    r2_t p_u[],
    r2_t b_u,
    double weight[],
    bool_t variable_R,
    bool_t affine_only, 
    bool_t debug,
    r4x4_t *S,
    double *f);
/* 
  Special case of {tf_calib_guess3_initial_camera_parameters}
  when the camera has faraway position but unknown orientation.

  Computes the camera projection parameters {cpar->S} and {cpar->f} that yield the
  best match between the given undistorted projected coordinates {p_u[0..n-1]}
  (whose barycenter is {b_u}) and the undistorted projected coordinates {p_u'[0..n-1]}
  computed from the world coordinates {p_w[k]} with that {S} and {f}; 
  constrained to the camera being very far away from the world barycenter {b_w}.

  If {variable_R} is FALSE, the procedure assumes that the rotation 
  sub-matrix of {S} is fixed, and computes only the elements {Tx,Ty,Tz} and {f}.
  Otherwise it computes {R} too. */

void tf_calib_guess3_compute_R_mu
  ( int32_t n,
    r3_t p_w[],
    r3_t b_w,
    r2_t p_u[],
    r2_t b_u,
    double weight[],
    bool_t affine_only, 
    bool_t debug,
    r4x4_t *S,
    double *mu );
/* Computes an approximate rotation matrix {R} and magnification
  factor {*mu} assuming that the camera is at infinity.  The {*mu}
  factor is the approximate mean value of {f/z_c} in the
  camera-to-undistorted mapping formula.  The matrix {R} is 
  stored as the rotation submatrix of {S}; element {S.c[0][0]} is set to 1,
  and the other elements are set to 0  */

void tf_calib_guess3_compute_mu_given_R
  ( int32_t n,
    r3_t p_w[],
    r3_t b_w,
    r2_t p_u[],
    r2_t b_u,
    double weight[],
    bool_t debug,
    r4x4_t *S,
    double *mu );
/* Computes an approximate magnification factor {*mu} assuming that
  {cpar} has the correct rotation matrix and the camera is at
  infinity.  The {*mu} factor is the approximate mean value of
  {f/z_c} in the camera-to-undistorted mapping formula. */

void tf_calib_guess3_compute_initial_affine_model
  ( int32_t n,
    r3_t p_w[],
    r3_t b_w,
    r2_t p_u[],
    r2_t b_u,
    double weight[],
    bool_t debug,
    r3_t *rh,
    r3_t *sh,
    r2_t *Th );
/* Computes an affine approximation {H d_w[k] + Th} to the mapping {d_w[k] --> (d_u[k].x,d_u[k].y)},
  where {d_w[k] = p_w[k] - b_w}, {d_u[k] = p_u[k] - b_u}, {H} is a {3 x 2} matrix with 
  rows {rh,sh}, and {Th} is a 2-vector.  */

void tf_calib_guess3_extract_camera_vectors
  ( r3_t rh, 
    r3_t sh, 
    bool_t affine_only,
    bool_t debug,
    r3_t *r, 
    r3_t *s, 
    r3_t *t,
    double *mu );
/* Computes the camera system axis vectors {*r,*s,*t} and the approximate
  magnification factor {*mu} from the two rows {rh,sh} of the
  matrix {H} computed by {tf_calib_guess3_initial_affine_model}. */

void tf_calib_guess3_compute_Tx_given_R_mu (int32_t n, r3_t b_w, r2_t b_u, r4x4_t *S, double mu, bool_t debug);

void tf_calib_guess3_compute_Ty_given_R_mu (int32_t n, r3_t b_w, r2_t b_u, r4x4_t *S, double mu, bool_t debug);

#endif
