/* Non-linear calibration refinement step - version 2. */
/* Last edited on 2022-10-20 05:53:34 by stolfi */

#ifndef tf_calib_refine2_H
#define tf_calib_refine2_H

#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include <values.h>
#include <r3.h> 
#include <r2.h> 
#include <r4x4.h> 
#include <tf_camera.h> 
#include <tf_calib.h> 

void tf_calib_refine2_gather_optimization_params
  ( int32_t nparams,
    double params[],
    tf_camera_specs_t *cspec,
    tf_calib_data_t * cdat,
    tf_camera_params_t *cpar,
    tf_optimization_choice_t *which );

void tf_calib_refine2_scatter_optimization_params
  ( int32_t nparams,
    double params[],
    tf_camera_specs_t *cspec,
    tf_calib_data_t * cdat,
    tf_camera_params_t *cpar,
    tf_optimization_choice_t *which );

double tf_calib_refine2_calc_Tz
  ( double zc_star,
    tf_calib_data_t * cdat,
    tf_camera_params_t *cpar );

double tf_calib_refine2_calc_zc_star 
  ( tf_calib_data_t * cdat,
    tf_camera_params_t *cpar );

#endif
