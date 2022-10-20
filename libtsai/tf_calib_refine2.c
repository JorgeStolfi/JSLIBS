/* See {fine2.h}. */
/* Last edited on 2022-10-20 05:53:42 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <jsfile.h>
#include <affirm.h>
#include <tf_camera.h>
#include <tf_calib.h>
#include <tf_matrix.h>
#include <tf_math.h>
#include <tf_calib.h>
#include <tf_calib_refine2.h>
 
/* !!! A conversao (f,Tz)<-->(logf,L) deveria estar em tf_camera_opt ou vice-versa !!! */

void tf_calib_refine2_gather_optimization_params
  ( int32_t nparams,
    double params[],
    tf_camera_specs_t *cspec,
    tf_calib_data_t * cdat,
    tf_camera_params_t *cpar,
    tf_optimization_choice_t *which ) 
{
    int32_t pos_param = 0;
    if (which->R) {
      r3_t R;
      tf_camera_matrix_to_euler_angles(&(cpar->S), &R);
      params[pos_param++] = R.c[0];
      params[pos_param++] = R.c[1];
      params[pos_param++] = R.c[2];
    }

    bool_t which_Tx_Ty = which->v_w[0] | which->v_w[1] | which->v_w[2] | which->R;

    if (which_Tx_Ty) {
      params[pos_param++] = cpar->S.c[1][0];
      params[pos_param++] = cpar->S.c[2][0];
    }
    
    if (which->f) {
      params[pos_param++] = tf_camera_safe_log(cpar->f);
    }

    bool_t which_L = which->v_w[0] | which->v_w[1] | which->v_w[2] | which->R | which->f;

    if (which_L) {
      double zc_star = tf_calib_refine2_calc_zc_star (cdat, cpar);
      params[pos_param++] = tf_camera_safe_log(zc_star/cpar->f);
    }
    if (which->kappa) {
       params[pos_param++] = tf_camera_stretch_param (cpar->kappa, &cspec->kappa);
    }
    if (which->sx) {
       params[pos_param++] = tf_camera_stretch_param (cpar->sx, &cspec->sx);
    }
    if (which->Cx) {
       params[pos_param++] = tf_camera_stretch_param (cpar->Cx, &cspec->Cx);
    }
    if (which->Cy) {
       params[pos_param++] = tf_camera_stretch_param (cpar->Cy, &cspec->Cy);
    }

    assert(pos_param == nparams);
}

void tf_calib_refine2_scatter_optimization_params
  ( int32_t nparams,
    double params[],
    tf_camera_specs_t *cspec,
    tf_calib_data_t * cdat,
    tf_camera_params_t *cpar,
    tf_optimization_choice_t *which ) 
{ 
  int32_t pos_param = 0;
  if (which->R) {
    r3_t R;
    R.c[0] = params[pos_param++];
    R.c[1] = params[pos_param++];
    R.c[2] = params[pos_param++];
    tf_camera_matrix_from_euler_angles(&R, &(cpar->S));
  }

  bool_t which_Tx_Ty = which->v_w[0] | which->v_w[1] | which->v_w[2] | which->R;

  if (which_Tx_Ty) {
    cpar->S.c[1][0] = params[pos_param++];
    cpar->S.c[2][0] = params[pos_param++];
    //fprintf(stderr,  "Tx %f, Ty = %f\n", cpar->S.c[1][0], cpar->S.c[2][0]);   
  }

  if (which->f) {
    cpar->f = tf_camera_safe_exp(params[pos_param++]);
  }

  bool_t which_L = which->v_w[0] | which->v_w[1] | which->v_w[2] | which->R | which->f;

  if (which_L) {
     double L = tf_camera_safe_exp(params[pos_param++]);
     double zc_star = L*cpar->f;
     cpar->S.c[3][0] = tf_calib_refine2_calc_Tz (zc_star, cdat, cpar);
  }

  if (which->kappa) {
     cpar->kappa = tf_camera_squeeze_param(params[pos_param++], &cspec->kappa);
  }
  if (which->sx) {
     cpar->sx = tf_camera_squeeze_param(params[pos_param++], &cspec->sx);
  }
  if (which->Cx) {
     cpar->Cx = tf_camera_squeeze_param(params[pos_param++], &cspec->Cx);
  }
  if (which->Cy) {
     cpar->Cy = tf_camera_squeeze_param(params[pos_param++], &cspec->Cy);
  }
 
  //fprintf(stderr,  "Numero de params final %d\n", pos_param);   
  assert(pos_param == nparams);
}

double tf_calib_refine2_calc_zc_star 
  ( tf_calib_data_t * cdat,
    tf_camera_params_t *cpar ) 
{

    int32_t i;
    double sum_mean = 0.0, sum_std = 0.0, sum_w = 0.0;

    double *z_v = (double *)malloc(cdat->np*sizeof(double));

    for (i = 0; i < cdat->np; i++) {
         /* map the known world coordinates to (predicted) undistorted sensor plane coordinates */
         r3_t p_c_predicted = tf_world_coords_to_camera_coords (cpar, cdat->world[i]);     
	 z_v[i] = p_c_predicted.c[2];
	 sum_mean += cdat->weight[i]*z_v[i];
	 sum_w += cdat->weight[i];
    }

    double mean = sum_mean/sum_w;

    for (i = 0; i < cdat->np; i++) {
         double v = ( z_v[i] - mean ); 
	 sum_std += cdat->weight[i] * (v * v);
    }

    double std = sqrt(sum_std/sum_w);

    double zc_star = mean + 0*std;

    return zc_star;
}


double tf_calib_refine2_calc_Tz
  ( double zc_star,
    tf_calib_data_t * cdat,
    tf_camera_params_t *cpar ) 
{
    int32_t i;
    double sum_mean = 0.0, sum_std = 0.0, sum_w = 0.0;
    double *dzc = (double *)malloc(cdat->np*sizeof(double));

    for (i = 0; i < cdat->np; i++) {

         double dzc_i = (cpar->S.c[3][1] * cdat->world[i].c[0] +
                         cpar->S.c[3][2] * cdat->world[i].c[1] +
                         cpar->S.c[3][3] * cdat->world[i].c[2] ); 

	 dzc[i] = dzc_i;      

         sum_mean += cdat->weight[i] * dzc[i];

  	 sum_w += cdat->weight[i];

    }

    double mean = sum_mean/sum_w;

    for (i = 0; i < cdat->np; i++) {

         double dzc_i = ( dzc[i] - mean ); 

         sum_std += cdat->weight[i] * (dzc_i * dzc_i);
    }

    double std = sqrt(sum_std/sum_w);

    double Tz = zc_star - (mean - 0*std);

    free(dzc); 

    return Tz;
}
