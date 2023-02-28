/* See {tf_calib_guess1.h}. */
/* Last edited on 2023-02-25 16:12:27 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <jsfile.h>
#include <rn.h>

#include <tf_camera.h>
#include <tf_matrix.h>
#include <tf_errors.h>
#include <tf_math.h>
#include <tf_calib.h>
#include <tf_calib_refine2.h>

#include <tf_calib_guess1.h>

void tf_calib_guess1_initial_camera_parameters 
  ( tf_calib_data_t * cdat, 
    tf_camera_specs_t *cspec, 
    r2_t p_i_dev, 
    tf_camera_params_t *cpar,
    bool_t debug
  )
{
  rm_t U;
  r3_t R; /*Rotation angles*/
  calibration_aux_data_t caux;
  
  tf_optimization_choice_t which;
  tf_select_no_parameters(&which);
  tf_select_all_variable_parameters(cspec, &which);
      
  /* Compute auxiliary data: */
  caux = tf_calib_guess1_create_auxiliary_calibration_aux_data_structure (cdat->np);
  tf_calib_guess1_compute_xd_yd_and_r_squared (cdat, cpar, caux);

  U = tf_calib_guess1_compute_U (cdat, caux);
  if (debug) {
    fprintf(stderr, "  Vector {U}\n");
    int32_t j;
    for (j = 0; j < 7; j++) { fprintf(stderr, "  U[%d] = %24.12f\n", j, U->c[j]);
      fprintf(stderr, "  \n");
    }
  }

  /* !!! should check whether {Tx,Ty,Tz} are fixed !!! */
  tf_calib_guess1_compute_Tx_and_Ty (cdat, cspec, cpar, caux, U);
  if (debug) {
    fprintf(stderr, "  Tx = %24.12f\n", cpar->S.c[1][0]);
    fprintf(stderr, "  Ty = %24.12f\n", cpar->S.c[2][0]);
  }
  
  if (which.sx)
    { tf_calib_guess1_compute_sx (cpar, U); }
  else
    { cpar->sx = LO(cspec->sx); }
  if (debug) {
    fprintf(stderr, "  sx = %24.12f\n", cpar->sx);
  }

  /* Recompute auxiliary data with new {sx}: */
  tf_calib_guess1_compute_xd_yd_and_r_squared (cdat, cpar, caux);

  /* !!! should check whether {R} is fixed !!! */
  tf_calib_guess1_compute_R_and_store_in_S (cpar, U, &R);
  if (debug) {
    fprintf(stderr, "  Rotation matrix {R}\n");
    tf_camera_matrix_print_rotation (&cpar->S, stderr);
  }

  /* !!! should check whether {f,Tz} are fixed !!! */
  tf_calib_guess1_compute_approximate_f_and_Tz (cpar, cdat, cspec, caux);
  if (debug) {
    fprintf(stderr, "  f =  %24.12f\n", cpar->f);
    fprintf(stderr, "  Tz = %24.12f\n", cpar->S.c[3][0]);
  }

  if (cpar->f < 0) {
    /*try the other solution for the orthonormal matrix*/
    cpar->S.c[1][3] = -cpar->S.c[1][3]; /*r3*/
    cpar->S.c[2][3] = -cpar->S.c[2][3]; /*r6*/
    cpar->S.c[3][1] = -cpar->S.c[3][1]; /*r7*/
    cpar->S.c[3][2] = -cpar->S.c[3][2]; /*r8*/
    if (debug) {
      fprintf(stderr, "  Reversing the {t} vector in {R}\n");
      tf_camera_matrix_print_rotation (&cpar->S, stderr);
    }

    tf_camera_matrix_to_euler_angles(&(cpar->S), &R);

    tf_calib_guess1_compute_approximate_f_and_Tz (cpar, cdat, cspec, caux);
    if (debug) {
      fprintf(stderr, "  f =  %24.12f\n", cpar->f);
      fprintf(stderr, "  Tz = %24.12f\n", cpar->S.c[3][0]);
    }

    if (cpar->f < 0) {
      fprintf (stderr, "error - possible handedness problem with data\n");
      //exit (-1);
    }
  } 

  tf_calib_guess1_compute_exact_f_Tz_kappa (cdat, cspec, p_i_dev, cpar);
  if (debug) {
    fprintf(stderr, "  f =     %24.12f\n", cpar->f);
    fprintf(stderr, "  Tz =    %24.12f\n", cpar->S.c[3][0]);
    fprintf(stderr, "  kappa = %24.12f\n", cpar->kappa);
  }

  /*disallocating structures*/
  tf_free_rm_structure (U);
  tf_calib_guess1_free_auxiliary_calibration_aux_data_structure (caux);
}

void tf_calib_guess1_compute_xd_yd_and_r_squared
  ( tf_calib_data_t * cdat, 
    tf_camera_params_t *cpar, 
    calibration_aux_data_t caux )
{   
  int32_t i;
  for (i = 0; i < cdat->np; i++) {
    caux->distorted_coords[i] = tf_image_coords_to_sensor_coords (cpar, cdat->image[i]);

    /*r_squared*/
    caux->r_squared[i] = SQR (caux->distorted_coords[i].c[0]) + SQR (caux->distorted_coords[i].c[1]);
  }
}

rm_t tf_calib_guess1_compute_U (tf_calib_data_t * cdat, calibration_aux_data_t caux)
{
  int32_t i;
  mat_rm_t M = tf_alloc_mat_rm (cdat->np, 7);
  rm_t b = tf_alloc_rm (cdat->np);
  rm_t U = tf_alloc_rm (7);

  for (i = 0; i < cdat->np; i++) {
    M->c[i*7]   =  caux->distorted_coords[i].c[1] * cdat->world[i].c[0];
    M->c[i*7+1] =  caux->distorted_coords[i].c[1] * cdat->world[i].c[1];
    M->c[i*7+2] =  caux->distorted_coords[i].c[1] * cdat->world[i].c[2];
    M->c[i*7+3] =  caux->distorted_coords[i].c[1];
    M->c[i*7+4] = -caux->distorted_coords[i].c[0] * cdat->world[i].c[0];
    M->c[i*7+5] = -caux->distorted_coords[i].c[0] * cdat->world[i].c[1];
    M->c[i*7+6] = -caux->distorted_coords[i].c[0] * cdat->world[i].c[2];
    b->c[i]     =  caux->distorted_coords[i].c[0];

    int32_t j;
    for (j = 0; j < 7; j++)
      fprintf(stderr, " %12.8f ", M->c[i*7+j]);
    fprintf(stderr, " = %12.8f\n", b->c[i]);    
  }

  tf_solve_system_mxn (M, U, b, cdat->weight);
  fprintf(stderr, "U solution\n");
  int32_t j;
  for (j = 0; j < 7; j++)
    fprintf(stderr, " %12.8f ", U->c[j]);
  fprintf(stderr, "\n");

  tf_free_mat_rm_structure (M);
  tf_free_rm_structure (b);

  return U;
}

void tf_calib_guess1_compute_Tx_and_Ty
( tf_calib_data_t * cdat, 
  tf_camera_specs_t *cspec, 
  tf_camera_params_t *cpar, 
  calibration_aux_data_t caux, 
  rm_t U )
{
  int32_t i, far_point;
  double Tx, Ty, Ty_squared, x, y;
  double distance, far_distance;
    
  /* first find the square of the magnitude of Ty */
  Ty_squared = 1 / (SQR (U->c[4]) + SQR (U->c[5]) + SQR (U->c[6]));

  /* find a point that is far from the image center */
  far_distance = 0;
  far_point = 0;
  for (i = 0; i < cdat->np; i++) {
    if ((distance = caux->r_squared[i]) > far_distance) {
      far_point = i;
      far_distance = distance;
    }
  }  
   
  /* now find the sign for Ty, start by assuming Ty > 0 */
  Ty = sqrt (Ty_squared);
  Tx = U->c[3] * Ty;

  x = (U->c[0] * Ty) * cdat->world[far_point].c[0] + 
    (U->c[1] * Ty) * cdat->world[far_point].c[1] +
    (U->c[2] * Ty) * cdat->world[far_point].c[2] + Tx;

  y = (U->c[4] * Ty) * cdat->world[far_point].c[0] + 
    (U->c[5] * Ty) * cdat->world[far_point].c[1] +
    (U->c[6] * Ty) * cdat->world[far_point].c[2] + Ty;

  /* flip Ty if we guessed wrong */
  if ( (SIGNBIT (x) != SIGNBIT (caux->distorted_coords[far_point].c[0])) || 
       (SIGNBIT (y) != SIGNBIT (caux->distorted_coords[far_point].c[1])) )
    Ty = -Ty;

  /* update the calibration constants */
  cpar->S.c[1][0] = U->c[3] * Ty;    /*Tx*/
  cpar->S.c[2][0] = Ty;              /*Ty*/
  fprintf(stderr, "Tx %f Ty %f\n", cpar->S.c[1][0], cpar->S.c[2][0]);
}

void tf_calib_guess1_compute_sx (tf_camera_params_t *cpar, rm_t U)
{
  cpar->sx = sqrt (SQR (U->c[0]) + SQR (U->c[1]) + SQR (U->c[2])) * fabs (cpar->S.c[2][0]);
}

void tf_calib_guess1_compute_R_and_store_in_S (tf_camera_params_t *cpar, rm_t U, r3_t *R)
{
  double r1, r2, r3, r4, r5, r6, r7;
  double sa, ca, sb, cb, sg, cg;

  r1 = U->c[0] * cpar->S.c[2][0] / cpar->sx;
  r2 = U->c[1] * cpar->S.c[2][0] / cpar->sx;
  r3 = U->c[2] * cpar->S.c[2][0] / cpar->sx;
  r4 = U->c[4] * cpar->S.c[2][0];
  r5 = U->c[5] * cpar->S.c[2][0];
  r6 = U->c[6] * cpar->S.c[2][0];

  /* use the outer product of the first two rows to get the last row */
  r7 = r2 * r6 - r3 * r5;

  /* now find the RPY angles corresponding to the estimated rotation matrix */
  R->c[2] = atan2 (r4, r1);

  sincos(R->c[2], &sg, &cg);

  R->c[1] = atan2 (-r7, r1 * cg + r4 * sg);
    
  R->c[0] = atan2 (r3 * sg - r6 * cg, r5 * cg - r2 * sg);

  sincos(R->c[0], &sa, &ca);
    
  sincos(R->c[1], &sb, &cb);

  fprintf(stderr, "Rx = %12.8f, Ry = %12.8f, Rz = %12.8f\n", R->c[0], R->c[1], R->c[2]);

  /* now generate a more orthonormal rotation matrix from the RPY angles */
  cpar->S.c[1][1] =  cb * cg;                  /*r1*/      
  cpar->S.c[1][2] =  cg * sa * sb - ca * sg;   /*r2*/
  cpar->S.c[1][3] =  sa * sg + ca * cg * sb;   /*r3*/
  cpar->S.c[2][1] =  cb * sg;                  /*r4*/ 
  cpar->S.c[2][2] =  sa * sb * sg + ca * cg;   /*r5*/
  cpar->S.c[2][3] =  ca * sb * sg - cg * sa;   /*r6*/ 
  cpar->S.c[3][1] = -sb;                       /*r7*/
  cpar->S.c[3][2] =  cb * sa;                  /*r8*/
  cpar->S.c[3][3] =  ca * cb;                  /*r9*/
 
}

void tf_calib_guess1_compute_approximate_f_and_Tz
  ( tf_camera_params_t *cpar, 
    tf_calib_data_t * cdat, 
    tf_camera_specs_t *cspec, 
    calibration_aux_data_t caux )
{
  /* !!! Should check whether {f} and/or {Tz} and/or {kappa} are fixed !!! */
  int32_t i;
  mat_rm_t M = tf_alloc_mat_rm (cdat->np, 2);
  rm_t b = tf_alloc_rm (cdat->np);
  rm_t U = tf_alloc_rm (2);

  for (i = 0; i < cdat->np; i++) {

    M->c[i*2] =  cpar->S.c[2][1] * cdat->world[i].c[0] + 
      cpar->S.c[2][2] * cdat->world[i].c[1] + 
      cpar->S.c[2][3] * cdat->world[i].c[2] + cpar->S.c[2][0];

    M->c[i*2+1] = -caux->distorted_coords[i].c[1];

    b->c[i]  = ( cpar->S.c[3][1] * cdat->world[i].c[0] + 
		 cpar->S.c[3][2] * cdat->world[i].c[1] + 
		 cpar->S.c[3][3] * cdat->world[i].c[2] ) * caux->distorted_coords[i].c[1];
     
    int32_t j;
    for (j = 0; j < 2; j++)
      fprintf(stderr, " %12.8f ", M->c[i*2+j]);
    fprintf(stderr, " = %12.8f\n", b->c[i]); 


  }

  tf_solve_system_mxn (M, U, b, cdat->weight);
  fprintf(stderr, "U solution\n");
  int32_t j;
  for (j = 0; j < 2; j++)
    fprintf(stderr, " %12.8f ", U->c[j]);
  fprintf(stderr, "\n"); 


  /* update the calibration constants */
  cpar->f = U->c[0];
  cpar->S.c[3][0] = U->c[1]; /*Tz*/
  cpar->kappa = 0.0; /* this is the assumption that our calculation was made under */

  tf_free_mat_rm_structure (M);
  tf_free_rm_structure (b); 
  tf_free_rm_structure (U); 
}

void tf_calib_guess1_compute_exact_f_Tz_kappa
 ( tf_calib_data_t * cdat,
   tf_camera_specs_t *cspec,
   r2_t p_i_dev,
   tf_camera_params_t *cpar )
{
    
  tf_optimization_choice_t which;
  tf_select_no_parameters (&which);
  which.f = TRUE;
  which.v_w[3] = TRUE;
  which.kappa = TRUE;

  tf_generic_optimization
    ( cdat, cspec, p_i_dev, FALSE, cpar, &which,
      tf_calib_refine2_gather_optimization_params,
      tf_calib_refine2_scatter_optimization_params,
      FALSE ); 
}

calibration_aux_data_t tf_calib_guess1_create_auxiliary_calibration_aux_data_structure (int32_t nmarks)
{
  calibration_aux_data_t caux = (calibration_aux_data_t)malloc(sizeof(struct _calibration_aux_data_t));   
   
  caux->nmarks = nmarks;
  caux->distorted_coords = (r2_t *)malloc(nmarks * sizeof(r2_t));
  caux->r_squared = rn_alloc(nmarks);

  return caux; 
}

void tf_calib_guess1_print_calibration_aux_data (calibration_aux_data_t caux, FILE *wr)
{
  int32_t i;
  fprintf(wr, "%s\n", "Calibration aux data:");
  for (i = 0; i < caux->nmarks; i++) {
    fprintf(wr, "caux->x[%d] -> %9.9f \n", i, caux->distorted_coords[i].c[0]);
    fprintf(wr, "caux->y[%d] -> %9.9f \n", i, caux->distorted_coords[i].c[1]);
    fprintf(wr, "caux->r_squared[%d] -> %9.9f \n", i, caux->r_squared[i]);
  }
  fprintf(wr, "\n");
}

void tf_calib_guess1_free_auxiliary_calibration_aux_data_structure (calibration_aux_data_t caux)
{
  free(caux->distorted_coords);
  free(caux->r_squared);
  free(caux);
}
