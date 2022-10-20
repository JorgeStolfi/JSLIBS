/* See {tf_calib_guess3.h}. */
/* Last edited on 2022-10-20 05:57:54 by stolfi */

#define _GNU_SOURCE
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

#include <jsfile.h>
#include <affirm.h>
#include <tf_camera.h>
#include <tf_matrix.h>
#include <tf_errors.h>
#include <tf_math.h>
#include <tf_calib.h>
#include <tf_calib_guess3.h>

void tf_calib_guess3_initial_camera_parameters 
  ( tf_calib_data_t * cdat, 
    tf_camera_specs_t *cspec, 
    bool_t distance_weights,
    tf_camera_params_t *cpar,
    bool_t affine_only, 
    bool_t debug)
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }

  tf_optimization_choice_t which;
  tf_select_no_parameters(&which);
  tf_select_all_variable_parameters(cspec, &which);
        
  int32_t nmarks = cdat->np;
  if (debug) {
    fprintf(stderr, "Computing the initial guess of {cpar}\n");
    tf_camera_specs_print (cspec, stderr);
    tf_show_which_parameters_are_selected(&which, stderr);
    tf_calib_data_print(cdat, stderr);
  }

  /* Some parameters must be known: */
  demand(! which.Cx, "cannot guess {Cx}");
  demand(! which.Cy, "cannot guess {Cy}");
  demand(! which.sx, "cannot guess {sx}");
 
  /* Compute the barycenter of world coordinates of all features: */
  r3_t *p_w = cdat->world;
  r3_t b_w; 
  tf_compute_world_barycenter(nmarks, p_w, cdat->weight, &b_w);

  /* Compute the observed undistorted coordinates of all marks features, and their barycenter: */
  assert(cpar->kappa >= LO(cspec->kappa) && cpar->kappa <= HI(cspec->kappa)); /* Paranoia. */
  r2_t *p_u = (r2_t *)notnull(malloc(nmarks*sizeof(r2_t)), "no mem");
  tf_compute_undistorted_obs_coordinates(nmarks, cdat->image, cpar, p_u);
  if (debug) {
    fprintf(stderr, "Undistorted projected mark coordinates:\n");
    tf_calib_data_print_image_points (cdat->np, p_u, stderr);
  }
  r2_t b_u;
  tf_compute_undistorted_coords_barycenter(nmarks, p_u, cdat->weight, &b_u);

  /* Compute the mark weights to use: */
  double* wt_loc = (double*)notnull(malloc(nmarks*sizeof(double)), "no mem");
  int32_t k;
  for (k = 0; k < nmarks; k++) {
    wt_loc[k] = cdat->weight[k];
    if (distance_weights) {
      /* Compute the camera-Z distance {zck} of mark from camera, using {cpar}: */
      double zck = cpar->S.c[3][0] + 
        cpar->S.c[3][1]*cdat->world[k].c[0] + 
        cpar->S.c[3][2]*cdat->world[k].c[1] + 
	cpar->S.c[3][3]*cdat->world[k].c[2];
      if (zck > 0.0) { wt_loc[k] /= zck; }
    }
  }

  /* Find the best-fitting camera at very large distance from {b_w}, ignoring the {f} range: */
  bool_t variable_R = TRUE; /* !!! CHECK !!! */
  tf_calib_guess3_compute_S_f_of_unrestricted_camera
    ( nmarks, p_w, b_w, p_u, b_u, wt_loc, variable_R, affine_only, debug, &(cpar->S), &(cpar->f) );
   
  /*disallocating structures*/
  free (p_u);
  free (wt_loc);
      
}

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
    double *f )
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }
  
  /* Initial least-squares fit.
     We look for 12 parameters

       P[0]  = Tz
       P[1]  = r7
       P[2]  = r8
       P[3]  = r9
       P[4]  = f*Tx
       P[5]  = f*r1
       P[6]  = f*r2
       P[7]  = f*r3
       P[8]  = f*Ty
       P[9]  = f*r4
       P[10] = f*r5
       P[11] = f*r6

     that provide best least-squares fit to the equations

      x_u[k]*(Tz) + x_u[k]*x_w[k]*r7 + x_u[k]*y_w[k]*r8 + x_u[k]*z_w[k]*r9 = 
        (f*Tx) + x_w[k]*(f*r1) + y_w[k]*(f*r2) + z_w[k]*(f*r3)

      y_u[k]*(Tz) + y_u[k]*x_w[k]*r7 + y_u[k]*y_w[k]*r8 + y_u[k]*z_w[k]*r9 = 
        (f*Ty) + x_w[k]*(f*r4) + y_w[k]*(f*r5) + z_w[k]*(f*r6)

     That is,

      [  E[0] E[1] ... E[11] ] P = 0;
      [  F[0] F[1] ... F[11] ] P = 0

    where

      E[0]  = x_u[k]             F[0]  = y_u[k]             
      E[1]  = x_u[k]*x_w[k]      F[1]  = y_u[k]*x_w[k] 
      E[2]  = x_u[k]*y_w[k]      F[2]  = y_u[k]*y_w[k] 
      E[3]  = x_u[k]*z_w[k]      F[3]  = y_u[k]*z_w[k] 
      E[4]  = -1                 F[4]  = 0  
      E[5]  = -x_w[k]            F[5]  = 0          
      E[6]  = -y_w[k]            F[6]  = 0          
      E[7]  = -z_w[k]            F[7]  = 0          
      E[8]  = 0                  F[8]  = -1                         
      E[9]  = 0                  F[9]  = -x_w[k]                    
      E[10] = 0                  F[10] = -y_w[k]                    
      E[11] = 0                  F[11] = -z_w[k]             

  */
       
  double P[8]; 
  tf_calib_guess3_compute_P(nmarks, p_w, b_w, p_u, b_u, distance_weights, weight, debug, P);

  /* Find the parameters {Tx,Ty.Tz,R,f} that best fit {P}: */
  tf_calib_guess3_extract_camera_parameters_from_P(rh, sh, affine_only, debug, S, f);
}

void tf_calib_guess3_compute_P
  ( int32_t nmarks,
    r3_t p_w[],
    r3_t b_w,
    r2_t p_u[],
    r2_t b_u,
    bool_t distance_weights,
    double weight[],
    bool_t debug,
    double P[] ) 
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }
  
  int32_t k;
  int32_t m = 2*nmarks;
  int32_t n = 12;
  int32_t p = 1;
  mat_rm_t A = tf_alloc_mat_rm (m, n);
  mat_rm_t B = tf_alloc_mat_rm (m, p);
  mat_rm_t X = tf_alloc_mat_rm (n, p);

  fprintf(stderr, "%s: System A,B:\n", __FUNCTION__);
  for (k = 0; k < nmarks; k++) {
    r3_t pwk = p_w[k];
    r2_t pik = p_i[k];
    double wtk = sqrt(weight[k]);
    int32_t i = 2*k;
    double *EA = &(A->c[i*n]);
    double *EB = &(B->c[i*p]);

    EA[0]  = wtk*puk.c[0];
    EA[1]  = wtk*puk.c[0]*pwk.c[0];
    EA[2]  = wtk*puk.c[0]*pwk.c[1];
    EA[3]  = wtk*puk.c[0]*pwk.c[2];
    EA[4]  = -wtk;
    EA[5]  = -wtk*pwk.c[0];
    EA[6]  = -wtk*pwk.c[1];
    EA[7]  = -wtk*pwk.c[2];
    EA[8]  = 0;
    EA[9]  = 0;
    EA[10] = 0;
    EA[11] = 0;                         

    i = 2*k + 1;
    double *FA = &(A->c[i*n]);
    double *FB = &(B->c[i*p]);

    FA[0]  = wtk*puk.c[1];
    FA[1]  = wtk*puk.c[1]*pwk.c[0];
    FA[2]  = wtk*puk.c[1]*pwk.c[1];
    FA[3]  = wtk*puk.c[1]*pwk.c[2];
    FA[4]  = 0;
    FA[5]  = 0;
    FA[6]  = 0;
    FA[7]  = 0;
    FA[8]  = -wtk;
    FA[9]  = -wtk*pwk.c[0];
    FA[10] = -wtk*pwk.c[1];
    FA[11] = -wtk*pwk.c[2];
 
  }

  for (i = 0; i < m; i++) {
    int32_t j;
    for (j = 0; j < n; j++) fprintf(stderr, " %12.8f ", A->c[i*n+j]);
    fprintf(stderr, " = ");  
    for (j = 0; j < p; j++) fprintf(stderr, " %12.8f ", B->c[i*p+j]);
    fprintf(stderr, "\n");  
  }
   
  /* Compute {X} such that {A X = B} in the least squares sense: */
  tf_solve_system_mxn_mxp_homogeneous (A, X, B);

  fprintf(stderr, "Solution X:\n");
  for (i = 0; i < n; i++) 
    { int32_t j;
      for (j = 0; j < p; j++) fprintf(stderr, " %12.8f ", X->c[i*p+j]);
      fprintf(stderr, "\n");
    }

  tf_free_mat_rm_structure (A);
  tf_free_mat_rm_structure (B);
  tf_free_mat_rm_structure (X);
}

void tf_calib_guess3_extract_parameters_from_P
  ( double P[], 
    bool_t affine_only,
    bool_t debug,
    r4x4_t *S,
    double *f )
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }

  /*
       P[0]  = Tz
       P[1]  = r7
       P[2]  = r8
       P[3]  = r9
       P[4]  = f*Tx
       P[5]  = f*r1
       P[6]  = f*r2
       P[7]  = f*r3
       P[8]  = f*Ty
       P[9]  = f*r4
       P[10] = f*r5
       P[11] = f*r6
  */  

  r3_t rh = (r3_t){{ P[5], P[6], P[7] }};
  r3_t sh = (r3_t){{ P[9], P[10], P[11] }};

  r3_t r, s, t; /* The three rows of {R}. */
  
  /* This needs to be improved: */
  double rabs = r3_norm(&rh);
  double sabs = r3_norm(&sh);
  (*f) = sqrt(rabs*sabs);
  if (affine_only) {
    r3_scale(1/(*f), &rh, &r);
    r3_scale(1/(*f), &sh, &s);
    t = (r3_t){{ P[1], P[2], P[3] }};
  }
  else {
    if (rabs > sabs)
      { r3_dir(&rh, r); 
        r3_decomp(&sh, r, NULL, s);
        r3_dir(s, s);      }
    else
      { r3_dir(&sh, s); 
        r3_decomp(&rh, s, NULL, r);
        r3_dir(r, r); 
      }
    r3_cross(r, s, t);
    r3_dir(t, t);
  }

  /* Set the {R} submatrix in {S}: */
  r4x4_ident(S);
  S->c[1][1] = r.c[0];
  S->c[1][2] = r.c[1];
  S->c[1][3] = r.c[2];
  S->c[2][1] = s.c[0];
  S->c[2][2] = s.c[1];
  S->c[2][3] = s.c[2];
  S->c[3][1] = t.c[0];
  S->c[3][2] = t.c[1];
  S->c[3][3] = t.c[2];
  
  /* Set {Tx,Ty,Tz}: */
  S->c[1][0] = P[4]/(*f);
  S->c[2][0] = P[8]/(*f);
  S->c[3][0] = P[0];
}
