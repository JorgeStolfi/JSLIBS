/* See {tf_calib_guess2.h}. */
/* Last edited on 2022-10-20 05:54:21 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <jsfile.h>
#include <affirm.h>
#include <tf_camera.h>
#include <tf_matrix.h>
#include <tf_errors.h>
#include <tf_math.h>
#include <tf_calib.h>
#include <tf_calib_guess2.h>

void tf_calib_guess2_initial_camera_parameters 
  ( tf_calib_data_t * cdat, 
    tf_camera_specs_t *cspec, 
    tf_camera_params_t *cpar,
    bool_t affine_only, 
    bool_t debug)
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }

  tf_optimization_choice_t which;
  tf_select_no_parameters(&which);
  tf_select_all_variable_parameters(cspec, &which);
        
  int32_t n = cdat->np;
  if (debug) {
    fprintf(stderr, "Computing the initial guess of {cpar}\n");
    tf_camera_specs_print(stderr, cspec);
    tf_show_which_parameters_are_selected(&which, stderr);
    tf_calib_data_print(stderr, cdat);
  }

  /* Some parameters must be known: */
  demand(! which.Cx, "cannot guess {Cx}");
  demand(! which.Cy, "cannot guess {Cy}");
  demand(! which.sx, "cannot guess {sx}");

  bool_t variable_R = which.R;
  bool_t variable_V = (which.v_w[0] || which.v_w[1] || which.v_w[2]);
 
  /* Compute the barycenter of world coordinates of all features: */
  r3_t *p_w = cdat->world;
  r3_t b_w; 
  tf_compute_world_barycenter(n, p_w, cdat->weight, &b_w);

  /* Compute the observed undistorted coordinates of all marks features, and their barycenter: */
  assert(cpar->kappa >= LO(cspec->kappa) && cpar->kappa <= HI(cspec->kappa)); /* Paranoia. */
  r2_t *p_u = (r2_t *)notnull(malloc(n*sizeof(r2_t)), "no mem");
  tf_compute_undistorted_obs_coordinates(n, cdat->image, cpar, p_u);
  if (debug) {
    fprintf(stderr, "Undistorted projected mark coordinates:\n");
    tf_calib_data_print_image_points (stderr, cdat->np, p_u);
  }
  r2_t b_u;
  tf_compute_undistorted_coords_barycenter(n, p_u, cdat->weight, &b_u);

  if ((! variable_R) && (! variable_V))
    {
      /* Camera position and orientation is fixed. */
      if (! which.f) { /* There is nothing left to guess! */ return; }
      /* Guess the best {f} only. */
      tf_calib_guess2_compute_f_of_fixed_camera
        (n, p_w, b_w, p_u, b_u, cdat->weight, &cpar->f, debug);
      /* Bring {f} into the given range: */
      tf_clip_to_range("f", &cpar->f, TRUE, &cspec->f);
    } 
  else if (! variable_V)
    { /* Camera has fixed position but variable orientation. */
      /* Guess best {R} and {f}, compute {Tx,Ty,Tz} from them: */
      r3_t v_w = (r3_t){{ LO(cspec->v_w[0]), LO(cspec->v_w[1]), LO(cspec->v_w[2]) }};
      tf_calib_guess2_compute_S_f_of_fixed_position_camera
	(n, p_w, b_w, p_u, b_u, cdat->weight, &v_w,  affine_only, debug, &cpar->S, &cpar->f);
      /* Bring {f} into the given range (is this the best we can do???): */
      tf_clip_to_range("f", &cpar->f, which.f, &cspec->f);
    }
  else 
    { /* Camera has variable position (and perhaps fixed orientation). */
      /* Find the best-fitting camera at very large distance from {b_w}, ignoring the {f} range: */
      tf_calib_guess2_compute_S_f_of_far_away_camera
        (n, p_w, b_w, p_u, b_u, cdat->weight, which.R, affine_only, debug, &cpar->S, &cpar->f);
    }
   
  /*disallocating structures*/
  free (p_u);
      
}

/* MAIN CASES */

void tf_calib_guess2_compute_f_of_fixed_camera
  ( int32_t n,
    r3_t p_w[],
    r3_t b_w,
    r2_t p_u[],
    r2_t b_u,
    double weight[], 
    double *f,
    bool_t debug)
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }
  
  demand(FALSE, "not implemented yet!");
}

void tf_calib_guess2_compute_S_f_of_fixed_position_camera
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
    double *f )
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }
  
  demand(FALSE, "not implemented yet!");
}

void tf_calib_guess2_compute_S_f_of_far_away_camera
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
    double *f )
{
  double mu;

  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }
  
  if (variable_R)
    { 
      tf_calib_guess2_compute_R_mu(n, p_w, b_w, p_u, b_u, weight, affine_only, debug, S, &mu);
    } 
  else 
    { 
      tf_calib_guess2_compute_mu_given_R(n, p_w, b_w, p_u, b_u, weight, debug, S, &mu);
    }
  if (debug) 
    {
      fprintf(stderr, "  Initial guess for rotation matrix {R}\n");
      tf_camera_matrix_print_rotation (S, stderr);
      fprintf(stderr, "  Initial guess of mu = %f\n", mu);
    }

  tf_calib_guess2_compute_Tx_given_R_mu(n, b_w, b_u, S, mu, debug);
  if (debug) { fprintf(stderr, "  Initial guess of Tx = %f\n", S->c[1][0]); }

  tf_calib_guess2_compute_Ty_given_R_mu(n, b_w, b_u, S, mu, debug);
  if (debug) { fprintf(stderr, "  Initial guess of Ty = %f\n", S->c[2][0]); }

  /* Choose a very large {Tz} and compute the matching {f}: */
  //S->c[3][0] = 100000.0 * r_w; /VERIFICAR/
  S->c[3][0] = 100000.0;
  if (debug) { fprintf(stderr, "  Initial guess of Tz = %f\n", S->c[3][0]); }

  /* Compute the camera Z coordinate of {b_w}: */
  double bz_c = S->c[3][0] + S->c[1][3]*b_w.c[0] + S->c[2][3]*b_w.c[1] + S->c[3][3]*b_w.c[2];
  demand(bz_c > 0.0, "point barycenter is behind faraway camera!");
  *f = mu*bz_c;
  if (debug) { fprintf(stderr, "  Initial guess of f = %f, mu = %f, bz_c = %f\n", *f, mu, bz_c); }

}

void tf_calib_guess2_compute_R_mu
  ( int32_t n,
    r3_t p_w[],
    r3_t b_w,
    r2_t p_u[],
    r2_t b_u,
    double weight[],
    bool_t affine_only,
    bool_t debug,
    r4x4_t *S,
    double *mu )
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }
  
  /* Find the best-fitting affine mapping from {p_w} to {p_u}: */
  r3_t rh, sh;
  r2_t d;
  tf_calib_guess2_compute_initial_affine_model(n, p_w, b_w, p_u, b_u, weight, debug, &rh, &sh, &d);

  /* Find the orthogonal unit vectors {r,s} that best fit {rh,sh}: */
  r3_t r, s, t;
  tf_calib_guess2_extract_camera_vectors(rh, sh, affine_only, debug, &r, &s, &t, mu);

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
}

void tf_calib_guess2_compute_mu_given_R
  ( int32_t n,
    r3_t p_w[],
    r3_t b_w,
    r2_t p_u[],
    r2_t b_u,
    double weight[],
    bool_t debug,
    r4x4_t *S,
    double *mu )
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }
  
  demand(FALSE, "not implemented yet!");
}

void tf_calib_guess2_compute_initial_affine_model
  ( int32_t n,
    r3_t p_w[],
    r3_t b_w,
    r2_t p_u[],
    r2_t b_u,
    double weight[],
    bool_t debug,
    r3_t *rh,
    r3_t *sh,
    r2_t *d ) 
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }
  
  int32_t i;
  mat_rm_t A = tf_alloc_mat_rm (n, 3);
  mat_rm_t B = tf_alloc_mat_rm (n, 2);
  mat_rm_t X = tf_alloc_mat_rm (3, 2);

  fprintf(stderr, "%s: System A,B:\n", __FUNCTION__);
  for (i = 0; i < n; i++) {
    double *Ai = &(A->c[i*3]);
    Ai[0] = p_w[i].c[0] - b_w.c[0];
    Ai[1] = p_w[i].c[1] - b_w.c[1];
    Ai[2] = p_w[i].c[2] - b_w.c[2];

    double *Bi = &(B->c[i*2]);
    Bi[0] = p_u[i].c[0] - b_u.c[0];
    Bi[1] = p_u[i].c[1] - b_u.c[1];
  }

  for (i = 0; i < n; i++) {
    int32_t j;
    for (j = 0; j < 3; j++) fprintf(stderr, " %12.8f ", A->c[i*3+j]);
    fprintf(stderr, " = ");  
    for (j = 0; j < 2; j++) fprintf(stderr, " %12.8f ", B->c[i*2+j]);
    fprintf(stderr, "\n");  
  }
   
  /* Compute {X} such that {A X = B} in the least squares sense: */
  tf_solve_system_mxn_mxp (A, X, B, weight);

  fprintf(stderr, "Solution X:\n");
  for (i = 0; i < 3; i++) 
    { int32_t j;
      for (j = 0; j < 2; j++) fprintf(stderr, " %12.8f ", X->c[i*2+j]);
      fprintf(stderr, "\n");
    }

  /* Extract {rh,sh} and compute the constant term {d}: */
  *d = b_u;
  for (i = 0; i < 3; i++)
    { int32_t j;
      for (j = 0; j < 2; j++) { d->c[j] -= (b_w.c[i]*X->c[i*2+j]); }
      rh->c[i] = X->c[i*2+0];
      sh->c[i] = X->c[i*2+1];
    }

  tf_free_mat_rm_structure (A);
  tf_free_mat_rm_structure (B);
  tf_free_mat_rm_structure (X);
}

void tf_calib_guess2_extract_camera_vectors
  ( r3_t rh, 
    r3_t sh, 
    bool_t affine_only,
    bool_t debug,
    r3_t *r, 
    r3_t *s, 
    r3_t *t,
    double *mu )
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }
  
  /* This needs to be improved: */
  double rabs = r3_norm(&rh);
  double sabs = r3_norm(&sh);
  (*mu) = sqrt(rabs*sabs);
  if (affine_only) {
    r3_scale(1/(*mu), &rh, r);
    r3_scale(1/(*mu), &sh, s);
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
  }
  r3_cross(r, s, t);
  r3_dir(t, t);
}

void tf_calib_guess2_compute_Tx_given_R_mu (int32_t n, r3_t b_w, r2_t b_u, r4x4_t *S, double mu, bool_t debug)
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }
  
  r3_t r = {{S->c[1][1], S->c[1][2], S->c[1][3]}};
  double d = r3_dot (&b_w, &r);
  fprintf(stderr, "Tx: Xw  = %f %f %f\n", r.c[0], r.c[1], r.c[2]);
  fprintf(stderr, "Tx: pw = %f %f %f\n", b_w.c[0], b_w.c[1], b_w.c[2]);
  fprintf(stderr, "Tx: d  = %f\n", d);
  fprintf(stderr, "Tx: p_u = %f %f\n", b_u.c[0], b_u.c[1]);
  fprintf(stderr, "mu: %f\n", mu);
  S->c[1][0] = b_u.c[0]/mu - d;
}

void tf_calib_guess2_compute_Ty_given_R_mu (int32_t n, r3_t b_w, r2_t b_u, r4x4_t *S, double mu, bool_t debug)
{
  if (debug) { fprintf(stderr, "Entering %s\n", __FUNCTION__); }
  
  r3_t s = {{S->c[2][1], S->c[2][2], S->c[2][3]}};
  double d = r3_dot (&b_w, &s);
  S->c[2][0] = b_u.c[1]/mu - d;
}
