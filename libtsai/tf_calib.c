/* See {tf_calib.h}. */
/* Last edited on 2023-02-25 16:11:40 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>

#include <jsfile.h>
#include <affirm.h>
#include <rn.h>
#include <rmxn.h>

#include <tf_lmdif.h>

#include <tf_camera.h>
#include <tf_matrix.h>
#include <tf_math.h>
#include <tf_errors.h>

#include <tf_calib.h>

/* Parameters controlling MINPACK's lmdif() optimization routine. */
/* See the file lmdif.f for definitions of each parameter.        */
#define REL_SENSOR_TOLERANCE_ftol    1.0E-5      /* [pix] */
#define REL_PARAM_TOLERANCE_xtol     1.0E-7
#define ORTHO_TOLERANCE_gtol         0.0
#define LMDIF_EPSFCN                 1.0E-12     /* Do not set to zero! */
#define LMDIF_MODE                   1           /* variables are scalled internally */
#define LMDIF_FACTOR                 100.0 

void tf_calib_generic
  ( tf_calib_data_t * cdat, 
    tf_camera_specs_t *cspec,
    r2_t p_i_dev, 
    bool_t use_cpar_dev_terms,
    tf_camera_params_t *cpar,
    tf_gather_params_proc_t gather_params,
    tf_scatter_params_proc_t scatter_params,
    bool_t debug
  )
  {
    /*Decide which parameters are variable*/
    tf_optimization_choice_t which;
    tf_select_no_parameters(&which);
    tf_select_all_variable_parameters(cspec, &which);
    if (debug) 
      { fprintf(stderr, "Guessed parameters:\n");
        tf_camera_params_print (cpar, stderr);
      }
    tf_generic_optimization 
      ( cdat, cspec, p_i_dev, use_cpar_dev_terms, cpar, &which, gather_params, scatter_params, debug );
    if (debug) 
      { fprintf(stderr, "Optimized parameters:\n");
    	tf_camera_params_print (cpar, stderr);
      }
  } 

void tf_compute_world_barycenter (int32_t n, r3_t p_w[], double weight[], r3_t *b_w)
{
  *b_w = (r3_t){{0,0,0}};
  double sum_w = 0.0;
  int32_t i;
    for (i = 0; i < n; i++) {
        b_w->c[0] += p_w[i].c[0]*weight[i];
	b_w->c[1] += p_w[i].c[1]*weight[i];
	b_w->c[2] += p_w[i].c[2]*weight[i];
	sum_w += weight[i]; 
    }

    r3_scale(1/sum_w, b_w, b_w);
}

void tf_compute_undistorted_coords_barycenter(int32_t n, r2_t p_u[], double weight[], r2_t *b_u)
{
  *b_u = (r2_t){{0,0}};
  double sum_w = 0.0;
  int32_t i;
    for (i = 0; i < n; i++) {
      b_u->c[0] += p_u[i].c[0]*weight[i];
	b_u->c[1] += p_u[i].c[1]*weight[i];
	sum_w += weight[i]; 
    }

    r2_scale(1/sum_w, b_u, b_u);
}

void tf_compute_undistorted_obs_coordinates(int32_t n, r2_t p_i[], tf_camera_params_t *cpar, r2_t p_u[])
{
  int32_t i;

  for (i = 0; i < n; i++) {
    r2_t p_d = tf_image_coords_to_sensor_coords (cpar, p_i[i]);
    p_u[i] = tf_dis_sensor_coords_to_und_sensor_coords (cpar, p_d); 
  }
}

void tf_generic_optimization_compute_error_terms
  ( double *err,
    int32_t nerr,
    tf_calib_data_t * cdat,
    tf_camera_specs_t *cspec,
    r2_t p_i_dev,
    bool_t use_cpar_dev_terms,
    tf_camera_params_t *cpar,
    tf_optimization_choice_t *which )
{
    int32_t i;

    int32_t kerr = 0;
    
    /* Mark position terms: */

    for (i = 0; i < cdat->np; i++) {

      /* Map the known world coordinates to (predicted) image coordinates using {cpar}: */
      r2_t e_i = tf_camera_compute_image_error(cpar, cdat->world[i], cdat->image[i]);

      /* Terms due to errors in the undistorted sensor coordinates */
      double wti = sqrt(cdat->weight[i]);
      err[kerr] = wti*e_i.c[0]/p_i_dev.c[0];
      kerr++;
      err[kerr] = wti*e_i.c[1]/p_i_dev.c[1];
      kerr++;
    }

    if (use_cpar_dev_terms) {
      /* Parameter strangeness terms: */

      if (which->sx) {
	assert(! interval_is_trivial (&(cspec->sx)) );
	assert(interval_is_finite  (&(cspec->sx)) );
	double avg = interval_mid(&cspec->sx);
	double dev = interval_rad(&cspec->sx)/3;
	double val = cpar->sx;
	err[kerr] = (val - avg)/dev;
	kerr++;
      }

      if (which->Cx) {
	assert(! interval_is_trivial (&(cspec->Cx)) );
	assert(interval_is_finite  (&(cspec->Cx)) );
	double avg = interval_mid(&cspec->Cx);
	double dev = interval_rad(&cspec->Cx)/3;
	double val = cpar->Cx;
	err[kerr] = (val - avg)/dev;
	kerr++;
      }

      if (which->Cy) {
	assert(! interval_is_trivial (&(cspec->Cy)) );
	assert(interval_is_finite  (&(cspec->Cy)) );
	double avg = interval_mid(&cspec->Cy);
	double dev = interval_rad(&cspec->Cy)/3;
	double val = cpar->Cy;
	err[kerr] = (val - avg)/dev;
	kerr++;
      }

      if (which->kappa) {
	assert(! interval_is_trivial (&(cspec->kappa)) );
	assert(interval_is_finite  (&(cspec->kappa)) );
	double avg = interval_mid(&cspec->kappa);
	double dev = interval_rad(&cspec->kappa)/3;
	double val = cpar->kappa;
	err[kerr] = (val - avg)/dev;
	kerr++;
      }
    
      if (which->f) {
	assert(! interval_is_trivial (&(cspec->f)) );
	assert(interval_is_finite  (&(cspec->f)) );
	/* Work in log scale: */
	double logf = tf_camera_safe_log(cpar->f);
	interval_t logf_range = tf_camera_interval_safe_log(&cspec->f);
	double avg = interval_mid(&logf_range);
	double dev = interval_rad(&logf_range)/3;
	double val = logf;
	err[kerr] = (val - avg)/dev;
	kerr++;
      }

      if (which->v_w[0] || which->v_w[1] || which->v_w[2]) {
	r3_t v_w = tf_camera_matrix_to_v_w (&(cpar->S));
	int32_t j;
	for (j = 0; j < 3; j++) {
	  if (which->v_w[j]) {
	    assert(! interval_is_trivial (&(cspec->v_w[j])) );
	    if (! interval_is_finite (&(cspec->v_w[j])) ) {
	      assert(interval_is_full (&(cspec->v_w[j])) );
	      err[kerr] = 0;
	    } else {
	      double avg = interval_mid(&(cspec->v_w[j]));
	      double dev = interval_rad(&(cspec->v_w[j]))/3;
	      double val = v_w.c[j];
	      err[kerr] = (val - avg)/dev;
	    }
	    kerr++;
	  }
	}
      }

      if (which->R) {
	int32_t j;
	for (j = 0; j < 3; j++) {
          assert(! interval_is_trivial (&(cspec->R[j])) );
	}

	if( (! interval_is_finite(&(cspec->R[0]))) || 
	    (! interval_is_finite(&(cspec->R[1]))) || 
	    (! interval_is_finite(&(cspec->R[2]))) ) {
	  /* If any is infinite, all must be full: */
	  for (j = 0; j < 3; j++) { 
	    assert( interval_is_full (&(cspec->R[j])) ); 
	    err[kerr] = 0; kerr++;
	  }
	} else {
	  /* All are finite: */
	  /* !!! Pensar melhor !!! */
	  r3_t R;
	  tf_camera_matrix_to_euler_angles (&(cpar->S), &R);
	  int32_t j;
	  for (j = 0; j < 3; j++) {
	    double avg = interval_mid(&(cspec->R[j]));
	    double dev = interval_rad(&(cspec->R[j]))/3;
	    double val = R.c[j];
	    err[kerr] = sin(val-avg)/dev; kerr++;
	  }
	} 
      }

    }

    assert(kerr == nerr);

}

void tf_generic_optimization_print_error_terms
  ( FILE *ferr,
    double *err,
    int32_t nerr,
    bool_t use_cpar_dev_terms,
    tf_calib_data_t * cdat,
    tf_optimization_choice_t *which )
{
    int32_t i;

    int32_t kerr = 0;
    
    /* Mark position terms: */

    for (i = 0; i < cdat->np; i++) {
      double wt_error_X = err[kerr];
      kerr++;
      double wt_error_Y = err[kerr];
      kerr++;
      fprintf(ferr, "  mark[%3d] X = %+15.8f  Y = %+15.8f\n", i, wt_error_X, wt_error_Y);
    }

    if (use_cpar_dev_terms) {
      /* Parameter strangeness terms: */

      if (which->sx) {
	double wt_error_sx = err[kerr];
	kerr++;
	fprintf(ferr, "  sx = %+15.8f\n", wt_error_sx);
      }

      if (which->Cx) {
	double wt_error_Cx = err[kerr];
	kerr++;
	fprintf(ferr, "  Cx = %+15.8f\n", wt_error_Cx);
      }

      if (which->Cy) {
	double wt_error_Cy = err[kerr];
	kerr++;
	fprintf(ferr, "  Cy = %+15.8f\n", wt_error_Cy);
      }

      if (which->kappa) {
	double wt_error_kappa = err[kerr];
	kerr++;
	fprintf(ferr, "  kappa = %+15.8f\n", wt_error_kappa);
      }
    
      if (which->f) {
	double wt_error_logf = err[kerr];
	kerr++;
	fprintf(ferr, "  logf = %+15.8f\n", wt_error_logf);
      }

      if (which->v_w[0] || which->v_w[1] || which->v_w[2]) {
	r3_t wt_error_v_w;
	int32_t j;
	for (j = 0; j < 3; j++) {
	  wt_error_v_w.c[j] = err[kerr];
	  kerr++;
	}
	fprintf(ferr, "  v_w: X = %+15.8f  Y = %+15.8f  Z = %+15.8f\n", wt_error_v_w.c[0], wt_error_v_w.c[1], wt_error_v_w.c[2]);
      
      }

      if (which->R) {
	r3_t wt_error_R;
	int32_t j;
	for (j = 0; j < 3; j++) {
	  wt_error_R.c[j] = err[kerr];
	  kerr++;
	} 
	fprintf(ferr, "  R:  X = %+15.8f  Y = %+15.8f  Z = %+15.8f\n", wt_error_R.c[0], wt_error_R.c[1], wt_error_R.c[2]);
      }
    }

    assert(kerr == nerr);
}

#define MAX_OPT_PARAMS 11

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
  ) 
{
  if (debug) {
    fprintf(stderr, "--- tf_generic_optimization (enter) ----------\n");
    tf_show_values_of_selected_parameters (cpar, which, stderr);
    fprintf(stderr, "--------------------------------------------\n");
  }

  /* no use optimizing parameters which are fixed: */
  tf_unselect_all_fixed_parameters(cspec, which);

  // demand(! which->Npx, "cannot calibrate Npx");
  // demand(! which->Npy, "cannot calibrate Npy");
  // demand(! which->dpx, "cannot calibrate dpx");
  // demand(! which->dpy, "cannot calibrate dpy");

  int32_t nparams = tf_generic_optimization_count_params(which); 
  int32_t nerrs   = tf_generic_optimization_count_error_terms(cdat, use_cpar_dev_terms,    which); 

  if (debug) {
    fprintf(stderr, "NPARAMS %d\n", nparams);
  }

  /* parameters needed by MINPACK's lmdif() */

  int32_t     m = nerrs;   /* num error terms in the goal function. */
  int32_t     n = nparams; /* num parameters. */
 
  
  double  ftol = REL_SENSOR_TOLERANCE_ftol;
  double  xtol = REL_PARAM_TOLERANCE_xtol;
  double  gtol = ORTHO_TOLERANCE_gtol;
  int32_t     maxfev = 1000*n;
  double  epsfcn = LMDIF_EPSFCN;

  int32_t     mode = LMDIF_MODE;
  double  factor = LMDIF_FACTOR;
  int32_t     nprint = 0;
  int32_t     info;
  int32_t     nfev;

  int32_t     ldfjac = m;
  double  x[MAX_OPT_PARAMS];
  double  diag[MAX_OPT_PARAMS];
  int32_t     ipvt[MAX_OPT_PARAMS];
  double  qtf[MAX_OPT_PARAMS];
  double  wa1[MAX_OPT_PARAMS];
  double  wa2[MAX_OPT_PARAMS];
  double  wa3[MAX_OPT_PARAMS];

  double *fvec;
  double *fjac;
  double *wa4;  

  fvec = rn_alloc(m);
  fjac = rmxn_alloc(m,n);
  wa4 =  rn_alloc(m); 


  if (debug) 
  { 
      fprintf(stderr, "Camera parameters after clipping:\n");
      tf_camera_params_print (cpar, stderr);
  }


  /* use the current calibration and camera constants as a starting point */
  gather_params(nparams, x, cspec, cdat, cpar, which);

  /* define optional scale factors for the parameters */
  int32_t i;

  if (mode == 2) {
    for (i = 0; i < nparams; i++)
      diag[i] = 1.0; /* some user-defined values */
  }

  auto void error (int32_t ns, int32_t np, double x[], double fvec[], int32_t *iflag);
  /* goal function for MINPACK's lmdif. */
    
  void error (int32_t ne, int32_t np, double x[], double fvec[], int32_t *iflag) {
    assert(ne == nerrs); assert(np == nparams);
    scatter_params(np, x, cspec, cdat, cpar, which); 
    if (debug) {
      fprintf(stderr, "--- tf_generic_optimization.error ----------\n");
      tf_show_values_of_selected_parameters (cpar, which, stderr);
    }
    tf_generic_optimization_compute_error_terms 
      ( fvec, ne, cdat, cspec, p_i_dev, use_cpar_dev_terms, cpar, which );
    if (debug) {
      tf_show_optimization_errors(fvec, cdat->weight, nerrs, stderr);
      fprintf(stderr, "--------------------------------------------\n");
    }
  }
    
  /* perform the optimization */
  lmdif (error,  m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,
	 diag, mode, factor, nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4);

  /* update the calibration and camera constants */
  scatter_params(nparams, x, cspec, cdat, cpar, which);

  if (debug) {
    fprintf(stderr, "--- tf_generic_optimization (exit raw) ----------\n");
    tf_show_values_of_selected_parameters(cpar, which, stderr);
    //tf_show_optimization_errors(fvec, cdat->weight, m, stderr);
    fprintf(stderr, "--------------------------------------------\n");
  }

  /* release allocated workspace */
  free (fvec);
  free (fjac);
  free (wa4);

}

int32_t tf_generic_optimization_count_params (tf_optimization_choice_t *which)
{
  /* count the variable parameters: */
  bool_t which_Tx_Ty = which->v_w[0] | which->v_w[1] | which->v_w[2] | which->R;
  bool_t which_L = which->v_w[0] | which->v_w[1] | which->v_w[2] | which->R | which->f;

  int32_t nparams =  (which->R*3) + (which_Tx_Ty*2) + (which_L) + (which->f)
                 + (which->kappa) + (which->sx) +  (which->Cx) + (which->Cy); 
  return nparams;
}

int32_t tf_generic_optimization_count_error_terms 
  ( tf_calib_data_t * cdat,
    bool_t use_cpar_dev_terms,
    tf_optimization_choice_t *which )
{
  int32_t nerrs = 2*cdat->np;
  if (use_cpar_dev_terms) {
    nerrs += (which->R*3) + (which->v_w[0]) + (which->v_w[1]) + (which->v_w[2])
              + (which->f) + (which->kappa) + (which->sx) +  (which->Cx) + (which->Cy); 
  }
  return nerrs;
}

void tf_clip_to_range (char *name, double *v, bool_t variable, interval_t *range)
{

  if (variable)    
    { double v_old = *v;
      *v = interval_project(range, *v);
      if (*v != v_old) {
         fprintf(stderr, "%s adjusted from %f ", name, *v);
         fprintf(stderr, "to %f\n", *v);
      }
    }  
  else
    { if (((*v) < LO(*range)) || ((*v) > HI(*range)))
        { fprintf(stderr, "%s = %f", name, (*v));
          fprintf(stderr, " is out of the range [ %f _ %f ]!\n", LO(*range), HI(*range));
          assert(FALSE);
	}
    }
}

void tf_show_optimization_errors
( double err[],
  double weight[],
  int32_t m,
  FILE *ferr )
{
  int32_t nmarks = m/2;
  assert(m == (2*nmarks));
  int32_t i;
  double sum = 0.0;
  fprintf(ferr, "----- tf_show_optimization_errors -----\n");
  for (i = 0; i < nmarks; i++) {
    fprintf(ferr, "err[%02d] X = %10.3f, Y = %10.3f, weight = %15.12f\n", i, err[2*i], err[2*i+1], weight[i]);
    sum += err[2*i]*err[2*i] + err[2*i+1]*err[2*i+1];
  }
  fprintf(ferr, "RMS error  = %10.6f\n", sqrt(sum/m));
  fprintf(ferr, "----------------------------\n");
}

void tf_show_values_of_selected_parameters
( tf_camera_params_t *cpar,
  tf_optimization_choice_t *which, 
  FILE *ferr )
{
    
  // if (which->Npx) {
  //   fprintf(ferr, "  .%-8s", "Npx"); 
  //   if (cpar != NULL) { fprintf(ferr, "= %24.16e", cpar->Npx); }
  //   fprintf(ferr, "\n");
  // }

  // if (which->Npy) {
  // fprintf(ferr, "  .%-8s", "Npy"); 
  // if (cpar != NULL) { fprintf(ferr, "= %24.16e", cpar->Npy); }
  // fprintf(ferr, "\n");
  // }

  // if (which->dpx) {
  // fprintf(ferr, "  .%-8s", "dpx"); 
  // if (cpar != NULL) { fprintf(ferr, "= %24.16e", cpar->dpx); }
  // fprintf(ferr, "\n");
  // }

  // if (which->dpy) {
  // fprintf(ferr, "  .%-8s", "dpy"); 
  // if (cpar != NULL) { fprintf(ferr, "= %24.16e", cpar->dpy); }
  // fprintf(ferr, "\n");
  // }

  if (which->f) {
    fprintf(ferr, "  .%-8s", "f"); 
    if (cpar != NULL) { fprintf(ferr, "= %24.16e", cpar->f); }
    fprintf(ferr, "\n");
  }

  if (which->kappa) {
    fprintf(ferr, "  .%-8s", "kappa"); 
    if (cpar != NULL) { fprintf(ferr, "= %24.16e", cpar->kappa); }
    fprintf(ferr, "\n");
  }

  if (which->sx) {
    fprintf(ferr, "  .%-8s", "sx"); 
    if (cpar != NULL) { fprintf(ferr, "= %24.16e", cpar->sx); }
    fprintf(ferr, "\n");
  }

  if (which->Cx) {
    fprintf(ferr, "  .%-8s", "Cx"); 
    if (cpar != NULL) { fprintf(ferr, "= %24.16e", cpar->Cx); }
    fprintf(ferr, "\n");
  }

  if (which->Cy) {
    fprintf(ferr, "  .%-8s", "Cy"); 
    if (cpar != NULL) { fprintf(ferr, "= %24.16e", cpar->Cy); }
    fprintf(ferr, "\n");
  }

  if (which->R) {
    fprintf(ferr, "  .%-8s", "R"); 
    if (cpar != NULL) {
      r3_t R;
      tf_camera_matrix_to_euler_angles(&(cpar->S), &R);
      r3_gen_print(ferr, &R, "%24.16e", "= ( ", " ", " )");
    }
    fprintf(ferr, "\n");
  }

  int32_t i;
  for (i = 0; i < 3; i++) {
    if (which->v_w[i]) {
      fprintf(ferr, "  .%-8s[%d] =", "v_w", i); 
      if (cpar != NULL) { fprintf(ferr, " %24.16e", cpar->S.c[i+1][0]); }
      fprintf(ferr, "\n");
    }
  }
}

void tf_select_no_parameters(tf_optimization_choice_t *which)
{ 
  (*which) = (tf_optimization_choice_t){
    // .Npx = FALSE,
    // .Npy = FALSE,
    // .dpx = FALSE,
    // .dpy = FALSE,
    .R =      FALSE,
    .v_w[0] = FALSE,
    .v_w[1] = FALSE,
    .v_w[2] = FALSE,
    .f =      FALSE,
    .kappa =  FALSE,
    .sx =     FALSE,
    .Cx =     FALSE,
    .Cy =     FALSE
  };
}

void tf_select_all_parameters(tf_optimization_choice_t *which)
{ 
  (*which) = (tf_optimization_choice_t){
    // .Npx = TRUE,
    // .Npy = TRUE,
    // .dpx = TRUE,
    // .dpy = TRUE,
    .R =      TRUE,
    .v_w[0] = TRUE,
    .v_w[1] = TRUE,
    .v_w[2] = TRUE,
    .f =      TRUE,
    .kappa =  TRUE,
    .sx =     TRUE,
    .Cx =     TRUE,
    .Cy =     TRUE
  };
}

void tf_select_all_variable_parameters(tf_camera_specs_t *cspec, tf_optimization_choice_t *which)
{
  demand(LO(cspec->Npx) ==  HI(cspec->Npx), "cannot set {which->Npx}");
  demand(LO(cspec->Npy) ==  HI(cspec->Npy), "cannot set {which->Npy}");
  demand(LO(cspec->dpx) ==  HI(cspec->dpx), "cannot set {which->dpx}");
  demand(LO(cspec->dpy) ==  HI(cspec->dpy), "cannot set {which->dpy}");
  if 
    ( (LO(cspec->R[0]) < HI(cspec->R[0])) || 
      (LO(cspec->R[1]) < HI(cspec->R[1])) || 
      (LO(cspec->R[2]) < HI(cspec->R[2]))
    ) { which->R = TRUE; }
  int32_t i;
  for (i = 0; i < 3; i++) {
     if (LO(cspec->v_w[i]) < HI(cspec->v_w[i])) { which->v_w[i] = TRUE; }
  }
  if (LO(cspec->f)      < HI(cspec->f))      { which->f =     TRUE; }
  if (LO(cspec->kappa)  < HI(cspec->kappa))  { which->kappa = TRUE; }
  if (LO(cspec->sx)     < HI(cspec->sx))     { which->sx =    TRUE; }
  if (LO(cspec->Cx)     < HI(cspec->Cx))     { which->Cx =    TRUE; }
  if (LO(cspec->Cy)     < HI(cspec->Cy))     { which->Cy =    TRUE; }
}

void tf_unselect_all_fixed_parameters(tf_camera_specs_t *cspec, tf_optimization_choice_t *which)
{
  demand(LO(cspec->Npx) ==  HI(cspec->Npx), "cannot set {which->Npx}");
  demand(LO(cspec->Npy) ==  HI(cspec->Npy), "cannot set {which->Npy}");
  demand(LO(cspec->dpx) ==  HI(cspec->dpx), "cannot set {which->dpx}");
  demand(LO(cspec->dpy) ==  HI(cspec->dpy), "cannot set {which->dpy}");
  if 
    ( (LO(cspec->R[0]) >= HI(cspec->R[0])) && 
      (LO(cspec->R[1]) >= HI(cspec->R[1])) && 
      (LO(cspec->R[2]) >= HI(cspec->R[2]))
    ) { which->R = FALSE; }
  int32_t i;
  for (i = 0; i < 3; i++) {
     if (LO(cspec->v_w[i])  >= HI(cspec->v_w[i]))  { which->v_w[i] = FALSE; }
  }

  if (LO(cspec->f)     >= HI(cspec->f))     { which->f =     FALSE; }
  if (LO(cspec->kappa) >= HI(cspec->kappa)) { which->kappa = FALSE; }
  if (LO(cspec->sx)    >= HI(cspec->sx))    { which->sx =    FALSE; }
  if (LO(cspec->Cx)    >= HI(cspec->Cx))    { which->Cx =    FALSE; }
  if (LO(cspec->Cy)    >= HI(cspec->Cy))    { which->Cy =    FALSE; }
}

void tf_show_which_parameters_are_selected(tf_optimization_choice_t *which, FILE *ferr)
{
  fprintf(ferr, "Selected parameters = (");
  // if (which->Npx)    { fprintf(ferr, " Npx"); }
  // if (which->Npy)    { fprintf(ferr, " Npy"); }
  // if (which->dpx)    { fprintf(ferr, " dpx"); }
  // if (which->dpy)    { fprintf(ferr, " dpy"); }
  if (which->R)      { fprintf(ferr, " R"); }
  if (which->v_w[0]) { fprintf(ferr, " v_w[0]"); }
  if (which->v_w[1]) { fprintf(ferr, " v_w[1]"); }
  if (which->v_w[2]) { fprintf(ferr, " v_w[2]"); }
  if (which->f)      { fprintf(ferr, " f"); }
  if (which->kappa)  { fprintf(ferr, " kappa"); }
  if (which->sx)     { fprintf(ferr, " sx"); }
  if (which->Cx)     { fprintf(ferr, " Cx"); }
  if (which->Cy)     { fprintf(ferr, " Cy"); }
  fprintf(ferr, " )\n");
}

void tf_write_cpar_and_errors
  ( char *name, 
    tf_calib_data_t * cdat,
    tf_camera_specs_t *cspec,
    r2_t p_i_dev, 
    bool_t use_cpar_dev_terms,
    tf_camera_params_t *cpar,
    char *out_dir )
{
  if (cpar != NULL) {

    tf_optimization_choice_t which;
    tf_select_no_parameters(&which);
    tf_select_all_variable_parameters(cspec, &which);
        
    /* Output file with given parameters: */
    char *cpar_fname = NULL;
    asprintf(&cpar_fname, "%s/%s.cpar", out_dir, name);
    FILE *cpar_file = open_write(cpar_fname, TRUE);
    free(cpar_fname);
    tf_camera_params_write(cpar_file, 0, cpar);
    fclose(cpar_file);

    /* Output file with errors for given parameters: */
    char *ferr_fname = NULL;
    asprintf(&ferr_fname, "%s/%s.ferr", out_dir, name);
    FILE *ferr = open_write(ferr_fname, TRUE);
    free(ferr_fname);

    int32_t nerr = tf_generic_optimization_count_error_terms (cdat, use_cpar_dev_terms, &which);
    double err[nerr];

    tf_generic_optimization_compute_error_terms
    ( err,
      nerr,
      cdat,
      cspec, 
      p_i_dev,
      use_cpar_dev_terms,    
      cpar,
      &which );

    tf_generic_optimization_print_error_terms 
    ( ferr, 
      err, 
      nerr,
      use_cpar_dev_terms,     
      cdat,
      &which );

    tf_calib_summarize_errors(cdat, cpar, ferr);
    fclose(ferr); 

    fprintf(stderr, "\n\n+++ %s_cpar +++++++++++++++++++++++\n", name);
    tf_camera_params_print (cpar, stderr);
  }
}

void tf_recompute_target_weights
  ( tf_camera_params_t *cpar,
    int32_t ntargets,
    r3_t p_w[],
    r2_t p_i[],
    double p_wgt[],
    r2_t dev_gud,
    r2_t dev_bad,
    bool_t debug )
{
  if (debug) { fprintf(stderr, "dev_gud = ( %8.4lf %8.4lf )\n", dev_gud.c[0], dev_gud.c[1]); }
  if (debug) { fprintf(stderr, "dev_bad = ( %8.4lf %8.4lf )\n", dev_bad.c[0], dev_bad.c[1]); }
  int32_t i;
  for (i = 0; i < ntargets; i++) {
    if (debug) { fprintf(stderr, "  mark %03d:\n", i); }
    r2_t e = tf_camera_compute_image_error(cpar, p_w[i], p_i[i]);
    if (debug) {
      fprintf(stderr, "    e = ( %8.4lf %8.4lf ) old weight = %10.8f\n", e.c[0], e.c[1], p_wgt[i]);
    }
    demand(p_wgt[i] <= 1.0, "confidence weight must not exceed 1");
    demand(p_wgt[i] >= 0.0, "confidence weight must be non-negative");
    
    double log_Pr_e_gud = 0.0; /* Probability of error {e} assuming it is an INLIER. */
    double log_Pr_e_bad = 0.0; /* Probability of error {e} assuming it is an OUTLIER. */
    int32_t axis;
    for (axis = 0; axis < 2; axis++) {
      /* Inlier's gaussian: */
      double sigma_gud = dev_gud.c[axis];
      double z_gud = e.c[axis]/sigma_gud;
      log_Pr_e_gud += -(z_gud*z_gud)/2 - log(sigma_gud*sqrt(2*M_PI));
      /* Outlier's gaussian: */
      double sigma_bad = dev_bad.c[axis];
      double z_bad = e.c[axis]/sigma_bad;
      log_Pr_e_bad += -(z_bad*z_bad)/2 - log(sigma_bad*sqrt(2*M_PI));
      if (debug) {
        fprintf(stderr, "    axis %d: z_gud = %8.4lf  z_bad = %8.4lf\n", axis, z_gud, z_bad);
      }
    }
    /* Subtract the largest from both logs to avoid underflow: */
    double C = fmax(log_Pr_e_gud, log_Pr_e_bad);
    double Pr_e_gud = exp(log_Pr_e_gud - C); 
    double Pr_e_bad = exp(log_Pr_e_bad - C);
    /* Multiply by the /a priori/ terms: */
    Pr_e_gud *= p_wgt[i];
    Pr_e_bad *= (1.0 - p_wgt[i]);
    if (debug) { 
      fprintf(stderr, "    Pr(e|gud) = %10.4e  Pr(e|bad) = %10.4e\n", Pr_e_gud, Pr_e_bad);
    }
    /* Bayes's formula: */
    double Pr_gud = Pr_e_gud/(Pr_e_gud + Pr_e_bad);
    p_wgt[i] = Pr_gud;
    if (debug) { fprintf(stderr, "    new weight = %10.8f\n\n", p_wgt[i]); }
  }
}
