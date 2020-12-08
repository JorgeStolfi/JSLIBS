/* Last edited on 2011-05-17 02:57:01 by stolfi */

#define PROG_NAME "test_calib_opt"
#define PROG_DESC "tests the camera calibration optimization algorithms"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include <float_image.h>
#include <affirm.h>
#include <r2x2.h>
#include <r2.h>
#include <pswr.h>
#include <jsfile.h>

#include <tf_calib_refine2.h>
#include <tf_calib_guess2.h>
#include <tf_camera.h>
#include <tf_camera_plot.h>
#include <tf_errors.h>
#include <tf_targets.h>
#include <tf_calib.h>
#include <sys/stat.h>
#include <sys/types.h>

void plot_cpars
  ( tf_calib_data_t *cdat,
    tf_camera_params_t *guess_cpar,
    tf_camera_params_t *opt1_cpar,
    tf_camera_params_t *opt2_cpar,
    tf_camera_params_t *true_cpar,
    char *out_dir
  );

tf_calib_data_t *read_calibration_data(char *in_dir, char *tag);

tf_camera_params_t *read_camera_parameters(tf_camera_specs_t *cspec, char *in_dir, char *name);

int main (int argc, char *argv[])
{
  char *camera_name = argv[1];
  char *in_dir = argv[2];  /* Directory where to find {p_i.txt,p_w.txt,p_wgt.txt,true.cpar}. */
  char *out_dir = argv[3]; /* Directory for all output files. */

  /* Shall we include the {cpar} Kalman error terms in the optimization? */
  bool_t use_dev = FALSE;

  fprintf(stderr, "Starting:\n");
  
  /* Get the fixed camera parameters: */
  tf_camera_specs_t *cspec;
  if (strcmp(camera_name,"optura") == 0) {
    cspec = tf_camera_specs_for_canon_optura ();
  }
  else if (strcmp(camera_name,"svga") == 0) {
    cspec = tf_camera_specs_for_povray_svga ();
  }
  else if (strcmp(camera_name,"hvga") == 0) {
    cspec = tf_camera_specs_for_povray_hvga ();
  }
  else if (strcmp(camera_name,"sony") == 0) {
    cspec = tf_camera_specs_for_sony_dv40 ();
  }
  else {
    fprintf(stderr, "error: define the camera specifications\n");
    exit(1);
  }
  fprintf(stderr, "Generic camera parameters:\n");
  tf_camera_specs_print(stderr, cspec);

  /* A priori uncertainty on mark image positions (for error computation): */
  r2_t q_dev = (r2_t){{ 1.0, 1.0 }};  /* Should be a parameter... */

  /*Output data directory*/
  mkdir(out_dir, 0777);

  /* Read a set of camera calibration data: */
  tf_calib_data_t *cdat = read_calibration_data(in_dir, "p");
  fprintf(stderr, "calibration data:\n");
  tf_calib_data_print(stderr, cdat); 
  fprintf(stderr, "-----------------------------\n");

  int nmarks = cdat->np;
  r2_t *q = cdat->image;
  r3_t *p = cdat->world;
  double *w = cdat->weight;

  /* Read the true {cpar}, if present: */
  tf_camera_params_t *true_cpar = read_camera_parameters(cspec, in_dir, "true");
  tf_write_cpar_and_errors("true", cdat, cspec, q_dev, use_dev, true_cpar, out_dir); 
  
  /* Read the guess {cpar}: */
  tf_camera_params_t *guess_cpar = read_camera_parameters(cspec, in_dir, "guess");
  tf_write_cpar_and_errors("guess", cdat, cspec, q_dev, use_dev, guess_cpar, out_dir); 

  bool_t debug = TRUE;

  tf_camera_params_t *opt1_cpar = tf_camera_params_copy(guess_cpar);

  /* First optimization pass: */
  tf_calib_generic
    ( cdat, cspec, q_dev, use_dev, opt1_cpar,
      tf_calib_refine2_gather_optimization_params,
      tf_calib_refine2_scatter_optimization_params,
      debug ); 
  tf_write_cpar_and_errors("opt1", cdat, cspec, q_dev, use_dev, opt1_cpar, out_dir); 

  /* Adjust weights based on position discrepancies: */
  /* Estimate the standard deviation {dev_bad} of outliers: */
  r2_t dev_bad;
  double cook1 = 5.0; /* Guessed ratio of outlier-dev/inlier-dev. */
  r2_scale(cook1, &q_dev, &dev_bad);
  /* Estimate the standard deviation {dev_gud} of inliers: */
  r2_t dev_gud; /* Standard deviation for inliers. */
  double cook2 = 0.10; /* Guessed relative effect of outliers on predicted {q}s. */
  r2_mix(1.0, &q_dev, cook2, &dev_bad, &dev_gud);
  tf_recompute_target_weights(opt1_cpar, nmarks, p, q, w, dev_gud, dev_bad, TRUE);

  /* Second optimization pass: */
  tf_camera_params_t *opt2_cpar = tf_camera_params_copy(opt1_cpar);
  tf_calib_generic
    ( cdat, cspec, q_dev, use_dev, opt2_cpar,
      tf_calib_refine2_gather_optimization_params,
      tf_calib_refine2_scatter_optimization_params,
      debug ); 
  tf_write_cpar_and_errors("opt2", cdat, cspec, q_dev, use_dev, opt2_cpar, out_dir); 

  plot_cpars(cdat, guess_cpar, opt1_cpar, opt2_cpar, true_cpar, out_dir);

  tf_calib_data_free(cdat); 

  tf_camera_params_free (opt1_cpar);
  tf_camera_params_free (opt2_cpar);
  tf_camera_params_free (true_cpar);
  tf_camera_params_free (guess_cpar);

  return 0;
}

void plot_cpars
  ( tf_calib_data_t * cdat,
    tf_camera_params_t *guess_cpar,
    tf_camera_params_t *opt1_cpar,
    tf_camera_params_t *opt2_cpar,
    tf_camera_params_t *true_cpar,
    char *out_dir
  )
  {
    /* Camera image dimensions and aspect ratio: */
    double Nx = guess_cpar->Npx;
    double Ny = guess_cpar->Npy;
    double aspect = Nx/Ny;

    /* Postscript figure sizeand margin, in pt: */
    double hSize = 400.0*sqrt(aspect);
    double vSize = 400.0/sqrt(aspect);
    double hMrg = 4;
    double vMrg = 3;

    /* Output Postscript file: */
    char *ps_prefix = NULL;
    asprintf(&ps_prefix, "%s/", out_dir);
    PSStream *ps = pswr_new_stream(ps_prefix, NULL, TRUE, NULL, NULL, FALSE, hSize + 2*hMrg, vSize + 2*vMrg);

    /* Plot the figures: */
    int nmarks = cdat->np;
    r2_t *q = cdat->image;
    r3_t *p = cdat->world;

    tf_camera_params_t *fig1[2] = { true_cpar, guess_cpar };
    tf_camera_params_t *fig2[2] = { true_cpar, opt1_cpar  };
    tf_camera_params_t *fig3[2] = { true_cpar, opt2_cpar  };

    tf_plot_cameras(ps, "true_guess",  nmarks, p, 2, fig1, q, hSize, vSize, hMrg, vMrg, Nx, Ny, 1.5, 0);
    tf_plot_cameras(ps, "true_opt1",   nmarks, p, 2, fig2, q, hSize, vSize, hMrg, vMrg, Nx, Ny, 1.5, 0);
    tf_plot_cameras(ps, "true_opt2",   nmarks, p, 2, fig3, q, hSize, vSize, hMrg, vMrg, Nx, Ny, 1.5, 0);

    pswr_close_stream(ps);
  }

tf_calib_data_t *read_calibration_data(char *in_dir, char *tag)
  {  char *world_coords_fname = NULL;
    asprintf(&world_coords_fname, "%s/%s_w.txt", in_dir, tag);
    char *image_coords_fname = NULL;
    asprintf(&image_coords_fname, "%s/%s_i.txt", in_dir, tag);
    char *point_weights_fname = NULL;
    asprintf(&point_weights_fname, "%s/%s_wgt.txt", in_dir, tag);
    tf_calib_data_t * cdat = tf_calib_data_read
      ( world_coords_fname, image_coords_fname, point_weights_fname );
    free(image_coords_fname);
    free(world_coords_fname);
    free(point_weights_fname);
    return cdat;
  }
 
tf_camera_params_t *read_camera_parameters(tf_camera_specs_t *cspec, char *in_dir, char *name)
  { char *cpar_fname = NULL;
    asprintf(&cpar_fname, "%s/%s.cpar", in_dir, name);
    FILE *cpar_file = open_read(cpar_fname, TRUE);
    tf_camera_params_t *cpar = tf_camera_specs_get_new_mean_params(cspec);
    tf_camera_params_read(cpar_file, cpar);
    if (cpar_file != stdin) { fclose(cpar_file); }
    free(cpar_fname);
    return cpar;
  }
