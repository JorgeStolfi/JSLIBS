/* Last edited on 2022-10-20 05:57:11 by stolfi */

#define PROG_NAME "test_calib_guess"
#define PROG_DESC "tests the initial camera calibration guess algorithms"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include <float_image.h>
#include <r2x2.h>
#include <r2.h>
#include <pswr.h>
#include <jsfile.h>
#include <tf_calib_refine2.h>
#include <tf_camera.h>
#include <tf_camera_plot.h>
#include <tf_errors.h>
#include <tf_targets.h>
#include <tf_calib.h>
#include <tf_calib_guess1.h>
#include <tf_calib_guess2.h>
#include <sys/stat.h>
#include <sys/types.h>

void plot_cpars
  ( tf_camera_specs_t *cspec,
    tf_calib_data_t * cdat,
    tf_camera_params_t *affine_cpar,
    tf_camera_params_t *guess_cpar,
    tf_camera_params_t *true_cpar,
    char *out_dir
  );
  
tf_calib_data_t *read_calibration_data(char *in_dir, char *tag);

tf_camera_params_t *read_camera_parameters(tf_camera_specs_t *cspec, char *in_dir, char *name);

int32_t main (int32_t argc, char *argv[])
{
  char *camera_name = argv[1];
  char *in_dir = argv[2]; /* Directory where to find {q.txt,pgt.txt,true.cpar}. */
  char *out_dir = argv[3]; /* Directory for all output files. */
  int32_t which_calib_guess = atoi(argv[4]); /* Version of {tf_calib_guessN} to use. */

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

  /* A priori uncertainty on image point positions (for error computation): */
  r2_t q_dev = (r2_t){{ 5.0, 5.0 }};  /* Should be a parameter... */

  /*Output data directory*/
  mkdir(out_dir, 0777);

  /* Read the calibration data: */
  tf_calib_data_t *cdat = read_calibration_data(in_dir, "p");
  fprintf(stderr, "calibration data:\n");
  tf_calib_data_print(stderr, cdat); 
  fprintf(stderr, "-----------------------------\n");

  /* Read the true {cpar}, if present: */
  tf_camera_params_t *true_cpar = read_camera_parameters(cspec, in_dir, "true");
  tf_write_cpar_and_errors("true", cdat, cspec, q_dev, FALSE, true_cpar, out_dir); 

  /* Compute the affine and guessed camera parameters: */
  tf_camera_params_t *affine_cpar = NULL;
  tf_camera_params_t *guess_cpar = NULL;
  if (which_calib_guess == 1)
    { guess_cpar = tf_camera_specs_get_new_mean_params(cspec);
      tf_calib_guess1_initial_camera_parameters(cdat, cspec, q_dev, guess_cpar, TRUE);
    }
  else if (which_calib_guess == 2)
    { affine_cpar = tf_camera_specs_get_new_mean_params(cspec);
      tf_calib_guess2_initial_camera_parameters(cdat, cspec, affine_cpar, TRUE, TRUE );
      guess_cpar = tf_camera_specs_get_new_mean_params(cspec);
      tf_calib_guess2_initial_camera_parameters(cdat, cspec, guess_cpar, FALSE, TRUE );
    }
  tf_write_cpar_and_errors("affine", cdat, cspec, q_dev, FALSE, affine_cpar, out_dir); 
  tf_write_cpar_and_errors("guess",  cdat, cspec, q_dev, FALSE, guess_cpar,  out_dir); 

  plot_cpars(cspec, cdat, affine_cpar, guess_cpar, true_cpar, out_dir);

  tf_calib_data_free(cdat); 

  tf_camera_params_free (guess_cpar);
  tf_camera_params_free (true_cpar);

  return 0;
}

void plot_cpars
  ( tf_camera_specs_t *cspec,
    tf_calib_data_t * cdat,
    tf_camera_params_t *affine_cpar,
    tf_camera_params_t *guess_cpar,
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

    /* Open the output Encapsulated Postscript file: */
    char *ps_prefix = NULL;
    asprintf(&ps_prefix, "%s/", out_dir);
    PSStream *ps = pswr_new_stream(ps_prefix, NULL, TRUE, NULL, NULL, FALSE, hSize + 2*hMrg, vSize + 2*vMrg);

    /* Plot the figures: */
    int32_t np = cdat->np;
    r2_t *q = cdat->image;
    r3_t *p = cdat->world;

    tf_camera_params_t *fig1[2] = { true_cpar, affine_cpar };
    tf_camera_params_t *fig2[2] = { true_cpar, guess_cpar };

    tf_plot_cameras(ps, "true_affine", np, p, 2, fig1, q, hSize, vSize, hMrg, vMrg, Nx, Ny, 1.0, 0);
    tf_plot_cameras(ps, "true_guess",  np, p, 2, fig2, q, hSize, vSize, hMrg, vMrg, Nx, Ny, 1.0, 0);

    /* Finish the postscript file: */
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
